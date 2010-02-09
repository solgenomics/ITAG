#!/usr/bin/env perl
use strict;
use warnings;
use autodie ':all';

use Cwd;

use Storable qw/ nstore retrieve /;
use Data::Dumper;
use List::Util qw/ shuffle max /;

use Getopt::Std;

use Set::CrossProduct;

use File::Path qw/ make_path remove_tree /;

use Path::Class;

use Bio::SeqIO;

use CXGN::Tools::List qw/ balanced_split balanced_split_sizes /;
use CXGN::ITAG::Pipeline;
use CXGN::ITAG::SeqSource::Fasta;

# auto-flush our output for more timely status messages
$| = 1;

our %opt;
getopts('g:c:Vj:',\%opt) or die 'invalid usage';
###########

my $genomic_chunks = $opt{g} || 300;
my $cdna_chunks    = $opt{c} || 30;
my $job_chunks     = $opt{j} || 1000;
my $batch_number   = 5;

###########

if( $opt{V} ) {
    die 2**20 * estimate_vmem(@ARGV);
}

# check the number of files expected in job output dirs under these setting
my $expected_job_outfiles =
    3 #control files
    + int($genomic_chunks * $cdna_chunks / $job_chunks ); #.xml output files
my $tempdir_files = int((13 * $genomic_chunks + 8 * $cdna_chunks) / $job_chunks); # tempdir files
if( 31999 < $expected_job_outfiles ) {
    die "invalid settings, these would make each job dir have $expected_job_outfiles files in it.\n";
}
if( 31999 < $tempdir_files ) {
    die "invalid settings, these would make each job temp dir have $tempdir_files files in it.\n";
}

$SIG{__DIE__} = \&Carp::confess;

my $cdna_dir    = dir( 'input/cdna' );
my $genomic_dir = dir( 'input/genomic' );
my $jobs_dir    = dir( 'jobs' );

my $cdna_files    = generate_cdna_input   ( 'cdna.fa',  $cdna_chunks, $cdna_dir, $jobs_dir );
my $genomic_files = generate_genomic_input( 'genomic.fa', $genomic_chunks, $genomic_dir, $jobs_dir );

# check all the vmem estimates for our job sizes 
print "checking memory requirement estimates...\n";
{ my $max_vmem = 5000;
  my $max_g = largest_file( @$genomic_files );
  my $max_c = largest_file( @$cdna_files    );
  my $vmem_est = estimate_vmem( $max_c, $max_g );
  $vmem_est > $max_vmem
      and die "vmem estimate $vmem_est for file $max_c vs $max_g is more than max $max_vmem. aborting.\n";
}

my $file_pairs = Set::CrossProduct->new( [ $genomic_files, $cdna_files ] );
my $job_sizes = balanced_split_sizes( $job_chunks, $file_pairs->cardinality );
my $numjobs = @$job_sizes;
my @jobs;
my $job_number = 0;
foreach my $jobsize ( @$job_sizes ) {
    my $jobdir  = $jobs_dir->subdir( ++$job_number );

    print "job $job_number/$numjobs - $jobsize pairs\n";

    my @pairs;    my $jobsize_copy = $jobsize;
    push @pairs, scalar $file_pairs->get
	while $jobsize_copy--; #< minimizing memory footprint

    my $donefile  = $jobdir->file('done');
    if( -f $donefile  ) {
	print "done, skipping.\n";
	next;
    }

    my $save_file = $jobdir->file('job.save');  #< file for saving the job record, so we can check if it's running
    my $job;
    if( -f $save_file ) {
	$job = eval {retrieve( $save_file )};
	warn $@ if $@;
    }

    if( $job && eval {$job->alive} ) {
	print "using existing, running job ".$job->job_id."\n";
    } else {
	$job = submit_job( $jobdir, $donefile, \@pairs );
	nstore( $job, $save_file );
    }

    push @jobs, $job;

    for my $j ( @jobs ) {
	eval { $j->alive };
	warn $@ if $@;
    }
}


sleep 10 while grep $_->alive, @jobs;

$_->cleanup for @jobs;

exit;

######### SUBS

# estimates the resources for each sequence name, used for giving
# resource hints to the resource manager (which, currently, is Torque)
sub estimate_vmem {
    my ( $cdna_file, $genomic_file ) = @_;

    my $cdna_size = _cached_file_size( $cdna_file )
	or die "sanity check failed, no size found for cdna file $cdna_file";
    my $genomic_size = _cached_file_size( $genomic_file )
	or die "sanity check failed, no size found for genomic file $genomic_file";

    my $megs = 2**20;
    my $vmem_est = sprintf('%0.0f',
			   (   200*$megs
			     + 100*($genomic_size + $cdna_size)
			     + '3e-03' * $genomic_size * $cdna_size
                           )/$megs
			  );
    #print "cdna $cdna_size, genomic size: $genomic_size => vmem $vmem_est M\n";

    return $vmem_est;
}

sub generate_genomic_input {
    my $source_file = file( shift );
    my $chunk_count = shift;
    my $input_dir   = dir( shift );

    my $donefile = $input_dir->file('done');
    unless( -f $donefile ) {
	print "generating genomic input files...";
	$input_dir->rmtree;
	$input_dir->mkpath;

	chunk_fasta_file( $source_file, $chunk_count, $input_dir );

	$donefile->openw->print("done ".time."\n");
	print "done.\n";
    } else {
	print "using cached genomic input.\n";
    }
    my @files;
    $input_dir->recurse( callback => sub { $_ = shift; _cached_file_size($_); push @files, $_ if -f && /\.fa$/ } );
    return \@files;
}

memoize('_cached_file_size');
sub _cached_file_size {
    -s $_[0]
}

sub chunk_fasta_file {
    my ( $file, $chunk_count, $dir, $make_name ) = @_;
    $dir = dir( $dir );
    $make_name ||= sub { my ($d,$chunk_num,$seq_list) = @_;
			 return $d->file( $chunk_num.'.fa')->stringify
		       };

    my $src = CXGN::ITAG::SeqSource::Fasta->new( file => $file );
    my $chunk_idx;
    foreach my $chunk ( @{ balanced_split( $chunk_count, [ shuffle $src->all_seq_names ] ) } ) {
	my $out = Bio::SeqIO->new( -format => 'fasta',
				   -file => '>'.$make_name->( $dir, ++$chunk_idx, $chunk ),
				  );
	foreach (@$chunk) {
	    my ($seq,undef) = $src->get_seq( $_ );
	    $out->write_seq( $seq );
	}
    }
}

sub generate_cdna_input {
    my $all_seqs = shift;
    my $chunk_count = shift;
    my $input_dir = dir( shift );

    my $donefile = $input_dir->file('done');
    if( -f $donefile ) {
	print "using cached cdna input.\n";
	return _cdna_file_list( $input_dir );
    }

    print "generating cdna input...";

    system "rm -rf $input_dir";
    $input_dir->mkpath;

    chunk_fasta_file( $all_seqs, $chunk_count, $input_dir );

    $donefile->openw->print("done\n");
    print "done.\n";

    return _cdna_file_list( $input_dir );
}

sub _cdna_file_list {
    my $input_dir = shift;
    my @files;
    $input_dir->recurse( callback => sub { $_ = shift; _cached_file_size($_); push @files, $_ if -f && /\.fa$/ } );
    return \@files;
}


sub submit_job {
    my ( $jobdir, $donefile, $pairs ) = @_;
    system "rm -rf $jobdir";
    $jobdir->mkpath;

    my $pairs_file = $jobdir->file('pairs');
    { local $Data::Dumper::Indent = 0;
      local $Data::Dumper::Terse  = 1;
      my $pairs_fh  = $pairs_file->openw;
      $pairs_fh->say( Dumper [ "$_->[0]", "$_->[1]" ] ) for @$pairs;
    }

    my %job_record = ( pairs    => $pairs_file,
		       dir      => $jobdir->stringify,
		       donefile => $donefile->stringify,
		       script => <<'EOS',
		       use Cwd 'abs_path';
		       use autodie ':all';
		       use File::Temp;
		       use File::Basename;

		       my $tempd;
		       my $tempd_count = 0;
		       my $pair_number = 0;
                       open my $pairs, $j->{pairs};
                       while( my $p = <$pairs> ) {
                         $p = eval $p;
                         my ($genomic,$cdna) = @$p;

			 # make a new tempdir every 300 runs of gth
			 if( !$tempd || $tempd_count++ > 300 ) {
			     $tempd = File::Temp->newdir( CLEANUP => 1 );
			     $tempd_count = 0;
			 }

			 # make tempfiles
			 my $cdna_s = "$tempd/c_".basename($cdna);
			 my $genomic_s = "$tempd/g_".basename($genomic);
			 symlink abs_path($cdna), $cdna_s unless -e $cdna_s;
			 symlink abs_path($genomic), $genomic_s unless -e $genomic_s;

                         my $outfile = $j->{dir}.'/'.++$out_cnt.'.xml';
                         #system "echo $cdna, $genomic > $outfile"; next;
			 system( 'gth',
						'-xmlout',
						'-force', #overwrite output
						-minalignmentscore => '0.97',
						-mincoverage       => '0.90',
						-seedlength        => 16,
						-o                 => $outfile,
						-species           => 'arabidopsis',
						-cdna              => $cdna_s,
						-genomic           => $genomic_s,
			                      );
		       }
		       open my $f, '>', $j->{donefile};
                       print $f "done\n";
EOS
	             );

    my $vmem_est = max map estimate_vmem( $_->[0], $_->[1] ), @$pairs;

    my $jobfile = $jobdir->file('job.dat');
    nstore \%job_record, $jobfile;

    my $job = CXGN::Tools::Run->run_cluster( 'itag_wrapper',
					     'perl',
					     '-MStorable=retrieve',
					     -E => q|my $j = retrieve($ARGV[0]); eval $j->{script}; die $@ if $@|,
					     $jobfile,
					     {
						 temp_base => $jobdir,
						 vmem      => $vmem_est,
					     },
	                                   );

    print "submitted as ".$job->job_id.", est. vmem $vmem_est\n";

    return $job;
}

# finds the largest of a list of files
sub largest_file {
    my $s = 0;
    my $f;
    for (@_) {
	my $ts = _cached_file_size($_);
	if( $ts > $s ) {
	    $s = $ts;
	    $f = $_;
	}
    }
    return $f;
}
