#!/usr/bin/env perl
use strict;
use warnings;
use autodie ':all';

use Cwd;

use Storable qw/ nstore /;
use Data::Dumper;
use List::Util qw/ shuffle max /;

use File::Path qw/ make_path remove_tree /;

use Path::Class;

use Bio::SeqIO;

use CXGN::Tools::List qw/ balanced_split /;
use CXGN::ITAG::Pipeline;

my $cdna_dir    = 'input/cdna';
my $genomic_dir = 'input/genomic';

my @cdna_files    = generate_cdna_input(    $cdna_dir    );
my @genomic_files = generate_genomic_input( 'genomic.fa', $genomic_dir );

# and now we generate an all-versus-all pairing of genomic and cdna files
my @file_pairs = map {
    my $g = $_;
    map [ $_, $g ], @cdna_files;
} @genomic_files;


my @jobs;
my $job_number = 0;

foreach my $pairs ( balanced_split( 200, @file_pairs ) ) {
    # each $pairs is an arrayref like [ [c,g],[c,g],[c,g], ... ]

    my $jobdir  = dir( 'jobs', ++$job_number );
    $jobdir->mkpath;

    print "job $job_number - $jobdir - ".@$pairs." pairs\n";

    if( -f  $jobdir->file('done') ) {
	print "done, skipping.\n";
	next;
    }

    my %job_record = ( pairs   => $pairs,
		       dir     => $jobdir->stringify,
		       script => <<'EOS',
		       use CXGN::Tools::Run;

		       my $out_cnt = 0;
		       foreach (@{$j->{pairs}}) {
                         my ($cdna, $genomic) = @$_;
			 CXGN::Tools::Run->run( 'gth',
						'-xmlout',
						-minalignmentscore => '0.90',
						-mincoverage       => '0.90',
						-seedlength        => 16,
						-o                 => $j->{dir}.'/'.++$out_cnt.'.xml',
						-species           => 'arabidopsis',
						-cdna              => $cdna,
						-genomic           => $genomic,
			                      );
		       }
		       unlink $g;
		       open my $f, '>', $j->{dir}.'/done';
EOS
	             );

    my $vmem_est = max map estimate_vmem( $_->[0], $_->[1] ), @$pairs;

    my $jobfile = "$jobdir/job.dat";
    nstore \%job_record, $jobfile;

    my $job = CXGN::Tools::Run->run_cluster( 'itag_wrapper',
					     'perl',
					     '-MStorable=retrieve',
					     -E => q|my $j = retrieve($ARGV[0]) or die; #eval $j->{script}; die $@ if $@|,
					     $jobfile,
					     {
						 temp_base => $jobdir,
						 vmem => $vmem_est,
					     },
	                                   );
    push @jobs, $job;

    $_->alive for @jobs;
}

sleep 10 while grep $_->alive, @jobs;

exit;

######### SUBS
# estimates the resources for each sequence name, used for giving resource hints to the job scheduler

sub estimate_vmem {
    my ( $cdna_file, $genomic_file ) = @_;

    my $cdna_size = -s $cdna_file
	or die "sanity check failed, no size found for cdna file $cdna_file";
    my $genomic_size = -s $genomic_file
	or die "sanity check failed, no size found for genomic file $genomic_file";


    $_ /= 1_000_000 for $cdna_size, $genomic_size;

    my $vmem_est = sprintf('%0.0f',
			       6 * ( $genomic_size + $cdna_size )
			   +   3 *   $genomic_size * $cdna_size
			  );
    #print "cdna size: $cdna_size, seq size: $seq_size => vmem $vmem_est M\n";

    return $vmem_est;
}

sub generate_genomic_input {
    my $source_file = file( shift );
    my $input_dir   = dir( shift );
    $input_dir->mkpath;

    my $src;
    if( $source_file =~ /\.gz$/ ) {
	open $src, "gunzip -c $source_file |";
    } else {
	open $src, '<', $source_file;
    }
    my $seqs = Bio::SeqIO->new( -fh => $src, -format => 'fasta' );
    while( my $seq = $seqs->next_seq ) {
	Bio::SeqIO->new( -format => 'fasta',
			 -file => '>'.$input_dir->file( $seq->id.'.fa')->stringify,
		       )
	          ->write_seq( $seq );
    }

    return grep /\.fa$/, $input_dir->children;
}

sub generate_cdna_input {
    my $input_dir = dir( shift );

    $input_dir->rmtree;
    $input_dir->mkpath;

    my $pipe = CXGN::ITAG::Pipeline->open( version => 1 );
    my $batch = $pipe->batch( 3 );
    my $an = $pipe->analysis( 'renaming' );

    my @files =
        # split each cdna file into many cdna files, one per cdna seq
        map {
	    my $cdna_file = file( $_ );
	    my $cdna_seqs = Bio::SeqIO->new( -format => 'fasta',
					     -file => "$cdna_file",
					    );

	    my @files;
	    my $cdna_count;
	    my $input_subdir = $input_dir->subdir($cdna_file->basename );
	    $input_subdir->mkpath;
	    while ( my $single_cdna_seq = $cdna_seqs->next_seq ) {
		my $single_cdna_file = $input_subdir->file( ++$cdna_count.'.fa' );
		Bio::SeqIO->new( -format => 'fasta',
				 -file   => ">$single_cdna_file",
				)
    		          ->write_seq( $single_cdna_seq );
		push @files, $single_cdna_file;
	    }

	    @files
        }
        # skip empty cdna files
        grep -s,
        # get all the cdna files
        map {
	    my (undef,undef,undef,$cdna_file) = $an->files_for_seq( $batch, $_ );
	    $cdna_file
        }
        # starting with a sorted sequence list
        sort $batch->seqlist;

    return @files;
}
