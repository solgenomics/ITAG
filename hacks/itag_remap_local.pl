#!/usr/bin/env perl
use strict;
use warnings;
use autodie ':all';

use Cwd;
use Data::Dumper;
use List::Util qw/ shuffle max /;
use File::Temp;
use Storable qw/ nstore /;


use CXGN::Tools::List qw/ balanced_split /;
use CXGN::Tools::Wget qw/ wget_filter /;
use CXGN::ITAG::Pipeline;

#######

#my $genomic = '/data/prod/public/tomato_genome/wgs/assembly/S_lycopersicum_scaffolds_20100122.fa.gz';
my $genomic = 'genomic.fa';

#######

my $pipe = CXGN::ITAG::Pipeline->open( version => 1 );
my $batch = $pipe->batch( 3 );
my $an = $pipe->analysis( 'renaming' );

my $job_number = 0;

my $unzipped_genomic = File::Temp->new;
{ open my $unzipped, "gunzip -c $genomic |";
  $unzipped_genomic->print($_) while <$unzipped>;
}
$unzipped_genomic->seek(0);
$unzipped_genomic->close;

 scalar <$unzipped_genomic>;


# get a list of all cDNA files that are nonzero size
my @cdna_files =
    grep -s,
    map {
	my (undef,undef,undef,$cdna_file) = $an->files_for_seq( $batch, $_ );
	$cdna_file
    }
    sort $batch->seqlist;

my $cdna_size;
my @jobs;
my @running_jobs;
foreach my $cdna_job ( balanced_split( 200, @cdna_files ) ) {
    # find the cdna result files for the seqs in this job

    my $jobdir  = getcwd()."/".++$job_number;
    mkdir $jobdir;
    print "job $job_number - $jobdir - ".@$cdna_job." seqs\n";
    if( -f "$jobdir/done" ) {
	print "done, skipping.\n";
	next;
    }
    my %job_record = ( cdna    => $cdna_job,
		       genomic => $genomic,
		       dir     => $jobdir,
		       script => <<'EOS',
		       use CXGN::Tools::Wget 'wget_filter';
		       use CXGN::Tools::Run;
		       my $g = wget_filter( $j->{genomic} => "$j->{dir}/genomic.fa", { gunzip => 1 } );

		       my $out_cnt = 0;
		       foreach my $c (@{$j->{cdna}}) {
			 CXGN::Tools::Run->run( 'gth',
						'-xmlout',
						-minalignmentscore => '0.90',
						-mincoverage       => '0.90',
						-seedlength        => 16,
						-o                 => $j->{dir}.'/'.++$out_cnt.'.xml',
						-species           => 'arabidopsis',
						-cdna              => $c,
						-genomic           => $g,
			                      );
		       }
		       unlink $g;
		       open my $f, '>', $j->{dir}.'/done';
EOS
	             );

    $cdna_size ||= -s $genomic
	or die "no cdna size";

    my $vmem_est = max map estimate_vmem($cdna_size,$_), @$cdna_job;

    my $jobfile = "$jobdir/job.dat";
    nstore \%job_record, $jobfile;

    my $job = CXGN::Tools::Run->run_cluster( 'itag_wrapper',
					     'perl',
					     '-MStorable=retrieve',
					     -E => q|my $j = retrieve($ARGV[0]) or die; eval $j->{script}; die $@ if $@|,
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

sub wget_size {
    my $url = shift;
    my $file = wget_filter($url);
    my $s = -s $file;
    unlink $file;
    return $s;
}

sub estimate_vmem {
    my ( $cdna_size, $seq_file ) = @_;

    my $seq_size = -s $seq_file
	or die "sanity check failed, no size found for seq file $seq_file";

    $_ /= 1_000_000 for $seq_size, $cdna_size;

    my $vmem_est = sprintf('%0.0f',
			       6 * ( $seq_size + $cdna_size )
			   +   3 *   $seq_size * $cdna_size
			  );
    #print "cdna size: $cdna_size, seq size: $seq_size => vmem $vmem_est M\n";

    return $vmem_est;
}
