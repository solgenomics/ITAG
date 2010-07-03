package CXGN::ITAG::Pipeline::Analysis::repeats;
use strict;
use warnings;
use English;
use Carp qw/ cluck croak confess /;

use Memoize;

use File::Copy;
use File::Temp qw/ tempdir /;
use File::Spec;

use List::MoreUtils qw/ part /;

use Bio::SeqIO;

use CXGN::Tools::Class qw/parricide/;
use CXGN::Tools::Run;

use CXGN::Tools::Wget qw/wget_filter/;

=head1 NAME

CXGN::ITAG::Pipeline::Analysis::repeats - pipeline analysis to produce
repeat-masked sequences, masked using the SGN repeats master set.

=head1 BASE CLASS(ES)

L<CXGN::ITAG::Pipeline::Analysis>

=cut

use base qw/CXGN::ITAG::Pipeline::Analysis/;

=head1 METHODS

implements all abstract methods defined in L<CXGN::ITAG::Pipeline::Analysis>,

=cut

sub locally_runnable { 1 }

sub run {
  my ($self,$batch) = @_;

  my $stringent_lib = wget_filter('cxgn-resource://itag_repeats_stringent'
				  => $self->cluster_temp('itag_repeats_stringent.seq')
				 );
  my $regular_lib   = wget_filter('cxgn-resource://itag_repeats_normal',
				  => $self->cluster_temp('itag_repeats_normal.seq'),
				 );
  #warn "got repeat lib file $repeat_lib_file\n";

  my $job_idx = 0;
  my @ops = map {
    #make a soft-masked version of the repeat-masked file
    my $seqname = $_;

    my ($seq_file) = $batch->pipeline->analysis('seq')->files_for_seq($batch,$seqname);
    $seq_file && -f $seq_file
      or die "expected sequence file '$seq_file' not found";

    my ($hardmasked_fasta,$softmasked_fasta,$hardmasked_stringent,$softmasked_stringent) = $self->files_for_seq($batch,$_);

    #now run repeatmasker with two sets of options on the cluster
    ( $self->run_repeatmasker( $seqname, $regular_lib,   $seq_file, $hardmasked_fasta,     $softmasked_fasta,     'nolow', $job_idx++ ),
      $self->run_repeatmasker( $seqname, $stringent_lib, $seq_file, $hardmasked_stringent, $softmasked_stringent, 0      , $job_idx++ ),
    )
  } $batch->seqlist;

  # wait for all of the cluster jobs to finish
  sleep 2 while grep $_->[2]->alive, @ops;

  # move all of the files into position, rolling back on failure
  $self->atomic_move( @ops );
}

sub run_repeatmasker {
    my ( $self, $seqname, $lib_file, $seq_file, $hard_file, $soft_file, $nolow, $job_idx ) = @_;

    # make a separate dir for each cluster job
    my $dir = $self->cluster_temp($job_idx);
    mkdir $dir or die "$! doing mkdir $dir";

    my $temp_seq_file = File::Spec->catfile($dir,"myseq");
    symlink($seq_file,$temp_seq_file) and -l $temp_seq_file
        or die "could not symlink $seq_file -> $temp_seq_file : $!";
    my $temp_hmf = "$temp_seq_file.hmf";
    my $temp_masked_file = "$temp_seq_file.masked";

    #my $rm = CXGN::Tools::Run->run_async( qq|sleep 5; echo ">foo" > $temp_seq_file; echo ACTGACTAGCTAGCTAC >> $temp_seq_file|,

    #run repeatmasker to generate the softmasked fasta
    my $rm = CXGN::Tools::Run->run_cluster
        (
         'itag_wrapper', 'RepeatMasker',        #<assume repeatmasker is in path
         '-q',
         $nolow ? ('-nolow') : (),
         '-xsmall',
         -maxsize => 2000000,
         -lib      => $lib_file,
	 -parallel => 4, #< double-proc the repeatmasker jobs, so they run fast in the time that they are not i/o bound
         $temp_seq_file,
         {
	  temp_base   => $dir,
          working_dir => $dir,
          procs_per_node => 4,
          on_completion => sub {
              my ($job) = @_;
              #cluck("completed $temp_seq_file\n");
              #filter the softmasked fasta to get a hardmasked version
              if ( -f $temp_masked_file ) {
                  open my $smf,$temp_masked_file or confess "$seqname: $! opening $temp_masked_file: ".$job->err;
                  open my $hmf,">$temp_hmf" or  confess "$seqname: $! opening $temp_hmf for writing";
                  while (my $line = <$smf>) {
                      unless($line =~ /^\s*>/) {
                          $line =~ s/[a-z]/N/g;
                      }
                      print $hmf $line;
                  }
              } else {
                  #if we get here, it made no masked seq file.  check the error
                  #output of RepeatMasker for a 'no repeats' pattern, and just
                  #copy the original sequence files as the masked ones
                  copy( $temp_seq_file => $temp_masked_file );
                  copy( $temp_seq_file => $temp_hmf         );
              }
          },
         }
        );


    return ( [$temp_masked_file => $soft_file, $rm],
             [$temp_hmf         => $hard_file, $rm],
           )
}

=head1 AUTHOR(S)

Robert Buels

=cut

sub DESTROY {
  return parricide(shift,our @ISA);
}

###
1;#do not remove
###
