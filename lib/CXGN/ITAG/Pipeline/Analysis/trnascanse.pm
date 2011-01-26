package CXGN::ITAG::Pipeline::Analysis::trnascanse;
use strict;
use warnings;
use English;
use Carp;

use File::Temp qw/tempdir/;
use File::Copy;

use Bio::SeqIO;
use Bio::Tools::tRNAscanSE;
use Bio::FeatureIO;

use CXGN::Tools::Class qw/parricide/;

=head1 NAME

CXGN::ITAG::Pipeline::Analysis::trnascanse - pipeline analysis for tRNAScanSE

=head1 BASE CLASS(ES)

L<CXGN::ITAG::Pipeline::Analysis>

=cut

use base qw/CXGN::ITAG::Pipeline::Analysis/;

=head1 METHODS

implements all abstract methods defined in L<CXGN::ITAG::Pipeline::Analysis>

=cut

sub locally_runnable {1}

=head2 run

  Usage: $an->run
  Desc : run this analysis, overrides the run() defined in the base
         class.
  Args : batch object
  Ret  : nothing meaningful
  Side Effects: runs this analysis

=cut

sub run {
  my ($self,$batch,%args) = @_;

  # mapping sequence names -> pairs of tempfiles and destinations to be fed to move
  # we assemble all these files first so the move will only happen if they all come out right
  my @ops = map {
    my $seqname = $_;

    #get the fasta file for this sequence name
    my ($seq_file) = $batch->pipeline->analysis('seq')->files_for_seq($batch,$seqname);
    $seq_file && -f $seq_file
      or die "expected sequence file '$seq_file' not found";

    my $temp_gff = $self->cluster_temp($seqname,"gff3");
    my $temp_raw = $self->cluster_temp($seqname,"raw");

    my ($raw_out,$gff3_out) = $self->files_for_seq($batch,$seqname);

    #run trnascan on it
    my $tse = CXGN::Tools::Run->run_cluster
        ( 'itag_wrapper', 'tRNAscan-SE', #<assume repeatmasker is in path
          $seq_file,
          { out_file => $temp_raw,
            working_dir => $self->cluster_temp($seqname),
            temp_base   => $self->cluster_temp($seqname),
            on_completion => sub {
                #convert the raw out to gff3
                stat $temp_raw; #< help nfs uncaching
                my $fi = Bio::Tools::tRNAscanSE->new(-file => $temp_raw);
                my $fo = Bio::FeatureIO->new(-file => ">$temp_gff", -format => 'gff', -version => 3);
                while (my $feature = $fi->next_prediction() ) {
                    $feature->primary_tag('tRNA');
                    my $f = Bio::SeqFeature::Annotated->new( -feature => $feature );
                    $fo->write_feature( $f );
                }
            },
          }
        );


    [$temp_raw => $raw_out,  $tse],
    [$temp_gff => $gff3_out, $tse]
  } $batch->seqlist;

  # wait for all the cluster jobs to finish
  sleep 2 while grep $_->[2]->alive, @ops;

  # move all the output files into position,
  # rolling back (unlinking targets) on failure
  $self->atomic_move( @ops );
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
