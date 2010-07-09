package CXGN::ITAG::Pipeline::Analysis::augustus;
use strict;
use warnings;
use English;
use Carp;

use File::Temp qw/tempdir/;
use File::Copy;

use CXGN::Tools::Class qw/parricide/;

=head1 NAME

CXGN::ITAG::Pipeline::Analysis::augustus - pipeline analysis to produce AUGUSTUS gene predicitions

=head1 BASE CLASS(ES)

L<CXGN::ITAG::Pipeline::Analysis>

=cut

use base qw/CXGN::ITAG::Pipeline::Analysis/;

=head1 METHODS

implements all abstract methods defined in L<CXGN::ITAG::Pipeline::Analysis>,

=cut

#this analysis is locally runnable
sub locally_runnable { 1 }

#and this is the routine to run it
sub run {
    my ($self,$batch) = @_;

    # mapping sequence names -> pairs of tempfiles and destinations to be fed to move
    # we assemble all these files first so the move will only happen if they all come out right
    my @ops = map {
        my $seqname = $_;

        #get the fasta file for this sequence name
        my ($seq_file) = $batch->pipeline->analysis('repeats')->files_for_seq($batch,$seqname);
        $seq_file && -f $seq_file
            or die "expected sequence file '$seq_file' not found";

        my $temp_seq = $self->cluster_temp( $seqname, "seq" )
            or confess "could not create seq temp file\n";

        #preprocess the seq file to remove the defline, since augustus puts it in the output otherwise
        open my $seq_in, $seq_file or die "$! opening '$seq_file'";
        open my $temp_seq_out,">$temp_seq" or die "$! opening $temp_seq for writing\n";
        while (<$seq_in>) {
            s/(>\S+).+/$1/ if />/;
            print $temp_seq_out $_;
        }
        close $temp_seq_out;
        close $seq_in;

        my $temp_gff   = $self->cluster_temp( $seqname, "gff3_1" );
        my $temp_gff_2 = $self->local_temp(   $seqname, "gff3_2" );
        #my $temp_gff_3 = "$tempdir/gff3_3";

        my ($gff3_out) = $self->files_for_seq($batch,$seqname);

        my $tse = CXGN::Tools::Run->run_cluster
            ( 'itag_wrapper', 'augustus',       #<assume augustus is in path
              #augustus --species=human --hintsfile= hints.E.gff --extrinsicCfgFile=extrinsic.ME.cfg genome.fa 
              '--species=tomato2',
              '--extrinsicCfgFile=extrinsic.ME.cfg',
              '--introns=off',
              '--gff3=on',
              $temp_seq,
              { out_file    => $temp_gff,
                temp_base   => $self->cluster_temp($seqname),
                working_dir => $self->cluster_temp($seqname),
                on_completion => sub {
                    #convert the dirty gff3 to valid gff3
                    my $wrong_seqname = $seqname;
                    $wrong_seqname =~ s/.$//;
                    CXGN::Tools::Run->run('itag_wrapper', 'gff3_reformat.pl',
                                          -A => 'ID|Parent=s/^g(\d+)\b/'.$seqname.'.g$1/',
                                          '-I',
                                          -L => "s/$wrong_seqname(\\s)/$seqname\$1/; \$_= \"###\\n\" if /^###/",
                                          $temp_gff,
                                          {
                                           out_file => $temp_gff_2 }
                                         );
                },
              }
            );

        [$temp_gff_2 => $gff3_out, $tse]
    } $batch->seqlist;

    # wait for all the cluster jobs to finish
    sleep 2 while grep $_->[2]->alive, @ops;

    # move all of the output files into position, rolling back on
    # failure
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
