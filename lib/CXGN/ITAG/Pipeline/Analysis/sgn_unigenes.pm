package CXGN::ITAG::Pipeline::Analysis::sgn_unigenes;
use strict;
use warnings;
use English;
use Carp;

use Bio::SeqIO;

use CXGN::Tools::Wget qw/ wget_filter /;

=head1 NAME

CXGN::ITAG::Pipeline::Analysis::sgn_unigenes - pipeline analysis to produce FILL IN SOMETHING HERE

=head1 BASE CLASS(ES)

L<CXGN::ITAG::Pipeline::Analysis::gth_base>

=cut

use base qw/CXGN::ITAG::Pipeline::Analysis::gth_base/;

=head1 METHODS

implements all abstract methods defined in L<CXGN::ITAG::Pipeline::Analysis>,
also overrides FILL IN OVERRIDES HERE

=cut

#this analysis is locally runnable

sub locally_runnable { 1 }

sub _source {
    'SGN_Unigenes'
}

sub _cdna_file {
    my $class = shift;
    wget_filter( 'cxgn-resource://sgn_unigenes_combined'
                 => $class->local_temp('sgn_unigenes_combined.fasta')
               );
}

sub _parse_mode {
    'alignments_merged'
}

=head1 AUTHOR(S)

Robert Buels

=cut

###
1;#do not remove
###
