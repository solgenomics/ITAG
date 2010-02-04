package CXGN::ITAG::SeqSource::TomatoAssembly;
use strict;
use warnings;
use English;
use Carp;

=head1 NAME

CXGN::ITAG::SeqSource::TomatoAssembly - ITAG pipeline sequence source
that provides tomato contig sequences from the new tomato assembly
from PRI Wageningen.

=head1 SYNOPSIS

my $src = CXGN::ITAG::SeqSource::TomatoAssembly->new( );

=head1 DESCRIPTION

=head1 BASE CLASS(ES)

L<CXGN::ITAG::SeqSourceI>

=cut

use base qw/CXGN::ITAG::SeqSourceI Class::Accessor/;

=head1 METHODS

Provides the methods specified in the SeqSourceI interface, using
tomato contig sequences.

=cut

#NOTE: Class::Accessor provides the new() method


use CXGN::DB::DBICFactory;

sub _chado {
    shift->{chado} ||= CXGN::DB::DBICFactory->open_schema('Bio::Chado::Schema');
}


sub all_seq_names {

    shift->_chado;

    # search for features of type contig with certain featureprops
    die 'todo';
}

sub get_seq {
  my ($self, $seqname) = @_;

  shift->_chado->resultset('...');
  die 'todo';
}

=head1 AUTHOR(S)

Robert Buels

=cut

###
1;#do not remove
###
