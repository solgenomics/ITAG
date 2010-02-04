package CXGN::ITAG::SeqSourceI;
use strict;
use warnings;
use English;
use Carp;

use base qw/ CXGN::DB::Ima /;

=head1 NAME

CXGN::ITAG::SeqSourceI - Interface for an object representing a source
of sequences for the ITAG annotation pipeline

=head1 SYNOPSIS

  package myseqsource;
  use base qw/CXGN::ITAG::SeqSourceI/;
  sub all_seq_names {
    my ($self) = @_;
    return qw/AC1234.1 AC12345.1/;
  }
  sub get_seq {
    my ($self,$seqname) = @_;
    return _lookup_name_and_get_bio_seq($seqname);
  }

=head1 DESCRIPTION

Interface specification for an object that can serve as a source of
sequences for the ITAG distributed annotation pipeline.

=head1 PROVIDED METHODS

=head2 new

  Usage: my $src = MySeqSource->new(%args);
  Desc : make a new seqsource object
  Args : hash-style list of options, passed through
         from CXGN::ITAG::Pipeline::create_batch()
  Ret  : a new seqsource object
  Side Effects: none

  The default implementation here just copies the args hash
  into $self, blesses that, and returns it.
  Override if you want something different.

=cut

sub new {
  my ($class,%args) = @_;

  my $self = bless \%args, $class;

  return $self;
}

=head1 ABSTRACT METHODS

=head2 all_seq_names

  Usage: my @seqs = $src->all_seq_names;
  Desc : get a list of all sequence names available from this source.
  Args : none
  Ret  : list of strings, which are seq names
  Side Effects: none

=cut

sub all_seq_names {
  confess 'all_seq_names not implemented';
}

=head2 get_seq

  Usage: my $seq = $src->get_seq($seqname);
  Desc : get a L<Bio::SeqI>-compliant object (such as a
         L<Bio::PrimarySeq>) for the given sequence name.
  Args : sequence name to fetch
  Ret  : undef if the name was not found, otherwise return
         a list of:
            sequence object (Bio::SeqI),
            a single contig arrayref like the ones
            returned by CXGN::BioTools::AGP::agp_contigs
              [ agp_hashref_line, ... ]

=cut

sub get_seq {
  my ($self,$seqname) = @_;
  confess 'get_seq not implemented, implement it in the child class';
}


=head1 AUTHOR(S)

Robert Buels

=cut

###
1;#do not remove
###
