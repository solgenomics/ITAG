package CXGN::ITAG::SeqSource::Fake;
use strict;
use warnings;
use English;
use Carp;

use Bio::PrimarySeq;

use CXGN::DB::Connection;
use CXGN::Tools::Class qw/parricide/;

=head1 NAME

CXGN::ITAG::SeqSource::Fake - ITAG pipeline fake sequence source

=head1 SYNOPSIS

my $src = CXGN::ITAG::SeqSource::Fake->new()

=head1 DESCRIPTION


=head1 BASE CLASS(ES)

L<CXGN::ITAG::SeqSourceI>

=cut

use base qw/CXGN::ITAG::SeqSourceI Class::Accessor/;

=head1 METHODS

Provides the methods specified in the SeqSourceI interface, using
fake sequences.

=cut

#Class::Accessor provides the new() method

sub all_seq_names {
  my ($self) = @_;

  return ('AC12345','AC12346','AC12347');
}

sub get_seq {
  my ($self, $seqname) = @_;

  my %seqs = ( AC12345 => 'ACTAGCTATATACTGCCTATAACTAGCTAGCATCGATCGATCC',
	       AC12346 => 'ACTGTGCTGCTCGTCGTAGTCGATCAGCTATGACTGATATACGTGCATCGA',
	       AC12347 => 'TTTACGATATAGCGTAGCTGCTGCTAGCATGATGATGCATGACATAATATATGCTAGC',
	     );
  my $newseq =  Bio::PrimarySeq->new(-id => $seqname,
				     -desc       => 'fake sequence',
				     -seq        => $seqs{$seqname},
				    );

  return ( $newseq,
	   [{ ostart => 1,
	      oend => $newseq->length,
	      objname => 'tomato-clone-'.($newseq->display_id),
	      partnum => 1,
	      linenum => 1,
	      type => 'F',
	      typedesc => 'finished',
	      ident => $newseq->display_id,
	      cstart => 1,
	      cend => $newseq->length,
	      orient => '+',
	    }],
	 );

}

=head1 AUTHOR(S)

Robert Buels

=cut

sub DESTROY {
  return parricide(shift, our @ISA);
}

###
1;#do not remove
###
