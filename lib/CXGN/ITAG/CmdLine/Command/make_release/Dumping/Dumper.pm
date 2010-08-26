package CXGN::ITAG::CmdLine::Command::make_release::Dumping::Dumper;
use Moose;
use namespace::autoclean;

use Carp;

use Moose::Util::TypeConstraints;
use MooseX::Types::Moose qw/ ArrayRef Int Str Bool /;
use MooseX::Types::Path::Class qw/ File /;

has 'file' => (
    is       => 'ro',
    required => 1,
    coerce   => 1;
   );

has 'desc' => (
    is => 'ro',
    isa => Str,
);

{ my $st = subtype 'seqtype', as Str,
      where { /^(genomic|protein|cdna|cds)$/ };

  has 'seq_type' => (
      is => 'ro',
      required => 1,
      isa => $st,
     );
}

sub dump {
    my ( $self, $app, $seqlist ) = @_;

    foreach ( @$seqlist ) {
        my ( $seq, $batch ) = @_;
        $self->dump_for_seq( $batch, $seq );
    }
}

sub dump_for_seq { confess 'this method is abstract' }

__PACKAGE__->meta->make_immutable;
1;
