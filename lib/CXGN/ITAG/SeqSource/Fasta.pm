package CXGN::ITAG::SeqSource::Fasta;
use strict;
use warnings;
use Carp;

use Path::Class::File ();
use Bio::Index::Fasta;

=head1 NAME

CXGN::ITAG::SeqSource::Fasta - ITAG pipeline sequence source that
provides sequences from an arbitrary fasta-format sequence file

=head1 SYNOPSIS

my $src = CXGN::ITAG::SeqSource::Fasta->new(
              file => '/home/rob/foo.fasta',
          );

=head1 DESCRIPTION

=head1 BASE CLASS(ES)

L<CXGN::ITAG::SeqSourceI>

=cut

use base qw/CXGN::ITAG::SeqSourceI Class::Accessor/;

=head1 METHODS

Provides the methods specified in the SeqSourceI interface, using
tomato BAC sequences.

By default, provides only finished tomato BAC sequences.  However,
if you pass unfinished => 1, will also provide unfinished BAC
sequences.

=cut

#Class::Accessor provides the new() method
__PACKAGE__->mk_accessors(qw/file/);
sub new {
    my $class = shift;
    my $self = $class->SUPER::new( @_ );

    $self->file
        or die "must provide 'file' argument to $class\n";

    # coerce the file arg into a Path::Class::File object
    croak "file '".$self->file."' does not exist"
        unless -f $self->file;
    $self->file( Path::Class::File->new( $self->file ) );

    return $self;
}

sub _index {
    my $self = shift;
    return $self->{index} ||= do {
        # make a tempfile
        my $tempfile = $self->{index_tempfile} = File::Temp->new;
        # write the index to it
        Bio::Index::Fasta->new( -filename => $tempfile->filename,
                                -write_flag => 1,
                              )
                         ->make_index( $self->file );

        # open the index in read mode and return the handle
        Bio::Index::Fasta->new( -filename => $tempfile );
    };
}

sub all_seq_names {
  my ($self) = @_;

  $self->{all_seq_names} ||= do {
      my @names;
      my $fh = $self->file->openr;
      while(<$fh>) {
          next unless /^\s*>\s*(\S+)/;
          push @names, $1;
      }
      \@names
  };

  return @{$self->{all_seq_names}};
}

sub get_seq {
  my ($self, $seqname) = @_;

  my $seq = $self->_index->fetch( $seqname );

  return (

          $seq,

	   [{ ostart => 1,
	      oend => $seq->length,
	      objname => $seq->display_id,
	      partnum => 1,
	      linenum => 1,
	      type => 'F',
	      typedesc => 'finished',
	      ident => $seq->id,
	      cstart => 1,
	      cend => $seq->length,
	      orient => '+',
	    }],
         );
}

=head1 AUTHOR(S)

Robert Buels

=cut

###
1;#do not remove
###
