package CXGN::ITAG::Pipeline::Analysis::seq;
use strict;
use warnings;
use English;
use Carp;

use Bio::SeqIO;

use CXGN::BioTools::AGP qw/agp_write/;

use CXGN::Tools::Class qw/parricide/;
use CXGN::Tools::Identifiers qw/identifier_namespace/;

use CXGN::Tools::List qw/str_in/;

=head1 NAME

CXGN::ITAG::Pipeline::Analysis::seq - pipeline analysis to produce a
fasta file of the sequences to be analyzed

=head1 BASE CLASS(ES)

L<CXGN::ITAG::Pipeline::Analysis>

=cut

use base qw/ CXGN::ITAG::Pipeline::Analysis CXGN::DB::Ima /;

=head1 METHODS

implements all abstract methods defined in L<CXGN::ITAG::Pipeline::Analysis>,
also overrides locally_runnable() and run()'

=cut

# #owner info for the sequence dumping analysis
# sub owner_info {
#   return { user  => 'itagpipeline',
# 	   group => 'users',
# 	   email => 'sgn-bugs@sgn.cornell.edu',
# 	 };
# }

#this analysis is locally runnable
sub locally_runnable { 1 }

# #this analysis does not depend on any others
# sub dependencies { () }

# #this analysis produces one file for each seq
# sub files_for_seq {
#   my ($self,$batch,$seqname) = @_;
#   return $self->assemble_filename($batch,$seqname,'vecscreened','fasta');
# }

=head2 run

  Usage: $an->run($batch)
  Desc : run this analysis, which fetches the sequences named in this
         batch's seqlist.  overrides the run() defined in the base class.
  Args : batch object
  Ret  : nothing meaningful
  Side Effects: writes a fasta file of the seq list to the filesystem

=cut

sub run {
  my ($self,$batch,%args) = @_;
  $batch = $self->_check_batch($batch);

  foreach my $seqname ($batch->seqlist) {
    my ($seqfile,$agpfile) = $self->files_for_seq($batch,$seqname);

    #write the sequence
    my ($seq,$contig) = _get_seq($seqname, $args{seq_source})
      or croak "could not fetch seq named '$seqname'";
    Bio::SeqIO->new(-file => ">$seqfile", -format => 'fasta')
	->write_seq($seq);

    #transform the global coordinates of its AGP lines to be based on
    #the contig instead of the whole chromosome assembly
    #and convert the SOL-style names into genbank accessions
    my $global_offset = $contig->[0]->{ostart} - 1;
    my $part_offset = $contig->[0]->{partnum} - 1;
    foreach my $member (@$contig) {
      $member->{objname} = $seq->id;
      $member->{ostart}  -= $global_offset;
      $member->{oend}    -= $global_offset;
      $member->{partnum} -= $part_offset;
    }

    #now write it into a little agp
    agp_write($contig,$agpfile);
  }
}

sub _get_seq {
  my ($ident,$force_source) = @_;

  return $force_source->get_seq($ident) if $force_source;

  my %sources = ( #tomato_contig => 'CXGN::ITAG::SeqSource::TomatoContigs',
		  genbank_accession => 'CXGN::ITAG::SeqSource::TomatoBACs',
		);
  my $type = identifier_namespace($ident);
  $sources{$type}
    or confess "I don't know how to look up sequences for identifiers of type '$type'";
  unless(ref $sources{$type}) {
    $sources{$type} = $sources{$type}->new;
  }
  return $sources{$type}->get_seq($ident);
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
