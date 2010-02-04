#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

#use Data::Dumper;

use Bio::SeqIO;
use CXGN::BioTools::AGP qw/agp_write/;

use CXGN::Tools::Identifiers qw/identifier_namespace/;

use CXGN::ITAG::Pipeline;
use CXGN::ITAG::SeqSource::TomatoContigs;
use CXGN::ITAG::SeqSource::TomatoBACs;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script fasta_file ...

  Make a similarly-named spec:agp file for each given fasta file.

  Options:

    none yet

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('',\%opt) or usage();



my $pipe = CXGN::ITAG::Pipeline->open( version => 0 );

foreach my $batch ( map {$pipe->batch($_)} $pipe->list_batches ) {

  foreach my $seqname ($batch->seqlist) {
    my ($seqfile,$agpfile) = CXGN::ITAG::Pipeline::Analysis->open('seq',pipeline => $pipe)->files_for_seq($batch,$seqname);

    #write the sequence
    my ($seq,$contig) = _get_seq($seqname)
      or croak "could not fetch seq named '$seqname'";
#     Bio::SeqIO->new(-file => ">$seqfile", -format => 'fasta')
# 	->write_seq($seq);

    #transform the global coordinates of its AGP lines to be based on
    #the contig instead of the whole chromosome assembly
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
  my ($ident) = @_;
  my %sources = ( tomato_contig => 'CXGN::ITAG::SeqSource::TomatoContigs',
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
