#!/usr/bin/env perl
use strict;
use warnings;
use autodie ':all';
use Path::Class;

use Carp;

use Getopt::Std;

use 5.10.0;

my @dirs = @ARGV
    or die "must provide a list of directories\n";

my $donefile = 'already_done';

my %already_done = do {
    if( -f $donefile ) {
	open my $d, $donefile;
	map { chomp; $_ => 1 } <$already_done>
    }
};

open my $files, "find @dirs -type f -and -name '*.xml' |";
open my $done_fh, ">>", $donefile;
while( my $xml_file = <$files> ) {
    chomp $xml_file;
    if( $already_done{$xml_file} ) {
	warn "file $xml_file already done, skipping.\n";
	next;
    }

    my $gff3_out_file = File::Temp->new;
    $gff3_out_file->close;
    GTH_XML_2_GFF3->_gthxml_to_gff3( $xml_file, "$gff3_out_file" );
    my $fh = file( "$gff3_out_file" )->openr;
    print while <$fh>;
    print $done_fh "$xml_file\n";
}

#use Data::Dumper;

# my %opt;
# getopts('',\%opt) or pod2usage(1);

package GTH_XML_2_GFF3;
use strict;
use warnings;
use English qw/ -no_match_vars /;

use Bio::SeqIO;
use Bio::FeatureIO;

sub _gthxml_to_gff3 {
    my ( $class, $outfile, $gff3_out_file ) = @_;

    my $pm = $class->_parse_mode;
    $pm eq 'alignments'
        and die "_parse_mode() cannot be 'alignments'.  consider using 'alignments_merged'";
    eval {
        my $gth_in = Bio::FeatureIO->new( -format => 'gthxml', -file => $outfile,
					  -mode => $pm,
					  -attach_alignments => 1,
					 );
        my $gff3_out = $class->_open_gff3_out( $gff3_out_file );
        while ( my $f = $gth_in->next_feature ) {
	
            # set each feature's source to the name of the gth subclass that's running this
            $class->_recursive_set_source( $f, $class->_source );

            # do additional processing on the feature if necessary
            # (can be implemented in subclasses)
            $class->_process_gff3_feature( $f );

            # make some ID and Parent tags in the subfeatures
            $class->_make_gff3_id_and_parent($f);

	    # stringify the supporting_alignment's to their IDs
	    $f->add_Annotation('supporting_alignment',$_)
		for map { $_->value->get_Annotations('ID') }
		    $f->remove_Annotations('supporting_alignment');

            # write the feature to the gff3 file
            $gff3_out->write_feature($f);
        }
    }; if( $EVAL_ERROR ) {
        #workaround for a gth bug.  will probably be fixed when we upgrade gth
        die $EVAL_ERROR unless $EVAL_ERROR =~ /not well-formed \(invalid token\)/;
        open my $gff3, '>', $gff3_out_file;
        print $gff3 <<EOF;
##gff-version 3
# no results.  genomethreader produced invalid output XML.
EOF
        open my $out, '>', $outfile;
    }
}


# does nothing here, but may be implemented in subclasses
sub _process_gff3_feature {
}

#get the source name to use in our GFF3
sub _source {
    'GTH'
}

#recursively set the source on a feature and its subfeatures
sub _recursive_set_source {
  my ($self,$feature,$newsource) = @_;
  $feature->source($newsource);
  $self->_recursive_set_source($_,$newsource) for $feature->get_SeqFeatures;
}

# object OR class method to
# open a gff3 outfile with the right version and sequence-region
sub _open_gff3_out {
  my ($self,$outfile) = @_;

  # handle for out merged output file
  return Bio::FeatureIO->new( -file => ">$outfile",
			      -format => 'gff',
			      -version => 3,
                            );
}


#take a feature hierarchy, manufacture ID and Parent tags to encode
#the hierarchical relationships, adding them to the features
sub _make_gff3_id_and_parent {
  my ($class,$feat,$parent_ID) = @_;

  $feat->add_Annotation('Parent',Bio::Annotation::SimpleValue->new(-value => $parent_ID))
    if defined $parent_ID;

  #make a unique id for this thing, keeping our id counters on a
  #per-analysis level
  if(my $idstem = $class->_feature_id($feat,$parent_ID)) {
    my $uniqid = $class->_unique_bio_annotation_id($idstem);
    $feat->add_Annotation('ID',Bio::Annotation::SimpleValue->new(-value => $uniqid));
    #recursively ID and Parent all the subfeatures, if any
    $class->_make_gff3_id_and_parent($_,$uniqid) for $feat->get_SeqFeatures;
  }
}

#given a stem, make a ID that's unique to this analysis
#by appending a number to the stem
my %uniq_id_ctrs;
sub _unique_bio_annotation_id {
  my ($class,$idstem)  = @_;
  $uniq_id_ctrs{$class} ||= {};
  return Bio::Annotation::SimpleValue->new(-value => $idstem.'_'.++$uniq_id_ctrs{$class}{$idstem});
}

sub _feature_id {
  my ($class,$feat,$parent_ID)  = @_;
  if($feat->type->name eq 'mRNA') {
    "${parent_ID}_AGS"
  } elsif ( $feat->type->name eq 'match') {
    #get the target name of the first subfeature's target
    my ($target_id) = (($feat->get_SeqFeatures)[0]->get_Annotations('Target'))[0]->target_id;
    $target_id.'_alignment'
  } elsif ( $feat->type->name eq 'region') {
    'PGL'
  } else {			#just name the feature for its source and type
    my $src = $feat->source;
    $src =~ s/GenomeThreader/GTH/; #< shorten the sources a little
    $src =~ s/(tom|pot)ato/$1/; #< shorten the sources a little
    $src.'_'.$feat->type->name;
  }
}

#can be overridden
sub _parse_mode {
  'both_merged';
}
