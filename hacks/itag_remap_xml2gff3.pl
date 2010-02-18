#!/usr/bin/env perl
use strict;
use warnings;
use autodie ':all';

use Carp;
use Getopt::Std;

use 5.10.0;


my ($outfile,@dirs) = @ARGV;
@dirs or die "must provide a list of directories\n";

my $donefile = 'xml2gff.already_done';

my $jobs = parallel_gthxml_to_gff3->new(
    donefile => $donefile,
    outfile  => $outfile,
    dirs     => \@dirs
);

$jobs->max_workers( 12 );
$jobs->run;

exit;


#######  PACKAGES ##########################


BEGIN {
package parallel_gthxml_to_gff3;
use Moose;
use MooseX::Types::Path::Class;
use File::Flock;
use Path::Class;
use autodie ':all';

use POSIX;
use DB_File;

with 'MooseX::Workers';

# file where we store which xml files we have already processed
has 'donefile'   => ( is => 'ro',
                      isa => 'Path::Class::File',
                      required => 1,
                      coerce => 1,
                     );

# append our gff3 output to this file
has 'outfile'    => ( is => 'ro',
                      isa => 'Path::Class::File',
                      required => 1,
                      coerce => 1,
                     );

# keep a hashref of what files were done when we started,
# parsed from our donefile on program start
has 'files_done' => ( is => 'ro',
                      isa => 'HashRef',
                      lazy_build => 1,
                      traits => ['Hash'],
                      handles => {
                          is_done => 'get'
                         },
                     ); sub _build_files_done {
                         my ( $self ) = @_;

                         return {} unless -f $self->donefile;

			 print "indexing donefile ... ";
                         tie my %done, 'DB_File', $self->donefile.'.index', O_CREAT|O_RDWR;
                         my $d = $self->donefile->openr;
                         $done{ $_ } = 1 while <$d>;
			 print "done.\n";
                         return \%done;
                     }

#arrayref of dirs we are searching for files
has 'dirs'  => ( is => 'ro',
                 isa => 'ArrayRef',
                 required => 1,
                );

# run a find process to find the XML files we want
has 'find_handle' => ( is => 'ro',
                       isa => 'FileHandle',
                       lazy_build => 1,
                      ); sub _build_find_handle {
                          open my $f, "ssh eggplant find ".join(' ',map dir($_)->absolute, @{shift->dirs})." -type f -and -name '*.xml' |";
                          return $f;
                      }

# wrap MX::Workers enqueue to just take an xml file name, and generate
# the worker for it
around enqueue => sub {
    my ( $orig, $self, $xml_file ) = @_;

    if( $self->is_done( $xml_file ) ) {
	print "$xml_file skipped\n";
        return;
    }
    print "$xml_file queued ...\n";

    $self->$orig(sub {
        my $gff3_out_file = File::Temp->new;
        $gff3_out_file->close;
        eval {
            GTH_XML_2_GFF3->gthxml_to_gff3( $xml_file, "$gff3_out_file" );
        };
        if( $@ ) {
            warn "$xml_file parse failed:\n$@";
        } else {
            # dump the converted results to our output file
            open my $gff3_fh,'<', "$gff3_out_file";
            { my $l = File::Flock->new( $self->outfile );
              my $out_fh = $self->outfile->open('>>');
	      while( <$gff3_fh> ) {
		  unless( /^##gff-version/ ) {
		      $out_fh->print($_);
		  }
	      }
            }

            # record this file as done
            { my $l = File::Flock->new( $self->donefile );
              $self->donefile->open('>>')->print("$xml_file\n")
	    }
        }
    });
};

sub run {
    my $self = shift;

    # run find first
    my $files = $self->find_handle;

    # print just one gff3 header in the outfile
    $self->outfile->open('>>')->print("##gff-version 3\n");

    # now queue the first set of jobs
    for (1 .. $self->max_workers ) {
        my $xml_file = <$files>;
        chomp $xml_file;
        last unless defined $xml_file;
        $self->enqueue($xml_file);
    }
    POE::Kernel->run;
}

#sub worker_done    { shift; warn join ' ', @_; }
sub worker_done {
    my $self = shift;
    for( 1, 2 ) {
	if( my $xml_file = $self->find_handle->getline ) {
            chomp $xml_file;
	    $self->enqueue( $xml_file );
	}
    }
}

# sub worker_manager_start { warn 'started worker manager' }
# sub worker_manager_stop  { warn 'stopped worker manager' }
# sub max_workers_reached  { warn 'max workers reached' }

sub worker_stdout  { shift; warn join ' ', @_; }
sub worker_stderr  { shift; warn join ' ', @_; }
sub worker_error   { shift; warn join ' ', @_; }
# sub worker_started { shift; warn join ' ', @_; }
# sub sig_child      { shift; warn join ' ', @_; }

no Moose;

package GTH_XML_2_GFF3;
use strict;
use warnings;
use English qw/ -no_match_vars /;

use Bio::SeqIO;
use Bio::FeatureIO;

sub gthxml_to_gff3 {
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
  } else {        	#just name the feature for its source and type
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

} # end of BEGIN block
