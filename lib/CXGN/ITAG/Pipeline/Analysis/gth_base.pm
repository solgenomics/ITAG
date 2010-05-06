package CXGN::ITAG::Pipeline::Analysis::gth_base;
use strict;
use warnings;
use English;
use Carp;

use Fatal qw/ open close chdir /;
use File::Basename qw/basename/;
use File::Path qw/ mkpath rmtree /;

use Bio::FeatureIO;
use Bio::SeqIO;

=head1 NAME

CXGN::ITAG::Pipeline::Analysis::gth_base - pipeline analysis base
class for genomethreader-based analyses

=head1 BASE CLASS(ES)

L<CXGN::ITAG::Pipeline::Analysis>

=cut

use base qw/CXGN::ITAG::Pipeline::Analysis/;

=head1 METHODS

=head2 locally_runnable

returns 1

=cut

#this analysis is locally runnable
sub locally_runnable { 1 }

#and this is the routine to run it
sub run {
    my ($self,$batch) = @_;


    # for each seq file, set up and run an analysis job on the
    # cluster, format its output, and queue it up to be copied into
    # place
    my @jobs;
    my @ops;
    foreach my $seqname ($batch->seqlist) {
        my ($outfile,$gff3_out_file) = map $self->cluster_temp($seqname,$_), 'xml','gff3';

        my $seqfile = $self->_seq_file($batch,$seqname);

        my $gs_est_job = $self->cluster_run_class_method( $batch,
							  $seqname,
							  'run_gth',
                                                          $seqname,
                                                          $seqfile,
							  $outfile,
                                                          $gff3_out_file,
							  #{ err_file => $errfile,
							  #},
                                                         );

        my ($out_dest,$gff_out_dest) = $self->files_for_seq($batch,$seqname);
	push @jobs, $gs_est_job;
	push @ops, [$gff3_out_file => $gff_out_dest ], [ $outfile => $out_dest ];

	# check for any failed jobs, so that we don't have to submit
	# all the jobs before we die.  submitting all the jobs could
	# take a very long time if there are many jobs, because the
	# job submission will block and wait for the torque server to
	# clear a bit if there are too many jobs in its queue.
	$_->alive for @jobs;
    }

    # wait for all the jobs to finish (also running their
    # on_completion hooks)
    sleep 10 while grep $_->alive, @jobs;

    #atomically move the results into position
    $self->atomic_move(@ops);
}

# class method, meant to be run on a cluster node through perl -e
sub run_gth {
    my ( $class, $seqname, $seqfile, $outfile, $gff3_out_file ) = @_;

    my $work_dir = $class->local_temp($seqname);
    mkdir $work_dir;

    local $SIG{TERM} = sub {
	rmtree( $work_dir );
	exit;
    };

    my $un_xed_seqs = $class->local_temp( $seqname, "un_xed.seq" );
    $class->make_un_xed_seqfile( $seqfile => $un_xed_seqs );

    my $cdna_file = Bio::SeqIO->new( -format => 'fasta',
                                     -file   => $class->_cdna_file,
                                   );

    # make xml outfile just a stub
    { open my $outfile_fh, '>', $outfile
          or die "$! writing outfile $outfile";
      $outfile_fh->print( <<EOF );
<?xml version="1.0" encoding="ISO-8859-1"?>
<!-- GenomeThreader XML output no longer provided, only GFF3 -->
EOF
    }

    open my $gff3_out_fh, '>', $gff3_out_file
        or die "$! writing gff3 file $gff3_out_file";
    $gff3_out_fh->print("##gff-version 3\n");


    my $got_seqregion;
    while( my $seq = $cdna_file->next_seq ) {

        my $tempdir = File::Temp->newdir( DIR => $work_dir );

        my $temp_cdna = "$tempdir/cdna";
        Bio::SeqIO->new( -format => 'fasta', -file => ">$temp_cdna" )
                  ->write_seq( $seq );

        my $temp_out  = "$tempdir/xml";

        my @cmd =
            ( 'gth',
              '-xmlout',
              '-force',
              -autointroncutout  => 500,
              '-gcmaxgapwidth'    => 40000, #< max gap in any alignment, this is 2x the largest intron size we have seen so far
              -minalignmentscore => '0.90',
              -mincoverage       => '0.90',
              -seedlength        => 16,
              -o                 => "$temp_out",
              -species           => 'arabidopsis',
              -cdna              => "$temp_cdna",
              -genomic           => $un_xed_seqs,
             );
        ### command: @cmd
        system @cmd
            and die "$! running @cmd";

        my $temp_gff3 = "$tempdir/gff3";

        #now convert the gthxml to gff3
        $class->_gthxml_to_gff3( "$temp_out", $un_xed_seqs, $seqname, "$temp_gff3" );

        #append the gff3 to the overall gff3 file
        open my $g, $temp_gff3 or die 'cannot open temp gff3??';
        while( <$g> ) {
            next if /^##gff-version/ || /^# no results/;
            if( /^##sequence-region/ ) {
                next if $got_seqregion;
                $got_seqregion = 1;
            }
            $gff3_out_fh->print( $_ );
        }
    }

    rmtree( $work_dir );
}


sub _gthxml_to_gff3 {
    my ( $class, $outfile, $un_xed_seqs, $seqname, $gff3_out_file ) = @_;

    my $pm = $class->_parse_mode;
    $pm eq 'alignments'
        and die "_parse_mode() cannot be 'alignments'.  consider using 'alignments_merged'";
    eval {
        my $gth_in = Bio::FeatureIO->new( -format => 'gthxml', -file => $outfile, -mode => $pm );
        my $seqlength = Bio::SeqIO->new( -format => 'fasta', -file => $un_xed_seqs )->next_seq->length;
        my $gff3_out = $class->_open_gff3_out( $seqname, $seqlength, $gff3_out_file );
        while ( my $f = $gth_in->next_feature ) {
	
            # set each feature's source to the name of the gth subclass that's running this
            $class->_recursive_set_source( $f, $class->_source );

            # do additional processing on the feature if necessary
            # (can be implemented in subclasses)
            $class->_process_gff3_feature( $f );

            # make some ID and Parent tags in the subfeatures
            $class->_make_gff3_id_and_parent($f);

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
  my $self = shift;
  die "_source needs to be implemented in ".ref($self);
}

#get the cdna file to run this against
sub _cdna_file {
  my $self = shift;
  die "_cdna_file needs to be implemented in ".ref($self);
}

#get the genomic file to run this against
sub _seq_file {
  my ($self,$batch,$seqname) = @_;
  my ($seq_file) = $batch->pipeline->analysis('seq')->files_for_seq($batch,$seqname);
  $seq_file && -f $seq_file
    or die "expected sequence file '$seq_file' not found";
  return $seq_file;
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
  my ($self,$seqname,$length,$outfile) = @_;

  # handle for out merged output file
  return Bio::FeatureIO->new( -file => ">$outfile",
			      -format => 'gff',
			      -sequence_region => Bio::SeqFeature::Generic->new( -seq_id => $seqname,
										 -start => 1,
										 -end => $length,
									       ),
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

sub make_un_xed_seqfile {
  my ($class,$seqfile,$un_xed_seqs) = @_;

  #don't make the file twice
  unless( -f $un_xed_seqs ) {
    open my $xfile, $seqfile
      or confess "could not open '$seqfile' for reading ($!)";
    open my $nfile, ">$un_xed_seqs"
      or confess "could not open un-xed-seqs file '$un_xed_seqs' for writing";
    while (my $line = <$xfile>) {
      unless($line =~ /^\s*[>#]/) {    #don't munge identifier or comment (comment?) lines
	$line =~ tr/X/N/;
      }
      print $nfile $line;
    }
    close $nfile;
    close $xfile;
  }

  return $un_xed_seqs;
}


=head1 AUTHOR(S)

Robert Buels

=cut

###
1;#do not remove
###
