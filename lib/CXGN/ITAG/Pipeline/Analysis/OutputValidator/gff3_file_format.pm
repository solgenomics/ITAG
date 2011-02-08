package CXGN::ITAG::Pipeline::Analysis::OutputValidator::gff3_file_format;
use strict;
use warnings;
use base 'CXGN::ITAG::Pipeline::Analysis::OutputValidator';

use Carp qw/ croak confess /;
use English qw/ -no_match_vars /;

use Data::Dumper;
use File::Basename;

use Storable qw/ nstore retrieve /;

use CXGN::Tools::List qw/ balanced_split /;
use File::Slurp qw/slurp/;
use CXGN::Tools::Run;
use CXGN::Tools::Wget qw/wget_filter/;

use CXGN::ITAG::Pipeline;

sub is_intensive {1}

sub validates_file {
  my ($self,$filename) = @_;
  return $filename =~ /\.gff3$/;
}

sub run_offline {
  my ($self,$analysis,$batch, $options ) = @_;

  warn "running gff3 validate on ".$analysis->tagname."\n";

  #run GFF3 validation
  my @gff3_files = $self->files_to_validate($analysis,$batch)
    or return; #< don't do any of this unless we have some gff3 files


  my @jobs =
      map $self->cluster_run_class_method( $analysis, '_run', $_, { queue => $options->{job_queue} } ),
          $self->_make_jobfiles( $analysis, $batch, \@gff3_files );

  # wait for jobs to finish
  sleep 2 while grep $_->alive, @jobs;

  my $failures = $self->_read_job_failures(\@jobs);

  $self->write_failures( $analysis, $batch, $failures );
}

sub _read_job_failures {
    my ( $self, $jobs ) = @_;

    return { map {
        my $failures = eval $_->out;
        die $@ if $@;
        %$failures;
    } @$jobs }
}

sub _make_jobfiles {
    my ( $self, $analysis, $batch, $files ) = @_;

    # split our files to validate into 100 or fewer batches, which we
    # will allocate to cluster jobs
    my $batches = balanced_split( 100, $files );

    # dispatch the cluster jobs, writing the arguments for each one in a
    # temp file
    my $jobnum;
    return map {
        my $files = $_;

        my $args = { self         => $self,
                     batch        => $batch,
                     analysis_tag => $analysis->tagname,
                     files        => $_,
                    };
        my $arg_file = $analysis->cluster_temp( 'gff3_validate_'.(++$jobnum).'.args' );
        nstore( $args, $arg_file ) or die "$! storing args in $arg_file\n";

        $arg_file
    } @$batches;
}


sub _run {
    my ( $class, $argfile ) = @_;

    my $args = Storable::retrieve( $argfile )
        or die "could not retrieve args from '$argfile'";

    my ( $self, $batch, $analysis_tag, $files ) = @{$args}{qw{ self batch analysis_tag files }};
    my $pipeline = $batch->pipeline;
    my $analysis = $pipeline->analysis( $analysis_tag );

    my $ontology_file = wget_filter( $batch->pipeline->feature_ontology_url
                                        => $analysis->local_temp('ontology_file')
                                       );

    my $val_config_file = $self->validator_config_file;

    my %failures;
    foreach my $file ( @$files ) {

        my $reportfile = $self->report_filename($analysis,$batch,$file);
        my $outbn = $reportfile;
        $outbn =~ s/\.report$// or die "report filename should end in .report ($reportfile)";

        CXGN::Tools::Run->run
          ( 'validate_gff3.pl',
	    -gff3_file => $file,
	    -out       => $outbn,
	    -db_type   => 'sqlite',
	    -config    => "$val_config_file",
	    -ontology_file => $ontology_file,
          );

        unless( slurp($reportfile) =~ /HAS BEEN VALIDATED/ ) {
            #record this failure
            $failures{$file} = time;
        }
    }

    # record the failures
    { local $Data::Dumper::Terse = 1;
      print Dumper( \%failures );
    }

    #try to unlink the ontology tempfile
    unlink $ontology_file;
}

sub cluster_run_class_method {
    my ($class, $analysis, $method, @args) = @_;

    my $run_args = pop @args if ref $args[-1];
    $run_args ||= {};

    $run_args->{temp_base}   ||= $analysis->cluster_temp;
    $run_args->{working_dir} ||= $analysis->cluster_temp;

    ref and die "only bare scalars can be passed as args to cluster_run_class_method"
        for @args;

    my $perl_string = "${class}->${method}(".join( ',', map "'$_'",
                                                   @args
                                                 )
                                        .')';


    return CXGN::Tools::Run->run_cluster
            ( itag_wrapper =>
	      perl =>
              '-M'.$class,
              -e => $perl_string,
              $run_args,
            );

}

{ my $config_file;
  sub validator_config_file {
      return $config_file if $config_file;
      $config_file = File::Temp->new;
      my $user = getpwuid( $UID );
      $config_file->print(<<EOF);  $config_file->close; return $config_file } }
# Configuration file for validate_gff3.pl
# Author: Payan Canaran (canaran\@cshl.edu)
# Copyright 2006 Cold Spring Harbor Laboratory

# DBI database access parameters. If these are filled in, the corresponding
# parameters need not be supplied in the command line.
# *** Parameters supplied in the command line override ones provided here ***

  datasource     
  username       
  password       

# The validation script performs a wide range of checks. Some of these
# checks are fatal and recorded as "errors". However, some of the checks
# can be configured to be recorded as "warnings". The following parameters
# allows for this customization.
#
# In order to customize a check, set it to "W" for warning and to "E" for
# error.

  same_multiple_tags        W

# Data is loaded to the database in chunks rather than individually.
# A buffer is used to store a set number of records before they are
# dumped into the database. This is done to improve loading speed. This
# parameter defines the maximum allowed number of records kept in
# this buffer at any time.

  max_buffer_size           600

# During processing a log file is generated to keep track of progress.
# For some processes, a message is recorded in the log file describing
# how many entries have been processed so far. This is not done for
# every record (as this would generate a very large log file). Instead
# logging is done every "line_number_chunk" records.

  line_number_chunk         10000

# A GFF3 file must have its sources coming from an ontology. The GFF3 specs
# state that, if an ontology (or ontologies) are provided using the
# ##feature-ontology directive, the ontology/ontologies are used. If none is
# provided, the default SOFA ontology is used. The following parameters
# describe the URL for this default ontology.

                            # As of this writing (April 2007) the following SOFA releases are available:
                            # <http://sourceforge.net/project/showfiles.php?group_id=72703&package_id=118027>
                            #
                            # SOFA 2.1 (sofa.obo revision 1.14) - 8/16/2006
                            #  <http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.14>
                            # SOFA 2.0 (sofa.obo revision 1.12) - 5/16/2005
                            #  <http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.12>
                            # SOFA 1.0 (sofa.obo revision 1.6 ) - 5/12/2004
                            #  <http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.6>

  default_sofa_url          http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.14

# This script uses a temporary directory for processing. This parameter,
# describes the location of the temporary directory.

  temp_dir                  /tmp/validate_gff3_$user

# SQLite cache_size pragma

  sqlite_cache_size 512000

# If the GFF3 file contatins an ontology directive that refers to a file that
# needs to be downloaded, the following parameters limit the size of the file 
# that can be downloaded and when the download agent should timeout. 
# Download size is specified in bytes.
# Timeout is specified in seconds.

   max_ontology_download_size 128000
   download_agent_timeout_sec 10
EOF


1;


