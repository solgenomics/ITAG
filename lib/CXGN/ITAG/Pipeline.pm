package CXGN::ITAG::Pipeline;
use strict;
use warnings;
use English;
use Carp;
use Memoize;

use File::Basename;
use File::Spec;
use File::Path;

use CXGN::Tools::List qw/all str_in/;
use CXGN::Tools::Run;

use CXGN::ITAG::Config;

#include this seqsource here since it's the default
use CXGN::ITAG::SeqSource::TomatoBACs;

use CXGN::ITAG::Pipeline::Analysis;
use CXGN::ITAG::Pipeline::Batch;

### INCLUDE EVERYTHING IN THE ANALYSIS DIR
use Module::Find ();
Module::Find::useall('CXGN::ITAG::Pipeline::Analysis');

=head1 NAME

CXGN::ITAG::Pipeline - class representing the ITAG annotation pipeline
at a specific version

=head1 SYNOPSIS

  #open an existing pipeline at version 2.  dies if it doesn't exist
  my $pipe = CXGN::ITAG::Pipeline->open(version => 2);

  #list the analyses that are present in there
  my @as = $pipe->list_analyses;


  ##########

  #make a new pipeline dir structure
  my $newpipe = CXGN::ITAG::Pipeline->create( version => 3,
                                              basedir => '/home/rob/testbase',
                                             );

=head1 DESCRIPTION

coming soon

=head1 BASE CLASS(ES)

L<Class::Accessor>

=cut

use base 'Class::Accessor';

=head1 SUBCLASSES

none yet

=head1 CLASS METHODS

=head2 open

  Usage: my $pipe = CXGN::ITAG::Pipeline->open(version => 1,
                                               basedir => $basedir,
                                              );
  Desc : get an object representing the state of the ITAG pipeline
         at a specific version, using the files in the given base dir
  Args : optional hash-style list of one or more parameters that
         differ from the defaults.  The params and defaults are:
           version  =>  most recent version (highest version number),
           basedir  =>  value of itag_pipeline_base in CXGN configuration
  Ret  : a new object representing the pipeline files at the given
         version and basepath

=cut

sub open {
  my($class, %fields) = @_;

  my $self = $class->_common_new(\%fields);

#  warn "open blessed $self";

  unless(defined $fields{version}) {
    croak "no pipelines found at base dir $fields{basedir}\n";
  }

  -d $self->_pipedir
    or croak "cannot open pipeline structure at version $fields{version}, basedir $fields{basedir}";

  #make sure we have a valid efeature_ontology_url and that it's accessible
  my $sofa_url = $self->feature_ontology_url; #<will die if it's not there

  return $self;
}

memoize('_default_base');
sub _default_base {
  CXGN::ITAG::Config->load->{'itag_pipeline_base'};
}

#given ref to arg hash, set defaults in the arg hash if they're not
#already there
sub _common_new {
  my ($class,$fields) = @_;

  #set defaults
  $fields->{basedir} ||= _default_base();
  $fields->{version} = _most_recent_pipeline_version($fields->{basedir}) unless defined $fields->{version};

  #validate
  croak "invalid version '$fields->{version}'" unless !defined $fields->{version} || $fields->{version} >= 0;
  croak "base dir '$fields->{basedir}' does not exist" unless -d $fields->{basedir};

  return __PACKAGE__->SUPER::new($fields);
}
sub new {
  croak "use open() or create() to make a pipeline object, not new()";
}

=head2 create

  Usage: my $pipe = CXGN::ITAG::Pipeline->create( version => 1);
  Desc : create a new pipeline file structure at the given version
         and basedir
  Args : optional hash-style list, same args as open, except version
         defaults to the next available version number, or 0 if no
         current versions exist
  Ret  : a CXGN::ITAG::Pipeline object for the new pipeline
  Side Effects: creates a new pipeline directory structure,
                copying any existing analysis definitions from the
                previous pipeline version, if the definitions exist
                dies if it already exists

=cut

sub create {
  my ($class,%fields) = @_;

  my $self = $class->_common_new(\%fields);
  $self->version(defined($self->version) ? $self->version + 1 : 0);

  -d $self->_pipedir
    and croak "cannot create, pipeline version $fields{version} at basedir $fields{basedir}, already exists";

  #find the previous version of the pipeline, if any
  my ($most_recent_pipe_version) = sort {$b<=>$a} $class->list_pipelines($fields{basedir});

  #now make our own dirs
  mkdir($self->_pipedir)
    or croak "could not make pipeline main dir ".$self->_pipedir.": $_";

  mkdir($self->_defsdir)
    or croak "could not make pipeline analysis defs dir ".$self->_defsdir.": $_";



  ## copy defs from previous pipeline version
  if(defined $most_recent_pipe_version) {
#    warn "copying over analysis files from $most_recent_pipe_version\n";

    #open the previous pipe. this will die if it fails
    my $oldpipe = $class->open( version => $most_recent_pipe_version,
			        basedir => $fields{basedir}
			      );

    ## copy the global defs file from the previous pipeline version
    CXGN::Tools::Run->run(qw/cp -a/, $oldpipe->_global_defs_filename, $self->_global_defs_filename);

    #for each of the analyses in the previous pipeline, copy over the
    #analysis def file if present
    foreach my $tag ($oldpipe->list_analyses) {
#      warn "copying def file for $tag\n";
      my $old_def_filename = $oldpipe->_analysis_def_filename($tag);
      if(-f $old_def_filename) {
	CXGN::Tools::Run->run(qw/cp -a/, $old_def_filename, $self->_analysis_def_filename($tag));
      }
    }
  }

  return $self;
}

#return latest pipeline version number, or undef if none found
sub _most_recent_pipeline_version {
  my ($basedir) = @_;

  my @pipedirs = __PACKAGE__->list_pipelines($basedir);

  return $pipedirs[-1] if @pipedirs;
  return;
}

sub _largest_existing_batchnum {
  my ($self) = @_;
  my @batches = $self->list_batches(include_deleted => 1);
  return @batches ? $batches[-1] : undef;
}

=head2 username

  Usage: CXGN::ITAG::Pipeline->username
  Desc : return the name of the pipeline system user.
         usable either as a class or object method
  Ret  : currently, 'itagpipeline'
  Args : none

=cut

sub username { 'itagpipeline' }

=head1 OBJECT METHODS

=head2 basedir

  Usage: my $dir = $pipe->basedir
  Desc : get the base directory for this pipeline object
  Args : none
  Ret  : string directory name

=head2 version

  Usage: my $ver = $pipe->version
  Desc : get the pipeline version for this pipeline
  Args : none
  Ret  : numerical pipeline version

=cut

__PACKAGE__->mk_accessors qw( basedir version );



#get the dir for our version of the pipeline
memoize '_pipedir';
sub _pipedir {
  my ($self) = @_;
  return File::Spec->catdir($self->basedir,sprintf('pipe%03d',$self->version));
}
memoize '_defsdir';
#get the dir for analysis defs in this pipeline
sub _defsdir {
  my ($self) = @_;
  return File::Spec->catdir($self->_pipedir,'analysis_defs');
}

memoize '_analysis_def_filename';
#get the def filename for a specific analysis tag
sub _analysis_def_filename {
  my ($self,$tag) = @_;
  return File::Spec->catfile($self->_defsdir,"$tag.def.txt");
}
memoize '_analysis_def_files';
#get a list of basenames of all analysis def files in this pipeline
sub _analysis_def_files {
  my ($self) = @_;
  return map {basename($_)} glob $self->_defsdir.'/*.def.txt';
}
#get the dir in which to put things that have been deleted from this
#pipeline
memoize '_deleted_subdir';
sub _deleted_subdir {
  my ($self) = @_;
  my $deldir = File::Spec->catdir($self->_pipedir,'.deleted');
  unless(-d $deldir) {
    mkdir $deldir
      or confess "could not create $deldir: $!";
  }
  return $deldir;
}
#get the full path to the global defs file
sub _global_defs_filename {
  my ($self) = @_;
  return File::Spec->catfile($self->_pipedir,'global.def.txt');
}
# returns a LIST Of the values associated with the given key
sub _global_defs_values {
  my ($self,$key) = @_;
  my $defs_file = $self->_global_defs_filename;
  my $p = eval{CXGN::ITAG::Tools::parse_kv_file($defs_file)}; #<might die
  $p ||= {};
  return unless $p->{$key};
  return @{$p->{$key}};
}

#return the filename of this pipeline's email notification log
sub email_log_filename {
  my ($self) = @_;
  return File::Spec->catfile($self->_pipedir,'.email_log');
}



=head2 list_analyses

  Usage: my @a = $pipe->list_analyses
  Desc : get the list of analyses (tags) that seem to be defined in this
         pipeline, based on the contents of the pipeline dir

         NOTE:

               this is not the same as list_defined_analyses() in the
               analysis base class.  that one lists all analyses that
               are currently defined in perl code or def files, while
               this method lists analyses that seem to be present in
               this specific (version of the) pipeline, based on the
               analysis directories and def files present on the
               filesystem

  Args : none
  Ret  : list of analysis tag names, if any

=cut

sub list_analyses {
  my ($self) = @_;
  return map {s/\..+$//; $_} $self->_analysis_def_files;
}

=head2 list_batches

  Usage: my @batchnums = $pipe->list_batches
  Desc : list the batch numbers that currently exist in this pipeline
  Args : none
  Ret  : possibly empty list of (zero-padded) batch numbers, sorted in
         ascending numerical order
  Side Effects: lists files in the filesystem

=cut

sub list_batches {
  CXGN::ITAG::Pipeline::Batch->list_batches(@_);
}

=head2 list_pipelines

  Usage: my @pvers = CXGN::ITAG::Pipeline->list_pipelines($basedir)
  Desc : list the versions of all the pipelines present in the given directory
  Args : (optional) base path to look in, defaults
  Ret  : list of pipeline version numbers in ascending order (e.g. (0,1,2))
  Side Effects: none

=cut

sub list_pipelines {
  my ($self,$basedir) = @_;
  $basedir ||= _default_base();

  return sort {$a <=> $b}
    map {
      m!(\d+)/$!;
      $1 + 0 #< get rid of leading zeros
    } glob "$basedir/pipe[0-9][0-9][0-9]/";
}


=head2 analysis

  Usage: my $eug = $pipe->analysis('eugene');
  Desc : get the CXGN::ITAG::Pipeline::Analysis object with the given tag
  Args : analysis tag name, batch number
  Ret  : analysis object, or undef if there is no such analysis

=cut

memoize('analysis');
sub analysis {
  my ($self,$tag,%args) = @_;
#  use Data::Dumper; print Dumper(\@_);
  $tag or croak 'must give analysis tag to analysis() method';
  my $a = eval { CXGN::ITAG::Pipeline::Analysis->open($tag, pipeline => $self) };
  if($EVAL_ERROR) {
    warn $EVAL_ERROR;
    return;
  }
  return $a;
}

=head2 analysis_def_file

  Usage: my $filename = $pipe->analysis_def_file('repeats');
  Desc : get the filename of the analysis's definition file
  Args : an analysis tag name
  Ret  : a string filename, which may or may not exist

=cut

sub analysis_def_file {
  my ($self,$tag) = @_;

  return $self->_analysis_def_filename($tag);
}


=head2 batch

  Usage: my $batch = $pipe->batch(1);
  Desc : get the Batch object for the given batch number
  Args : batch number
  Ret  : a L<CXGN::ITAG::Pipeline::Batch> object for the batch,
         or undef if the batch does not exist
  Side Effects: none

=cut

memoize('batch');
sub batch {
  my ($self,$batchnum) = @_;
  defined $batchnum or confess 'must give a batch number as an argument';
  $batchnum =~ /\D/ and confess "batch number must be numeric";
  $batchnum >= 1 or confess "batch number must be >= 1";

  return $self->{batches}[$batchnum] ||= CXGN::ITAG::Pipeline::Batch->open(pipeline => $self, batchnum => $batchnum);
}

=head2 create_batch

  Usage: $pipe->create_batch
  Desc : create a new pipeline batch
  Args : optional hash-style list as:
           batchnum => batch number,
                       defaults to next available batch number,
           seqlist  => [arrayref of seqs to include in the batch],
                       defaults to a list of sequences that have not been
                       batched yet,
           seq_source => L<CXGN::ITAG::SeqSourceI>-implementing object to use for
                        fetching sequences.  Defaults to
                        L<CXGN::ITAG::SeqSource::TomatoBACs>,
                        instantiated with default options,
           size     => number of seqs to include in the batch,
                       default 100,
           no_bacs_in_contigs => if set, do not consider clone seqs eligible for
                                 membership in the new batch if they have been
                                 previously annotated as part of a contig.
                                 Default false.
           include_all => do not exclude sequences that are already part
                          of another batch.  Default 0
  Ret  : a L<CXGN::ITAG::Pipeline::Batch> object for the new pipeline batch,
         or undef if no batch could be created
  Side Effects: creates a new batch directory structure in the file system,
                puts sequence list in place to start the batch

=cut

sub create_batch {
  my ($self,%args) = @_;

  $args{size} ||= 100;

  return CXGN::ITAG::Pipeline::Batch->create(pipeline => $self, %args);
}

=head2 delete_batch

  Usage: $pipe->delete_batch(1);
  Desc : deletes a pipeline batch
  Args : optional batch number to delete, defaults to latest available batch
  Ret  : nothing meaningful
  Side Effects: dies on error

=cut

sub delete_batch {
  my ($self,$batchnum) = @_;
  $batchnum+0 eq $batchnum
    or croak "batch number must be numeric\n";

  my @batches = $self->list_batches;
  my $batch = do {
    if($batchnum) {
      str_in($batchnum+0,@batches)
	or croak "batch number $batchnum does not exist";
      $self->batch($batchnum);
    } else {
      $self->batch($batches[-1]);
    }
  };

  $batchnum = $batch->batchnum;
  $batch->delete;
  delete $self->{batches}[$batchnum]; #delete from cache too
}


=head2 valid_file_exts

  Usage: my @exts = $pipe->valid_file_exts
  Desc : get a list of valid file extensions for files
         in this pipeline
  Args : none
  Ret  : list of strings, possibly empty

=cut

sub valid_file_exts {
  my ($self) = @_;
  return $self->_global_defs_values('file_exts');
}

=head2 feature_ontology_url

  Usage: my $url = $pipe->feature_ontology_url
  Desc : get a URL for the OBO-format ontology file that
         should be used to validate the GFF3 files in this
         version of the pipeline
  Args : none
  Ret  : a URL string
  Side Effects: dies on error

=cut

sub feature_ontology_url {
  my ($self) = @_;
  my ($url) = $self->_global_defs_values('feature_ontology_url');
  #warn "got url '$url'";
  $url or die "no feature_ontology_url defined in global defs file ".$self->_global_defs_filename."\n";
  return $url;
}

########### UTILITY SUBCLASSES

=head1 MAINTAINER

Robert Buels

=head1 AUTHOR

Robert Buels, E<lt>rmb32@cornell.eduE<gt>

=head1 COPYRIGHT & LICENSE

Copyright 2009 Boyce Thompson Institute for Plant Research

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

###
1;#do not remove
###

