package CXGN::ITAG::Pipeline::Analysis;
use strict;
use warnings;
use English qw( -no_match_vars );

use Carp;

use Digest::MD5;

use File::Basename;
use File::Copy ();
use File::NFSLock ();
use File::Path;
use File::Spec;
use File::Temp qw/tempdir/;

use Memoize qw/memoize/;
use Module::Find ();

use List::MoreUtils qw/ any all /;
use CXGN::Tools::List qw/ str_in distinct/;
use CXGN::ITAG::Tools;

use CXGN::ITAG::Config;

=head1 NAME

CXGN::ITAG::Pipeline::Analysis - base class for ITAG pipeline
analyses.  These classes are only used by the L<CXGN::ITAG::Pipeline>
object

=head1 SYNOPSIS

coming soon

=head1 DESCRIPTION

coming soon

=cut

use base qw/Class::Accessor::Fast/;

=head1 SUBCLASSES

all the ITAG pipeline analyses

=head1 CLASS METHODS

=head2 open, aopen

  Usage: my $an = CXGN::ITAG::Analysis->open('my_analysis',
                                             pipeline => $pipe);
  Desc : open an existing analysis in a pipeline
  Args : tag name of analysis to open,
         hash-style list as:
           pipeline => L<CXGN::ITAG::Pipeline> object this analysis belongs to
  Ret  : a new analysis object
  Side Effects: dies if there was an error opening the analysis

=cut

memoize('aopen');
sub open { shift->aopen(@_) }
sub aopen {
  my ($class,$tag,%args)  = @_;
  my $apackage = _tag2pkg($args{pipeline},$tag);

  my $self = bless \%args,$apackage;

  my $def_file = $self->_filename('def');
  $def_file && -f $def_file
    or die "no def file '$def_file' found for analysis $tag\n";

  my $p = eval {$self->_parse_kv_file('def')};

  return $self;
}

#check for the existence and proper permissions of our analysis dir.
#return a list of errors if anything's wrong
sub _check_dir {
  my ($self,$batch) = @_;
  $batch = $self->_check_batch($batch);

  my @errors;

  #check that our dir exists
  my $dir_exists = -d $self->work_dir($batch)
    or push @errors,'analysis directory does not exist for analysis '.$self->tagname.', batch '.$batch->batchnum.', pipeline '.$self->pipeline->version.'.  Do you need to recreate the analysis dirs for this batch?';

  unless( !$dir_exists || $ENV{CXGNITAGPIPELINEANALYSISTESTING} ) {
    #check the permissions on our dir, die if they're wrong
    my ($dir_uid,$dir_gid) = (stat($self->work_dir($batch)))[4,5];
    my $exp_string =  "expecting grp ".$self->owner_info->{group}.", user ".$self->owner_info->{user}."\n";

    my $pipeuser_name = $self->pipeline->username;
    my $pipeuser_uid = getpwnam($pipeuser_name)
      or push @errors,("pipeline user $pipeuser_name does not exist on this system.");

    my $groupname = $self->owner_info->{group}
      or push @errors,("no owner group found");
    my $gid = getgrnam($groupname)
      or push @errors,("group '$groupname' for analysis ".$self->tagname." doesn't exist on this system.  aborting");
    $dir_uid == $pipeuser_uid && $dir_gid == $gid
      or push @errors,("incorrect permissions on ".$self->tagname." analysis dir ".$self->work_dir($batch).", $exp_string");

  }

  return @errors;
}

=head2 create_dirs

  Usage: my $an = CXGN::ITAG::Pipeline::Analysis->create('foo',batch => $mybatch);
  Desc : create or recreate the directory structure for an analysis,
         also setting the correct permissions
  Args : hash-style list as:
           batch => batch object this analysis belongs to
  Ret  : a new analysis object
  Side Effects: dies if the dir structure already exists

=cut

sub create_dirs {
  my ($self,$batch) = @_;
  $batch = $self->_check_batch($batch);

  my $dir = $self->work_dir($batch);

  -d $dir || mkdir $dir
    or die "could not make analysis dir $dir: $!\n";

  # skip dir creating if running under unit tests
  return $self if $ENV{CXGNITAGPIPELINEANALYSISTESTING};

  # look up the user and group names and ids we are supposed to have
  my $group_name = $self->owner_info->{group} or die 'no group name defined?';
  my (undef,undef,$group_id) = getgrnam( $group_name )
    or croak "could not find GID for group ".$self->owner_info->{group};
  my $user_name = $self->pipeline->username or die 'no pipeline username defined?';
  my (undef,undef,$user_id) = getpwnam($self->pipeline->username)
    or croak "could not find UID for user ".$self->pipeline->username;

  return if _dir_perms_ok( $user_id, $group_id, $dir );

  chown $user_id, $group_id, $dir;
  chmod 0775, $dir;

  _dir_perms_ok( $user_id, $group_id, $dir )
    or croak  "could not chown $dir to $user_name.$group_name";

  return $self;
}
sub _dir_perms_ok { #little helper function that checks user and group perms
  my ($uid,$gid,$dir) = @_;
  my @s = stat $dir;
  return $s[4] == $uid && $s[5] == $gid  &&  ($s[2] & 07777) == 0775;
}


=head1 OBJECT METHODS

=head2 owner_info

  Usage: my $owning_user = $analysis->owner();
  Desc : get info about the owner of this analysis.  unless overridden,
         reads this info from the analysis definition file.
  Args : none
  Ret  : hashref as:
         {   user  => 'itagimperial',
             group => 'itagimperial',
             email => 'dbuchan@imperial.ac.uk',
         }
  Side Effects: dies if the analysis definition file cannot be read
                or is not properly formatted

=cut

sub owner_info {
  my ($self) = @_;
  my $h =  {map {my ($v) = $self->_kv_file_value(def => "owner_$_"); $_ => $v} qw/user group email/};
  $h->{group} ||= $h->{user};

  return $h;
}

=head2 dependencies

  Usage: my @deps = $analysis->dependencies;
  Desc : get the list of analyses this one depends on for input
  Args : none
  Ret  : list of tagnames of analyses this one depends on for
         its input
  Side Effects: none

=cut

sub dependencies {
  my ($self) = @_;
  my @deps = $self->_kv_file_value(def => 'depends_on');
  return grep {$_ ne $self->tagname} @deps;
}

=head2 files_for_seq

  Usage: my @files = $analysis->files_for_seq($batch,$seqname);
  Desc : get the list of expected filenames this analysis should
         produce for a sequence with the given name
  Args : batch object, a single sequence name
  Ret  : a list of files that this analysis should produce for
         the given sequence in the given batch
  Side Effects: none

=cut

memoize 'files_for_seq';
sub files_for_seq {
  my ($self,$batch,$seqname) = @_;
  UNIVERSAL::isa($batch,'CXGN::ITAG::Pipeline::Batch')
      or confess "must give batch object to files_for_seq";
  my @filedefs = map {
    my ($desc,$ext) = split /:/,$_;
    $desc && $ext
      or die "error parsing produces_files line in def file for ".$self->tagname.", disabling analysis\n";
    $desc =~ /^[a-z0-9_]+$/
      or die "'$desc' is not a valid file description.  Must contain only lower-case letters, numbers, and underscores\n";
    $self->assemble_filename($batch,$seqname,$desc,$ext)
  } $self->output_files_specs;

  return @filedefs;
}

=head2 output_files_specs

  Usage: my @specs = $an->output_files_specs;
  Desc : get the output specification strings from this analysis's definition file
  Args : none
  Ret  : list of strings like 'foo:txt','bar:fasta', etc
  Side Effects: none

=cut

memoize 'output_files_specs';
sub output_files_specs {
  my ($self) = @_;
  return $self->_kv_file_value(def => 'produces_files');
}


=head2 expected_output_files

  Usage: my @files = $an->expected_output_files($batch);
  Desc : get a list of filenames expected to be in the output
         dir of this analysis when it's done
  Args : L<CXGN::ITAG::Pipeline::Batch> object
  Ret  : list of string filenames

=cut

memoize 'expected_output_files';
sub expected_output_files {
  my ($self,$batch) = @_;
  $batch = $self->_check_batch($batch);
  my @seqlist = $batch->seqlist;
#   use Data::Dumper;
#   warn 'got seqlist '.Dumper(\@seqlist);
  return map {
#    warn "getting files for '$_'\n";
    $self->files_for_seq($batch,$_)
  } @seqlist;
}

=head2 pipeline

  Usage: my $pipe = $an->pipeline;
  Desc : get/set the pipeline object this analysis is
         associated with
  Args : (optional) new pipeline to set
  Ret  : the current L<CXGN::ITAG::Pipeline> object
  Side Effects: none

=cut

__PACKAGE__->mk_accessors('pipeline');


=head2 status

  Usage: my $stat = $an->status($batch);
  Desc : get the status of this analysis, based on the files
         in the filesystem
  Args : L<CXGN::ITAG::Pipeline::Batch> object, or batch number
  Ret  : one of ('not_ready','ready','running','done','error','disabled')
  Side Effects: looks in the filesystem
  Example:

     print join ' ',$an->tagname,'is in',$an->status,"status\n";

=cut

# we cache the result of the status method for a couple of seconds
memoize 'status',
  NORMALIZER =>  sub {
    my ($self,$batch) = @_;
    # status depends on:
    join ' ',( $self,   #< the analysis tagname
	       ref $batch ? $batch->batchnum : $batch, #< the batch number
	       (stat $self->_filename('control',$batch))[9] || '', #< the last update of its control file
	     )
  };

sub _list_work_dir {
  my ($self,$batch) = @_;
  my $d = $self->work_dir($batch);
  return unless -d $d;
  opendir my $dh, $d or confess "$! reading dir $d";
  return map File::Spec->catfile($d,$_), grep !/^\./, readdir $dh;
}

sub is_disabled {
    shift->_kv_file_value( def => 'disabled');
}

sub status {
  my ($self,$batch) = @_;
  $batch = $self->_check_batch($batch);
  #warn('status called for '.$self->tagname."\n");

  return 'disabled' if $self->is_disabled;

  return 'error' if $self->errors($batch);

  my ($running) = do {
    if(-f $self->_filename('control',$batch) ) {
      $self->_kv_file_value( control => 'running', $batch );
    }
  };
  return 'running' if $running;

  #are there any files in there?
  my %special_files = map {$_ => 1} $self->_adir_special_files($batch);
  if(my @files_in_dir = grep !$special_files{$_},
                        $self->_list_work_dir($batch)
    ) {

    #are all the required files there?
    my @exp_outfiles = $self->expected_output_files($batch);
    my %files_in_dir = map {$_ => 1} @files_in_dir;
    if( @exp_outfiles && all {$files_in_dir{$_}} @exp_outfiles ) {
      #do they need validation checks to be run on them?
      return 'validating' if any  {
	$_->needs_update($self,$batch)
      } $self->_list_validator_packages;

      return 'done';
    } else { #don't have all the required files
      $self->_store_error("analysis directory is not empty, but not all specified output files are present");
      return 'error';
    }
  } else {
    #OK, no files
    #is it ready to run?
    my @deps = $self->dependencies;
    my $class = ref $self or die "sanity check failed";
    return 'ready'
      if !@deps || all {
	  my $atag = $_;
	  #if the dependency is an analysis that exists, check whether that analysis is done
	  my $dep_a = eval { $self->aopen( $atag, pipeline => $self->pipeline) };
	  $dep_a && $dep_a->status($batch) eq 'done'
      } @deps;

    return 'not_ready';
  }

  die 'sanity check failed, this point should not be reached';
}


=head2 force_validate

  Usage: $an->force_validate( $batch )
  Desc : force an analysis in a certain batch to be considered validated.
  Args : batch object
  Ret  : nothing meaningful
  Side Effects: dies on error

=cut

sub force_validate {
    my ( $self, $batch ) = @_;


    for my $validator ( $self->_list_validator_packages ) {
	$validator->force_validate( $self, $batch );
    }
}

=head2 uncached_status

  Usage:  same as status()
  Desc :  sames as status() above, but does not cache any of the
          status information.  You'll want to use this instead of
          status() if you are writing something that calls status()
          multiple times on an analysis, for example, if you're
          polling for completion.
  Args :  same as status()
  Ret  :  same as status()

=cut

sub uncached_status {
  my $self = shift;
  Memoize::flush_cache('status');
  return $self->status(@_);
}


=head2 tagname

  Usage: my $tagname = $a->tagname;
  Desc : get the tagname for an analysis object, or analysis class.
         This method can be used as either a class method OR
         an object method
  Args : none
  Ret  : this analysis's tag name, which is the last
         part of its package name.
  Side Effects: none
  Example:

    my $a = CXGN::ITAG::Pipeline::Analysis::foo->open(batch => $mybatch);
    print $a->tagname,"\n"; #prints 'foo'

=cut

memoize 'tagname';
sub tagname {
  my ($self) = @_;
  my $packname = ref($self) || $self;
  my @s = split /::/,$packname;
  $s[-1] =~ /^[a-z]\w+$/
    or confess "$s[-1] is not a valid tagname for an analysis.  tagnames must be all lower-case letters, numbers, and underscores\n";
  return $s[-1];
}

=head2 locally_runnable

  Usage: $an->run if $an->locally_runnable;
  Desc : get whether this analysis is runnable locally
         on CXGN systems.  usable as either a class or object method.
  Args : none
  Ret  : false.  override if you want yours to return true.

=cut

sub locally_runnable {
  return 0;
}


=head2 run

  Usage: my $result = $an->run( $batch )
  Desc : if this analysis is locally runnable, execute the analysis.
         this should produce the correct output files in the correct places.
         override this if your analysis is locally runnable.
  Args : L<CXGN::ITAG::Pipeline::Batch> object, hashref of options as:
         {  job_queue => `qsub` queue for submitting analysis jobs, uses
                         the system's default if unset
         }
  Ret  : nothing meaningful
  Side Effects: varies

=cut

sub run {
  confess 'run() is abstract, it must be implemented in your subclass';
}

=head2 resource_estimate

  Usage: my $res_hashref = $an->resource_estimate( $batch, $seqname );
  Desc : get an estimate of resources required to to run the analysis
         on the given sequence in the given batch
  Args : batch object or number, sequence name
  Ret  : possibly empty hashref of the form:
          { vmem => total memory in megabytes,
          }
  Side Effects: may download dataset files or create temp dirs

  The base implementation simply returns no information.  Concrete
  analysis classes should override this.

=cut

sub resource_estimate {
    my ( $self, $batch, $seqname ) = @_;
    return {}
}

=head2 errors

  Usage: my @errors = $an->errors
  Desc : get descriptive strings for errors present in an analysis
  Args : optional current batch object.  if not provided,
         does not check anything related to any specific batch.
  Ret  : (possibly empty) list of strings
  Side Effects: none

=cut

sub errors {
  my ($self,$batch) = @_;
  eval {$self->_parse_kv_file('def')};

  $self->{stored_errors} ||= [];
  my @errors = distinct @{$self->{stored_errors}};

  if($batch) {
    $batch = $self->_check_batch($batch);
    push @errors, $self->_check_dir($batch);
    push @errors, $self->output_validation_errors($batch);
  }

  return @errors;
}

sub _store_error {
  my ($self,$e) = @_;
  $self->{stored_errors} ||= [];
  push @{$self->{stored_errors}},$e;
}


=head2 output_validation_errors

  Usage: my @errors = $an->output_validation_errors;
  Desc : object method to validate the output files of an analysis
  Args : none
  Ret  : list of strings, each of which describes an individual
         error in this analysis's output.  method provided
         here just checks the manifest against the files in the dir

=cut

sub output_validation_errors {
  my ($self,$batch) = @_;

  my @errors;

  #check for any unexpected files in the dir
  my %expected_files = map {$_ => 1} $self->expected_output_files($batch), $self->_adir_special_files($batch);
  foreach my $file ( $self->_list_work_dir($batch) ) {
    $expected_files{$file}
      or do{ my $bn = basename($file); push @errors, "unexpected file $bn, not specified in ".$self->tagname." def file"; next};
  }

  #check the results of all the output validators
#    warn "checking $outfile...";
  foreach my $validator_package ($self->_list_validator_packages) {
    #      warn "with $validator_package\n";
    push @errors, $validator_package->errors($self,$batch);
  }

  return @errors;
}

sub _list_validator_names {
  my ($self) = @_;
  return map _valname($_), $self->_list_validator_packages;
}
sub _valname {
    return $1 if $_[0] =~ /::([^:]+)$/;
    return;
}
sub _list_validator_packages {

    # list validation modules from more to less specific
    my @modules = ( Module::Find::useall( ref(shift).'::OutputValidator' ),
                    Module::Find::useall( __PACKAGE__.'::OutputValidator' ),
                   );
    # remove duplicates (can't have 2 md5sum.pm for example)
    my %seen;
    return grep !$seen{_valname($_)}++, @modules;

}

=head2 run_intense_validation

  Usage: $an->run_intense_validation($batch);
  Desc : run computationally-intensive validation routines
         on this analysis.  unless the force flag is true,
         skips analyses that report they do not need to be
         run again yet
  Args : batch object,
        (optional) hashref of options:
        { force     => boolean, default false,
          job_queue => queue for submitting Torque/SGE jobs,
                       uses system default if not set,
        }
  Ret  : nothing
  Side Effects: dies on failure

=cut

sub run_intense_validation {
  my ($self,$batch,$options) = @_;
  #warn "intense_validate on ".$self->tagname." in batch ".$batch->batchnum.", with force '$force'\n";
  my %options = %{$options ||= {}};

  my $force = delete $options{force};

  #run all the validators that are intensive
  foreach my $valpackage (grep {$_->is_intensive} $self->_list_validator_packages) {
    #print $self->tagname." need to run ".$valpackage->name."?\n";
    next unless $force ||  $valpackage->needs_update($self,$batch);
    #print "yes.\n";
    $valpackage->run_offline($self,$batch,\%options);
  }
}

#each validation has a name, which is a dir inside the cache dir
#in there, you can have val_failed files


=head2 check_manifest

  Usage: print 'yep' if $an->check_manifest;
  Desc : check the upload manifest in the dir,
         if included
  Args : none
  Ret  : ('none') if no manifest,
         list of error strings if mismatches found,
         empty list if manifest checks out OK
  Side Effects: looks in the filesystem

=cut

sub check_manifest {
  my ($self,$batch) = @_;
}

#convert batch numbers to batch objects if necessary,
#check that this batch object is from the same pipeline
memoize '_check_batch';
sub _check_batch {
  my ($self,$batch) = @_;

  $batch or confess "must give a batch number or object to _check_batch";

  if(ref $batch) {
    $batch->isa('CXGN::ITAG::Pipeline::Batch')
      or croak 'batch object must be a CXGN::ITAG::Pipeline::Batch';
  } elsif($batch) {
    $batch = $self->pipeline->batch($batch);
  }

  return $batch;
}

=head2 work_dir

  Usage: my $axdir = $an->work_dir($batch);
  Desc : get the full path to the dir this analysis's files reside in
  Args : L<CXGN::ITAG::Pipeline::Batch> object this work dir will be part of
  Ret  : a path string

=cut

memoize 'work_dir';
sub work_dir {
  my ($self,$batch) = @_;
  $batch = $self->_check_batch($batch);
  return File::Spec->catdir($batch->dir,$self->tagname);
}

=head2 cache_dir

  Usage: my $cache_dir = $an->cache_dir($batch)
  Desc : get the (hidden) directory to use for caching
         things relating to this analysis and batch,
         like validation results
  Args : batch object
  Ret  : path to a directory
  Side Effects: creates the directory if it does not already exist

=cut

memoize 'cache_dir';
sub cache_dir {
  my ($self,$batch) = @_;
  my $work_dir = $self->work_dir($batch);
  my $cache_dir = File::Spec->catdir($work_dir,'.itag_pipeline_cache');
  -d $cache_dir || mkdir $cache_dir;
#    or croak "$! creating dir $cache_dir";
  return $cache_dir;
}


=head2 assemble_filename

  Usage: my $fn = $an->assemble_filename($seqname,$desc,$ext);
  Desc : given a sequence name and file type (extension),
         assemble a valid ITAG pipeline filename for it.
         uses assemble_filename() from L<CXGN::ITAG::Tools>
  Args : batch object, sequence name, description string,
        file extension, (optional) file version
  Ret  : a filename string
  Side Effects: none

=cut

sub assemble_filename {
  my ($self,$batch,$seqname,$desc,$ext,$filever) = @_;
  $batch = $self->_check_batch($batch);
  $filever ||= 1;
  $seqname && $desc && $ext or confess 'must provide both a sequence name and extension to assemble_filename()';
  return CXGN::ITAG::Tools::assemble_filename({ seqname      => $seqname,
						analysis_tag => $self->tagname,
						ext          => $ext,
						desc         => $desc,
						pipe_ver     => $self->pipeline->version,
						batch_num    => $batch->batchnum,
						file_ver     => $filever,
						dir          => $self->work_dir($batch),
					      }
					     );
}


=head2 atomic_move

  Usage: $self->atomic_move([$tempfile,$destination],[$tempfile2,$destination2],...);


  Desc : do a series of move() operations quasi-atomically, meaning if
         one of them fails, then the whole batch is rolled back (the targets
         of the ones that succeeded are unlinked if present) and then the
         script dies
  Args : list of arrayref operations as [file => dest],[file => dest],...
         ignores extra elements in each arrayref, if there are more than 2
  Ret  : nothing meaningful
  Side Effects: does some filesystem moves.  will die if some fail
  Example:

=cut

sub atomic_move {
  my ($self,@ops) = @_;

  #move all the analyzed files into position, attempting to roll it
  #back if it failes
  eval {
    foreach (@ops) {
      my ( $source, $target ) = @$_;
      File::NFSLock::uncache( $source );
      sleep 30 unless -f $source; #< wait 30 seconds for NFS to catch
                                  #  up if any of the files are missing
      File::Copy::move( $source, $target )
	    or confess "Failed ($!) mv $source -> $target\n";
    }
  }; if( $EVAL_ERROR ) {
    unlink $_->[1] foreach @ops;
    die $EVAL_ERROR;
  }
}


=head2 local_temp

  Usage: my $tempdir = $self->local_temp('mytemp.seq');

  Desc : first time a package calls it, makes a temporary dir for use
         in storing tempfiles when running this analysis locally.  returns
         that temp dir, and File::Spec->catfile()s any arguments you give it
  Args : (optional) filename inside that temp dir
  Ret  : if no arguments, creates a temp dir (if on
  Side Effects: may create a directory
  Example:

=cut

sub local_temp {
  my ($self,@args) = @_;

  return $self->_cached_tempdir( File::Spec->tmpdir, @args );
}

=head2 cluster_temp

  Usage: my $tempdir = $self->cluster_temp('mytemp.seq');

  Desc : first time a package calls it, makes a CLUSTER-ACCESSIBLE
         temporary dir for use in storing tempfiles when running this
         analysis locally.  returns that temp dir, and
         File::Spec->catfile()s any arguments you give it
  Args : (optional) filename inside that temp dir
  Ret  : if no arguments, creates a temp dir (if on
  Side Effects: may create a directory
  Example:

=cut

{ my $cluster_dir;
  sub cluster_temp {
      my ($self,@args) = @_;
      $cluster_dir ||= CXGN::ITAG::Config->load->{'cluster_shared_tempdir'}
	  or die 'no cluster_shared_tempdir configuration variable defined!';
      return $self->_cached_tempdir( $cluster_dir, @args );
  }
}

#$File::Temp::DEBUG = 1;
sub _cached_tempdir {
  my ($self, $basedir, @args) = @_;
  my $class = ref($self) || $self;

  our %temp_dirs; #< temp dirs cache
  my $dir = $temp_dirs{$basedir}{$class} ||= do {
    my $aname = $self->tagname;
    tempdir( File::Spec->catdir( $basedir, "itag-analysis-$aname-XXXXXX"), CLEANUP => 1 )
      or confess "could not create tempdir\n";
  };

  unless(@args) {
    return $dir;
  } else {
      my $file = File::Spec->catfile($dir,@args);
      pop @args;
      my $dir = File::Spec->catdir($dir,@args);
      -d $dir
          or mkpath($dir)
          or croak "making dir '$dir'";
      return $file;
  }
}


sub cluster_run_class_method {
    my ($self, $batch, $seqname, $method, @args) = @_;
    my $class = ref $self or confess 'cluster_run_class_method is an object method';
    $batch = $self->_check_batch($batch);
    my $run_args = pop @args if ref $args[-1];
    $run_args ||= {};

    $run_args->{temp_base}   ||= $class->cluster_temp;
    $run_args->{working_dir} ||= $class->cluster_temp;

    ref and die "only bare scalars can be passed as args to cluster_run_class_method"
        for @args;

    my $resource_estimate = $self->resource_estimate( $batch, $seqname );
    my $perl_string = "${class}->${method}(".join( ',', map "'$_'",
                                                   @args
                                                 )
                                        .')';

    my %merged_run_args = ( %$resource_estimate, %$run_args );

    return CXGN::Tools::Run->run_cluster
            ( itag_wrapper =>
	      perl =>
              '-M'.$class,
              -e => $perl_string,
              \%merged_run_args,
            );

}


##### PRIVATE METHODS

#lists package names for all available analyses that have defined
#packages.  if a pipeline is passed, creates perl packages for all the
#def files in the pipeline that don't already have them.
sub _an_pkgs {
  my ($pipe) = @_;
  $pipe or confess 'must provide a pipeline object';

  #beginning at the given root package, recursively list all
  #sub-package names that end in all lower-case.
  #recall: %Foo::Bar:: is the symbol table for the package Foo::Bar.
  #it is a hash whose keys are all defined symbols (variables,
  #packages, etc) in a package
  memoize 'list_lc_packages';
  sub list_lc_packages {
    my ($root_packagename) = @_;
    no strict 'refs'; #< using symbolic refs for symbol table names
    my @sub_packages = map {
      s/::$//;
      $root_packagename.'::'.$_
    } grep /::$/,keys %{"${root_packagename}::"};
#    warn "$root_packagename has subs: ".join(',',@sub_packages),"\n";
    my @lc_packs = grep {/::[a-z]\w+$/} @sub_packages;
#    warn "$root_packagename has subs: ".join(',',@sub_packages),"\n";
    return @lc_packs,map {list_lc_packages($_)} @sub_packages;
  }

  my %existing_packages = map {$_ => 1} list_lc_packages(__PACKAGE__);

  #now create packages for all the analyses that seem to be present
  #in our pipeline dir, but have no packages
  foreach my $tag ( grep {!$existing_packages{__PACKAGE__.'::'.$_}}
		    $pipe->list_analyses
		  ) {
    my $newpackage = __PACKAGE__.'::'.$tag;
    my $thispackage = __PACKAGE__;

    no strict 'refs';
    $::{$tag.'::'} = *{'::'.$tag.'::'};
    @{$newpackage.'::ISA'} = ($thispackage);

#     eval 'our @'.$newpackage."::ISA = ('$thispackage')";
#     warn "made new package $newpackage\n";
    $existing_packages{$newpackage} = 1;
  }

  return keys %existing_packages;
}

#look up the package that goes with a given analysis tag
sub _tag2pkg {
  my ($pipe,$tag) = @_;

  #if there's a perl module for it, use that as its package
  if(my @matching_packages = grep { $_ =~ /::$tag$/ } _an_pkgs($pipe)) {
    if(@matching_packages > 1) {
      confess "multiple matches found for analysis tag '$tag':\n",
	map {"  - $_\n"} @matching_packages;
    }
    return $matching_packages[0];
  }
  else {
    die "no analysis found with tag '$tag'\n"
  }
  confess 'this point should never be reached';
}

#get the name of one of one of this analysis's special files, most are
#inside this dir
memoize '_filename';
sub _filename {
  my ($self,$file,$batch) = @_;
  $batch &&= $self->_check_batch($batch);

  my %files = (
	       control  => 'control.txt',
	       manifest => 'manifest.txt',
	       md5sum   => 'md5sums.txt',
	       def      => $self->pipeline->_analysis_def_filename($self->tagname),
	      );

  $files{$file} or croak "unknown file tag '$file'";
  return ($files{$file} =~ m!/!) ? $files{$file} : File::Spec->catfile($self->work_dir($batch),$files{$file});
}

#list of special files found inside of analysis dirs
sub _adir_special_files {
  my ($self,$batch) = @_;
  return map {$self->_filename($_,$batch)} qw/control manifest md5sum/;
}

# =head2 _kv_file_value (private)

#   Usage: my @content = $self->_kv_file_value('keyname');
#   Desc : get the contents of the definition file for the
#          given keyname.  always returns empty if there was a parse error
#   Args : file shortname, like 'def',
#          the name of the key to fetch,
#          optional batch number or object, if required to find the file
#   Ret  : list of values found in the definition file
#   Side Effects: none
#   Example:

# =cut

#memoize '_kv_file_value';
sub _kv_file_value {
  my ($self,$file,$keyname,$batch) = @_;

  my $p = $self->_parse_kv_file($file,$batch)
    or warn "could not parse $file for analysis ".$self->tagname.":\n",map {"- $_\n"} $self->errors;

  return $p && $p->{$keyname} ? @{$p->{$keyname}} : ();
}

# =head2 _parse_kv_file (private)

#   Usage: my $ctl = $self->_parse_kv_file('control');
#   Desc : parse a whitespace-separated key-value file
#   Args : file shortname and batch, fed to _filename(),
#   Ret  : hashref of parsed key-value pairs, as
#        { key => [val,val,val...],
#          ...
#        }, or undef for parse error
#   Side Effects: sets $an->error() on error

# =cut

memoize '_parse_kv_file', NORMALIZER =>
  sub { my ($self,$file,$batch) = @_;
	my $filename = $self->_filename($file,$batch)
	  or die "unknown file type '$file'";
	return $self.$file.(stat($filename))[9].($batch || '');
      };
sub _parse_kv_file {
  my ($self,$file,$batch) = @_;

  my $filename = $self->_filename($file,$batch)
    or die "unknown file type '$file'";

  my $got_errors;
  my $err = sub { $self->_store_error(shift); $got_errors = 1};

  my $p = eval{ CXGN::ITAG::Tools::parse_kv_file($filename) };
  if( $EVAL_ERROR ) {
    $err->($EVAL_ERROR);
    $p = {};
  }

  #check def file
  if( $file eq 'def' ) {
    my %valid_keys = map {$_ => 1}
      qw/owner_group owner_user owner_email depends_on produces_files disabled/;

    foreach my $rkey (qw/owner_user owner_email depends_on produces_files/) {
      defined $p->{$rkey}
	or $err->($self->tagname." def file required key $rkey not found");
    }

    #unknown keys are not allowed in def files
    foreach (keys %$p) {
      $valid_keys{$_}
	or $err->($self->tagname." def file contains unknown key '$_'");
    }

    #check the dependencies list for self-reference
    if( exists $p->{depends_on} ) {
      grep {$_ eq $self->tagname} @{$p->{depends_on}}
	and $err->($self->tagname." def file lists self in depends_on");
    }

    #check the extensions specified in the pipeline def file
    if( exists $p->{produces_files} ) {
      @{$p->{produces_files}} > 0
	or $err->($self->tagname." def file analysis must produce at least one file per input sequence, produces_files must have at least one value");

      my %valid_exts = map {$_=>1} $self->pipeline->valid_file_exts
	or $err->("no valid file extensions defined.  is the pipeline global definitions file present?");
      foreach my $def (@{$p->{produces_files}}) {
	my ($desc,$ext) = split /:/,$def;
	$desc && $ext
	  or $err->($self->tagname." def file: error parsing produces_files line");
	$desc =~ /^[a-z0-9_]+$/
	  or $err->("'$desc' is not a valid file description.  Must contain only lower-case letters, numbers, and underscores");

	$valid_exts{$ext}
	  or $err->("in produces_files value '$def', $ext is not a defined file extension.  If you wish to add it, please add it to the PipelineGeneral page on the ITAG wiki and send an email to the ITAG list");
      }
    }
  }

  #check control file and set defaults
  if( $file eq 'control' ) {
    if( exists $p->{running} ) {
      @{$p->{running}} == 1
	or $err->("in control file, running key can only have one value");

      $p->{running}->[0] == 0 || $p->{running}->[0] == 1
	or $err->("in control file invalid value for running key: '$p->{running}->[0]'");

    } else {
      $p->{running} = [0];
    }
  }

  return if $got_errors;
  return $p;
}


=head1 AUTHOR(S)

Robert Buels

=cut

sub DESTROY {
  my $class = shift;
  our %temp_dirs;
  if( my $dir = $temp_dirs{$class} ) {
    system 'rm', '-rf', $dir;
  }
}


###
1;#do not remove
###
