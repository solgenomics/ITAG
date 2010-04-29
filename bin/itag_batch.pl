#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

#use Data::Dumper;

use File::Spec;

use Module::Find;

use CXGN::Tools::List qw/str_in/;
use CXGN::Tools::Script;

use CXGN::ITAG::Pipeline;

###### CONFIGURATION

our $default_batch_size = 10_000;
our %global_opt;

###### /CONFIGURATION

#misc constants
use constant GB_VER_ACC_PATTERN => qr/^[A-Z_]+\d+\.\d+$/;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script -d dir [global opts] command [cmd opts] cmd_arg ...

  Script to manage batches in the ITAG distributed annotation
  pipeline.

  Author:  Robert Buels <rmb32\@cornell.edu>

  Commands:

   Run $FindBin::Script help <command> for help on an individual
   command.

    new,create         Create a new batch of sequences
    rm,del             Delete a batch.
    lb                 print a list of batches
    lbs                print the list of sequences in a given batch
    la                 list analyses
    status             query the status of analyses in a batch
    recreate           Check a batch and recreate its analysis
                       directory structure if necessary.
    run                Run local analyses on the batch.
    validate           Run quick validation scan of result files.
    intense_validate   Run computationally intensive validation
                       procedures on result files.

  Global Options:

  -d <dir>
      REQUIRED: base directory where all the ITAG files are kept.

  -v <num>
     pipeline version to operate on.  Defaults to most recent
     pipeline version.

  -L skip script-level locking, run even if another itag_batch.pl
     job is running.  THIS CAN BE DANGEROUS.

  -Q <queue name>
     job queue to use for submitting cluster jobs.  defaults to unset,
     which will use whatever default queue the system's `qsub` is
     currently using

  Examples:

    #create a new batch in the current pipeline version,
    #using the next available batch num
    $FindBin::Script -d /data/shared/tomato_genome/itagpipeline create

    #delete batch number 2 in pipeline version 4
    $FindBin::Script -d /data/shared/tomato_genome/itagpipeline -v 4 del 2

EOU
}
*HELP_MESSAGE = \*usage;

####################################################
# README:  STRUCTURE OF THIS SCRIPT
#
# from this point forware, this script consists of:
# 1.) command definitions and subroutines
# 2.) utility subroutines
# 3.) main code that parses arguments and executes one of
#     the defined commands
#
# TO ADD A COMMAND
#
# 1.) add a call to register_command with the command's
#     name and a help message (see examples below)
# 2.) add a subroutine with the same name as the command.
#     Its @ARGV will contain everything that appears after
#     the command name on the command line
# 3.) That's all!  Just run itag_batch.pl  yourcommand
#     to test it.

####### COMMAND SUBROUTINES

register_command( 'help', <<EOF );
Usage: $FindBin::Script help command_name

   Get help on a command.

EOF
#command to display another command's help message
sub help {
  my ($name,$message) = @_;
  unless($name) {
    ($name) = @ARGV or usage();
  }

  our %help;
  $help{$name} or die "no help available for '$name'";
  print "ERROR: $message\n" if $message;
  print $help{$name};
  exit;
}

my $seq_source_list = join '', map { /::([^:]+)$/; "  - $1\n" }
    Module::Find::findsubmod CXGN::ITAG::SeqSource;

register_command(['new','create'],<<EOF);
Usage: $FindBin::Script [global options] new [options] seq_source_name

Available seq sources are:
$seq_source_list
To see the options available for a given seq source, run
`perldoc CXGN::ITAG::SeqSource::<name>`.

Options:

  -s <size>
    maximum number of sequences to put in a new batch.
    Default $default_batch_size.

  -f include all sequences, regardless of whether they
     are already in another batch

  -C do not include BACs that belong to a contig that
     is in another batch

EOF
sub new {
  my %opts;
  getopts('s:fC',\%opts) or help('new');
  $opts{s} ||= $default_batch_size;

  my $pipe = _open_pipe();

  my ($seq_source_name,@seq_source_args) = @ARGV;
  $seq_source_name or help('new');

  $seq_source_name = 'CXGN::ITAG::SeqSource::'.$seq_source_name;
  eval "require $seq_source_name";
  $@ and die "could not load seq source $seq_source_name:\n$@";

  my $seq_source = $seq_source_name->new(@seq_source_args);

  my $newbatch = $pipe->create_batch( size => $opts{s},
                                      seq_source => $seq_source,
				      $opts{C} ? (no_bacs_in_contigs => 1) : (),
				      $opts{f} ? (include_all => 1) : (),
				    )
    or die "Batch creation failed, aborting.\n";

  print "created new batch number ".$newbatch->batchnum." for pipeline ".$pipe->version."\n";
  print "running `seq` analysis to fetch sequences...";
  $pipe->analysis('seq')->run( $newbatch, seq_source => $seq_source );
  print "done\n";
}

sub _open_pipe {
  return CXGN::ITAG::Pipeline->open( version => $global_opt{v},
				     basedir => $global_opt{d},
				   );
}

register_command(['del','rm'] => <<EOF);
Usage:  $FindBin::Script [global options] del batchnum

  Delete a pipeline batch.  Does not actually delete the files,
  but moves them into a .deleted/ subfolder, from whence you
  could probably recover them if necessary.

EOF
sub del {
  my %opts;
  getopts('',\%opts) or help('del');

  my $batchnum = shift @ARGV;
  my $pipe = _open_pipe;

  my $batch = $pipe->batch($batchnum)
    or die "batch $batchnum does not exist\n";

  print "In pipeline version ".$pipe->version.", are you sure you want to delete batch $batchnum ? (yes/no)\n";
  my $response = <STDIN>;
  unless($response =~ /^\s*yes\s*$/i) {
    print "Aborted.\n";
    return;
  }
  print "batch $batchnum moved to .deleted\n";

  #delete the batch (moves it to <pipedir>/.deleted
  $batch->delete;
}

register_command('recreate' => <<EOF);
Usage:  $FindBin::Script [global opts] recreate <batchnum> <batchnum> ...

  Recreate a pipeline batch, which makes sure that
  everything in a batch is as it should be.  Dirs
  are never deleted, just created if they are not
  there.  Permissions are also set as specified
  in the analysis def files.

  If you pass a batch number of 'all', it will do
  this for all batches in a pipeline.

EOF
sub recreate {
  lock_script() or die "recreate cannot run at the same time as other itag_batch.pl operations\n";

  my %opts;
  getopts('',\%opts) or help('recreate');

  my @batches = @ARGV;
  my $pipe = _open_pipe;

  @batches = $pipe->list_batches if str_in('all',@batches);

  foreach my $batchnum (@batches) {
    my $batch = $pipe->batch($batchnum)
      or die "batch $batchnum does not exist\n";
    $batch->recreate;
  }

  unlock_script();
}


register_command('validate' => <<EOF);
Usage:  $FindBin::Script [global opts] validate <batchnum> [ <analysis_tag> ]

  Run non-intensive validation procedures on the given batch, also
  checking for cached results of intensive validation procedures.  If
  a batch number of 'all' is given, runs validation procedures on all
  batches.

  If an analysis tag is given, runs validation only on results from
  that analysis.

  OPTIONS

   none

EOF
sub validate {
  my ($batchnum,$atag) = @ARGV
    or help('validate');

  my $pipe = _open_pipe;

  my @batchnums = $batchnum eq 'all' ? $pipe->list_batches : ($batchnum);

  my @atags = ($atag) || $pipe->list_analyses;

  #now just go through each analysis in each batch and print its errors
  foreach my $batch (map {$pipe->batch($_)} @batchnums) {
    my @batch_error_strings = map {
      my $analysis = $pipe->analysis($_);
      if(my @errors = $analysis->errors($batch)) {
	' - '.$analysis->tagname.":\n", map "   * $_\n",@errors;
      } else {
	()
      }
    } @atags;
    if(@batch_error_strings) {
      print "batch ".$batch->batchnum.":\n", @batch_error_strings;
    } else {
      print "batch ".$batch->batchnum." contained no errors\n";
    }
  }
}

register_command('intense_validate' => <<EOF);
Usage:  $FindBin::Script [global opts] intense_validate [options] <batchnum> [ <analysis_tag> ]

  Run computationally-intensive validation procedures on the batch.
  If a batch number of 'all' is given, runs validation procedures on
  all batches.

  If an analysis tag is given, runs intense validation only on results
  from that analysis.

  Unless the -f (force) option is given, only runs validation
  procedures on files that are new or have changed since the last
  intensive validation run.

  OPTIONS

  -f   force running validation on all files, regardless of whether
       they are new

  -p <num>
       use up to <num> processors to run validation locally

  -q <queue>
       use the given Torque/PBS/GridEngine queue for running
       validation jobs.

EOF
sub intense_validate {
  lock_script() or die "intense_validate cannot run at the same time as other itag_batch.pl operations\n";

  my %opt;
  getopts('fp:q:',\%opt) or help('intense_validate');
  $opt{p} ||= 1;

  my ($batchnum,$atag) = @ARGV
    or help('intense_validate');

  help('intense_validate','cannot give both -p and -q')
    if $opt{p} && $opt{q};

  my $pipe = _open_pipe;

  my @batchnums = $batchnum eq 'all' ? $pipe->list_batches : ($batchnum);

  #divide the validation by batch and analysis
  my @vjobs = #< array of [batch, analysis], each of which 
    map {
      my $batch = $pipe->batch($_);
      map {
	my $an = $pipe->analysis($_);
	[$batch,$an]
      } $atag || $pipe->list_analyses
    } @batchnums;

  ### now run the validation jobs
  if( $opt{q} || $opt{p} > 1) {
    #^ if we've been given -q or -p args, do some parallelization,
    #  running and controlling other copies of this script 

    my $simultaneous_job_limit = $opt{p} || 1_000_000;

    my @running_jobs;
    do {
      #grep out any jobs that are done
      @running_jobs = grep {$_->alive} @running_jobs;

      #start some more jobs if we need to
      while ( @vjobs && @running_jobs < $simultaneous_job_limit ) {
	my @run_args = do {
	  my %pass_opts;	#< pass along any global options we got
	  $pass_opts{"-$_"} = $global_opt{$_} foreach keys %global_opt;

	  my ($b,$an) = @{pop @vjobs};

	  ( $0,
	    %pass_opts,
	    'intense_validate',
	    $opt{f} ? '-f' : (),
	    $b->batchnum, $an->tagname
	  )
	};

#	use Data::Dumper;
#	warn "using run args ",Dumper(\@run_args);

	push @running_jobs,
	  $opt{q}     ? CXGN::Tools::Run->run_cluster(@run_args, {queue => $opt{q}})
	    : CXGN::Tools::Run->run_async(@run_args);
      }

      sleep 1;			#< wait a sec before we check again
    } while (@running_jobs);

  } else {
    #^ otherwise, just run the validation on each of our jobs in turn
    foreach (@vjobs) {
      my ($batch,$analysis) = @$_;
      $analysis->run_intense_validation($batch,{ force => $opt{f}, job_queue => $global_opt{Q} });
    }
  }

  unlock_script();
}

register_command('run' => <<EOF);
Usage:  $FindBin::Script [global opts] run <batchnum> <analysis_name>

  Run one or more locally runnable analyses on the given pipeline
  batch.  If you give 'all' as the batch name, will run analyses in
  all batches.  If you give 'all' as the analysis name, will run all
  analyses that it can, in the proper sequence according to their
  dependencies.

  So to run absolutely everything possible, you would run
     $FindBin::Script [global opts] run all all

EOF
sub run {
  my %opts;
  getopts('',\%opts) or help('run');

  lock_script() or die "run cannot run at the same time as other itag_batch.pl operations\n";

  my $batchnum = lc shift(@ARGV)
    or help('run');

  my $pipe = _open_pipe;

  my @batches;
  if($batchnum eq 'all') {
    @batches = map $pipe->batch($_),$pipe->list_batches;
  } elsif( $batchnum + 0 == $batchnum ) {
    @batches = ( $pipe->batch($batchnum) or die "batch '$batchnum' not found\n" );
  } else {
    die "invalid batch name, must be a number, or 'all'\n";
  }

  #figure out the analyses to try to run, instantiate them
  my @alist = @ARGV;
  if(@alist == 1 && $alist[0] eq 'all') {
    @alist = grep {
      $_->locally_runnable
    } map {
      $pipe->analysis($_)
	or die "analysis '$_' does not exist\n";
    } $pipe->list_analyses
  } elsif( @alist) {
    @alist = map {
      $pipe->analysis($_)
	or die "analysis '$_' does not exist\n";
    } @alist
  } else {
    die "must specify analysis to run, or 'all'\n";
  }

  #make sure all the specified analyses are locally runnable
  foreach ( @alist ) {
    die $_->tagname." is not locally runnable.  Do you need to make a perl module for it?\n"
      unless $_->locally_runnable;
  }

  #now run all the analyses we can, in ready order
  my %run_records; #< list of things we actually ran
  foreach my $batch (@batches) {
    while( my @ready = grep {$_->uncached_status($batch) eq 'ready'} @alist ) {
      foreach my $r (@ready) {
	warn "running ".$r->tagname." in batch ".$batch->batchnum."\n";
	$r->run($batch);
	$run_records{$batch->batchnum}{$r->tagname} = 1;
      }
    }
  }

  unless($ARGV[0] && $ARGV[0] eq 'all') {
    foreach my $a (@alist) {
      foreach my $batch (@batches) {
	unless($run_records{$batch->batchnum}{$a->tagname}) {
	  warn "WARNING: Did not run analysis '".$a->tagname."' in batch ".$batch->batchnum.", it was never ready.  Current status: '".$a->status($batch)."'\n";
	  if($a->status($batch) eq 'error') {
	    warn map "- $_\n",$a->errors($batch);
	  }
	}
      }
    }
  }

  unlock_script();
}


####### LB #########
register_command( 'lb', <<EOF );
Usage: $FindBin::Script [global options] lb

List currently existing pipeline batch numbers, zero-padded to 3 digits.

EOF
#command to list pipeline batches
sub lb {
  my $pipe = _open_pipe;
  help('lb') if @ARGV;
  print map "$_\n", $pipe->list_batches;
}


register_command('lbs',<<EOF );
Usage: $FindBin::Script [global options] lbs  batch_num batch_num ...
       OR   $FindBin::Script [global options] lbs  all

Print list of sequences identifiers in the given batch.

EOF
sub lbs {
  my $pipe = _open_pipe;
  my @batchnums = @ARGV
    or help('lbs');
  my @all_batches = map $_+0, $pipe->list_batches;
  if( @batchnums == 1 && $batchnums[0] eq 'all' ) {
    @batchnums = @all_batches;
  }

  #validate the given batch numbers
  foreach (@batchnums) {
    str_in($_+0,@all_batches)
      or die "batch '$_' does not exist in this pipeline version\n";
  }

  #now list all the sequences
  foreach my $bn (@batchnums) {
    my $batch = $pipe->batch($bn);
    print map "$_\n", $batch->seqlist;
  }
}

register_command('la',<<EOF);
Usage: $FindBin::Script [global options] la

Outputs a tab-separated list of:

analysis tag \t user:group \t email \t output specs

EOF
sub la {
  my $pipe = _open_pipe;
  foreach my $a ( map $pipe->analysis($_), $pipe->list_analyses ) {
    my $o = $a->owner_info;
    print join "\t", ( $a->tagname,
		       "$o->{user}:$o->{group}",
		       $o->{email},
		       join(',',$a->output_files_specs),
		     );
    print "\n";
  }
}

register_command('status',<<EOF);
Usage: $FindBin::Script [global options] status [options] batch_num optional_analysis_name

List the current status of an analysis or multiple analyses in a batch.

Options:

    -e   print text of any errors present

EOF
sub status {
  my %opts;
  getopts('e',\%opts) or help('status');

  my $pipe = _open_pipe;
  my ($batchnum,$aname) = @ARGV;
  $batchnum or help('status');
  $batchnum =~ /\D/
    and die "batch number must be numeric\n";

  my $batch = $pipe->batch($batchnum)
    or die "batch number $batchnum not found\n";

  my @atags = $aname ? ($aname)
                     : $pipe->list_analyses;
  foreach my $atag (@atags) {
    my $a = $pipe->analysis($atag)
      or die "analysis $atag not found!\n";
    my $stat = $a->status($batch);
    print "$atag\t\t$stat\n";

    if($opts{e} && $stat eq 'error') {
      print map "  -  $_\n", $a->errors($batch);
    }
  }
}



######### UTILITY SUBROUTINES

#register a command with its help text
sub register_command {
  my ($opnames,$helptext) = @_;
  $opnames = [$opnames] unless ref $opnames;
  our %help;
  our %commands;
  foreach my $op (@$opnames) {
    $help{$op}     = $helptext;
    $commands{$op} = \&{"$opnames->[0]"};
  }
}

# make lock and unlock subs that are disabled by -L
sub lock_script {
  return 1 if $global_opt{L};
  CXGN::Tools::Script::lock_script(@_);
}
sub unlock_script {
  return 1 if $global_opt{L};
  CXGN::Tools::Script::unlock_script(@_);
}

############# MAIN CODE

#lock_script() or die "$FindBin::Script already running, only run one at a time.\n";

#parse the global options
getopts('d:v:LQ:',\%global_opt) or usage();

#now execute the command
my $opname = shift @ARGV
  or usage;
our %commands;
my $op = $commands{$opname}
  or die "unknown command '$opname'\n";
unless($opname eq 'help') {
  -d $global_opt{d} or die "pipeline base directory '$global_opt{d}' does not exist\n";
}
$op->();  #< execute the operation subroutine

#unlock_script();
