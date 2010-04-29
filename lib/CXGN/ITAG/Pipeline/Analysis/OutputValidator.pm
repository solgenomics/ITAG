### superclass for all the different types of output validation
package CXGN::ITAG::Pipeline::Analysis::OutputValidator;
use strict;
use warnings;

use Carp;
use Storable ();
use File::Basename;
use Memoize;
use List::MoreUtils qw/any/;

sub readable_name {
  shift->name;
}

sub name { #< the validator name of this class is given by its package name
  my ($self) = @_;
  my ($name) = (ref($self) || $self)  =~ /::([^:]+)$/;
#  $name =~ s/_/ /g;
  return $name;
}

sub is_intensive { 1 } #< returns true if this should be run offline

memoize 'cache_dir';
sub cache_dir {
  my ($self,$analysis,$batch) = @_;
  my $cache_dir = File::Spec->catdir( $analysis->cache_dir($batch), $self->name );
  -d $cache_dir || mkdir $cache_dir;
#    or croak "$! creating dir $cache_dir";
  return $cache_dir;
}
memoize 'failures_filename';
sub failures_filename {
  my ($self,$analysis,$batch) = @_;
  return File::Spec->catfile($self->cache_dir($analysis,$batch),
			     'failures.dat'
			    );
}
# list the result filename and time of each failed file
sub failures {
  my $file = shift->failures_filename(@_);
  return {} unless -f $file;
  return Storable::retrieve( $file );
}
sub write_failures {
  my ($self,$analysis,$batch,$failures) = @_;
  my $file = $self->failures_filename($analysis,$batch);
  if( %$failures ) {
    Storable::nstore( $failures, $file );
  }
  elsif( -f $file ) {
    unlink $file or confess "$! unlinking $file";
  }
}

memoize 'report_filename';
sub report_filename {
  my ($self,$analysis,$batch,$filename) = @_;
  my $bn = basename($filename);
  return File::Spec->catfile($self->cache_dir($analysis,$batch),
			     $bn.'.report');
}

sub needs_update {
  my ($self,$analysis,$batch) = @_;

  my @files = $self->files_to_validate($analysis,$batch);
  return @files && any {
      my $file = $_;
      my $reportfile = $self->report_filename($analysis,$batch,$file);
      my $rs = (stat($reportfile))[9];
      !$rs || $rs < (stat($file))[9];
  } @files
}

sub files_to_validate {
  my ($self,$analysis,$batch) = @_;
  $analysis or confess 'must pass analysis object to files_to_validate()';
  $batch or confess 'must pass batch object to files_to_validate()';
  return grep { $self->validates_file($_) && -f } $analysis->expected_output_files($batch);
}

#return true if this validator should be run on this file
sub validates_file { 0 };

sub errors { #< return an array of errors from the cache, or from the
             #run_online() method if this validation is not intensive
  my ($self,$analysis,$batch) = @_;
  $analysis or croak 'must provide analysis object';
  $batch or croak 'must provide batch object';

  if( $self->is_intensive() ) {
    my @errors;

    my $failures = $self->failures($analysis,$batch);
    foreach my $file ($analysis->expected_output_files($batch)) {
      if ( $failures->{$file} && $failures->{$file} >= (stat($file))[9]) {
	my $bn = basename($file);
	if( -f $self->report_filename($analysis,$batch,$file) ) {
	  push @errors, $self->readable_name." validation failed for file $bn, see [report file]";
	} else {
	  push @errors, $self->readable_name." validation failed for file $bn";
	}
      }
    }
    return @errors;
  }
  else {
    return $self->run_online($analysis,$batch);
  }
}

#default run_offline uses the run_online() method to run the analysis
sub run_offline {
  my ($self,$analysis,$batch,$options) = @_;

  my @errors = $self->run_online($analysis,$batch);
  #use Data::Dumper;
  #warn "got errors ".Dumper(\@errors);

  #hash the errors by basename
  my %errors;
  foreach my $error (@errors) {
    my ($file,$error) = split /\s*:\s*/,$error,2;
    push @{$errors{$file}},$error;
  }

  my @files = $self->files_to_validate($analysis,$batch)
    or return; #< don't do any of this unless we have some files

  my %failures;
  my $failtime = time;
  foreach my $file (@files) {
    my $bn = basename($file);
    #warn "validating $bn...\n";
    my $outbn = $self->cache_dir($analysis,$batch)."/$bn";
    my $reportfile = $self->report_filename($analysis,$batch,$file);
    #next if -f $reportfile && (stat($reportfile))[9] >= (stat($file))[9]; #< skip if the report is new enough

    CORE::open(my $er,">$reportfile") or die "failed to open validation failure report file $reportfile ($!)";
    if($errors{$bn}) {
      #record this failure and time in the hash
      $failures{$file} = $failtime;
      #print a report
      print $er "$bn: $_\n" for @{$errors{$file}};
    }
    close $er;
  }

  $self->write_failures( $analysis, $batch, \%failures );
}

sub run_online  {}



1;
