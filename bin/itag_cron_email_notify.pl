#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use Mail::Sendmail;

use CXGN::IndexedLog;
use CXGN::ITAG::Pipeline;
use CXGN::Tools::List qw/distinct/;

#use Data::Dumper;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script -d itag_dir  email_address

  Scans ITAG pipeline repository for changes, emails the given email
  address about them if it hasn't already.

  Uses .email_log files in the pipeline repository to log what it's
  already emailed about.

  Options:

  -d <dir>
    REQUIRED: base directory to search for ITAG pipelines.

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('d:',\%opt) or usage();
$opt{d} or usage();

my ($itag_email) = @ARGV;
die "must provide an email address\n"
  unless $itag_email;
die "invalid email address '$itag_email'\n"
  unless $itag_email =~ /^\w+@[\w\.]+\w$/;

my $pipe = CXGN::ITAG::Pipeline->open( basedir => $opt{d} );

#open our email log
our $log = CXGN::IndexedLog->open( File => $pipe->email_log_filename );

#now go through all the things we might send an email about, check
#that we haven't already sent an email about it yet, then send a
#single email telling about all the things that happened to the
#pipeline

our @email_descriptions;
our @email_summaries;
our @log_appends;

sub send_if_needed {
  my ($tagwords, $summary, $message, $logmessage) = @_;

  $logmessage ||= $message;
  $logmessage = "$tagwords $logmessage";
  my %entry = $log->lookup( content => $tagwords );

  return if %entry && $entry{content} eq $logmessage;

  push @email_descriptions, $message;
  push @email_summaries, $summary;
  push @log_appends, $logmessage;
}

#send about a new pipeline
send_if_needed('NEW_PIPELINE '.$pipe->version, 'new pipeline version', "NEW PIPELINE VERSION ".$pipe->version." created");

#send about new batches
foreach my $bn ($pipe->list_batches) {
  my $batch = $pipe->batch($bn);
  my $seq_cnt = $batch->seqlist;
  send_if_needed("NEW_BATCH $bn", 'new batch', "NEW BATCH $bn created containing $seq_cnt sequences");

  #send about new analysis state changes
  foreach my $tag ($pipe->list_analyses) {
    my $an = $pipe->analysis($tag);
    my $status = $an->status($batch);
    my $pipever = $pipe->version;
    if( $status eq 'ready' ) {
      send_if_needed("ASTATE $bn:$tag",'analyses ready', "analysis $tag is ready to run in batch $bn","ready");
    }
    elsif( $status eq 'error' ) {
      send_if_needed("ASTATE $bn:$tag",'errors', "analysis $tag has problems in batch $bn, see http://sgn.cornell.edu/sequencing/itag/status_html.pl?pipe=$pipever&batch=$bn for details","error");
    }
  }
}

if( @email_descriptions ) {
  my $action_string = join '', map "* $_\n", @email_descriptions;
  my $email_body = <<EOM;
This is the ITAG pipeline management system, informing you of the following:

$action_string

Regards,

The ITAG Pipeline Managment System
EOM
  sendmail( From => 'itagpipeline@upload.sgn.cornell.edu',
	    To => $itag_email,
	    Subject => 'status - '.join(', ',distinct @email_summaries),
	    Body => $email_body,
	  );
#  warn "sent email:\n".$email_body;
}
# else {
#   warn "no changes detected, no email sent\n";
# }

$log->append($_) foreach @log_appends;
