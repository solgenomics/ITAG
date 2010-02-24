#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

#use Data::Dumper;

use CXGN::DB::GFF::Versioned;

use CXGN::ITAG::Release;
use CXGN::Tools::Script qw/lock_script unlock_script/;

use CXGN::Config;
use CXGN::IndexedLog;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script

  Reload whatever ITAG GBrowse databases need reloading.  Looks at
  what ITAG releases exist on the FTP site, and syncs our complement
  of GBrowse DBs to reflect them.

  Options:

    none yet

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('',\%opt) or usage();

lock_script() or die "Don't run more than one $FindBin::Script at once!\n";

# find all ITAG release directories
#my @releases = find_published_releases();
my @releases = CXGN::ITAG::Release->find();

unless(@releases) {
  warn "no releases found, exiting.\n";
  exit;
}

# if they are regular or pre, make sure they are loaded
foreach my $one_time_release (grep !$_->is_devel_release, @releases) {
  my $bn = $one_time_release->dir_basename;
  unless( $one_time_release->has_gbrowse_dbs ) {
    warn "versioned database $bn does not exist, loading...\n";
    load_release($one_time_release);
    warn "finished loading $bn\n";
  }
  else {
    warn "versioned database $bn already exists\n";
  }
}

# if they are devel, check the timestamps of the files and do a versioned reload if necessary
if(my ($devel_release) = grep $_->is_devel_release, @releases) {
  #open our processing log
    my $conf = CXGN::Config->load_locked;
    my $dbh = CXGN::DB::Connection->new({ config => $conf });
    my $log = CXGN::IndexedLog->open( DB => $dbh, 'itag_loading_log' );

    my %f = get_fileset($devel_release);
    my $modtime = (stat( $f{genomic_gff3}))[9];

    my $log_string = "LOADED_DEVEL_RELEASE ".$devel_release->dir.",$modtime";

    unless( $log->lookup(content => $log_string) ) {
	load_release($devel_release);
	$log->append($log_string);
    }
}

unlock_script();

########### SUBROUTINES #######

sub load_release {
  my ($release) = @_;

  my %f = get_fileset($release);

  my $db_name_root = $release->release_name;

  #load each of the three gbrowse DBs: genomic and protein
  foreach my $type ('genomic','protein') {
    warn "loading $type for $db_name_root...\n";
    my $bdb = CXGN::DB::GFF::Versioned->new( -db => $release->dbname_root($type), -user => '' );
    $bdb->load_new($f{$type.'_seq'},$f{$type.'_gff3'});
  }
}

sub get_fileset {
  my $release = shift;
  my %files = qw(genomic_gff3  combi_genomic_gff3
		 genomic_seq   genomic_fasta
		 protein_gff3     functional_prot_gff3
		 protein_seq      protein_fasta
		);
  #lookup the filename for each of these
  $_ = $release->get_file_info($_)->{file} foreach values %files;

  while( my ($name,$file) = each %files) {
    -f $file or warn "no $name found in release dir ".$release->dir;
  }

  while( my ($name,$file) = each %files) {
    -f $file or die "no $name found in release dir ".$release->dir;
  }

  return %files;
}


