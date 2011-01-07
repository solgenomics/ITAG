package CXGN::ITAG::Tools;
use strict;
use warnings;
use English;
use Carp;

use Memoize;

use File::Spec;
use File::Basename;
use Data::Dumper;

use CXGN::ITAG::Config;
use CXGN::Tools::List qw/all/;

=head1 NAME

CXGN::ITAG::Tools - collection of miscellaneous functions dealing with the
ITAG (International Tomato Annotation Group) distributed annotation
pipeline

=head1 SYNOPSIS

coming soon

=head1 DESCRIPTION

coming soon

=head1 FUNCTIONS

All functions below are EXPORT_OK.

=cut

use base qw/Exporter/;

BEGIN {
  our @EXPORT_OK = qw(
		      parse_filename
		      assemble_filename
		      parse_kv_file
		      parse_release_dirname
		      assemble_release_dirname
		      assemble_release_tag
		      find_published_releases
		     );
}
our @EXPORT_OK;

=head2 parse_filename

  Usage: my $parsed = parse_filename($file);
  Desc : parse an ITAG pipeline-formatted filename into its component parts
  Args : filename to parse
  Ret  : hashref containing the info from the filename,
         { seqname      => versioned sequence name,
           analysis_tag => name of the analysis,
           ext          => file extension, e.g. 'seq',
           desc         => contents description string,
           pipe_ver  	=> pipeline version number,
           batch_num 	=> batch number this file is from,
           file_ver  	=> file version number,
           dir       	=> dirname of the file, if any,
           basename  	=> basename of the file, including ext,
         }
         or undef if the name could not be parsed
  Side Effects: none

=cut

sub parse_filename {
  my ($name,$dirname) = fileparse($_[0]);
  my $orig_name = $name;

  my %parsed;
  foreach my $field (qw/ext file_ver batch_num pipe_ver desc analysis_tag/) {
    #pop the last thing off the end
    $name =~ s/\.([^\.]+)$//;
    $1 or return;
    $parsed{$field} = $1;
  }

  #cut the itag off of pipe_ver and numericize it
  $parsed{pipe_ver} =~ s/^itag//
    or return;
  $parsed{pipe_ver} += 0;

  #cut the batch off the batch num and numericize it
  $parsed{batch_num} =~ s/^batch//
    or return;
  $parsed{batch_num} += 0;

  #cut the v off the file ver and numericize it
  $parsed{file_ver} =~ s/^v//
    or return;
  $parsed{file_ver} += 0;

  $parsed{seqname} = $name;

  $parsed{dir} = $dirname if $dirname;

  #validate
  _valid_parsed(\%parsed)
    or return;

  return \%parsed;
}

sub _valid_parsed {
  my ($parsed) = @_;

  #validate the parsed name
  foreach my $field (qw/seqname desc analysis_tag ext pipe_ver batch_num/) {
    defined($parsed->{$field})
      or return;
  }
  $parsed->{pipe_ver} >= 0 && $parsed->{pipe_ver} =~ /^\d+$/
    or return;

  $parsed->{batch_num} >= 1 && $parsed->{batch_num} =~ /^\d+$/
    or return;

  $parsed->{file_ver} >= 1 && $parsed->{batch_num} =~ /^\d+$/
    or return;

  $parsed->{desc} =~ /^[a-z0-9_]+$/
    or return;

  $parsed->{analysis_tag} =~ /^[a-z][a-z0-9_]+$/
    or return;

  return 1;
}

=head2 assemble_filename

  Usage: my $filename = assemble_filename($parsed_filename);
  Desc : assemble a proper ITAG pipeline filename from info in
         a hashref
  Args : hashref containing the fields in the filename, as:
         { seqname  => versioned sequence name,
           analysis_tag => name of the analysis,
           ext      => file extension, e.g. 'seq',
           desc     => contents description string,
           pipe_ver => pipeline version number,
           batch_num => batch number the file is from,
           file_ver => pipeline version number,
           dir      => dirname of the file, if any
         }
  Ret  : filename string
  Side Effects: dies if the hashref you pass in can't be used

=cut

sub assemble_filename {
  my ($parsed) = @_;

  _valid_parsed($parsed)
    or croak "input must be a hashref, and must contain at least seqname, analysis_tag, desc, pipe_ver, and ext";

  return File::Spec->catfile( $parsed->{dir} ? ($parsed->{dir}) : (),
			      join('.',
				   $parsed->{seqname},
				   $parsed->{analysis_tag},
				   $parsed->{desc},
				   sprintf('itag%03d',$parsed->{pipe_ver}),
				   sprintf('batch%03d',$parsed->{batch_num}),
				   'v'.$parsed->{file_ver},
				   $parsed->{ext},
				  ),
			    );

}

=head2 parse_kv_file

  Usage:
  Desc : parse a key-value format file, similar to the format
         used for CXGN configuration files (see L<CXGN::Configuration>),
         but with a special command, 'include <filename>',
         which allows one to include other similar files.
  Args : filename to parse
  Ret  : hashref representation of the file's contents, as:
         { key_name =>  [ val, val, val, ...],
           ...
         }
  Side Effects: dies if cannot open file, or if circular includes found
  Example:

    a file containing:

      foo  bar baz bloo
      monkeys
      bonobos jump

    would be parsed as

      { foo => ['bar','baz','bloo'],
        monkeys => [],
        bonobos => ['jump'],
      }

=cut

memoize 'parse_kv_file', NORMALIZER =>
  sub { my ($filename,$visited) = @_;
	no warnings 'uninitialized';
	return $filename.','.(stat($filename))[9].','.$visited;
# 	return $nstr;
      };
	
sub parse_kv_file {
  my ($filename,$visited) = @_;
  my $thisdir = dirname($filename);

  CORE::open my $fh,$filename
    or croak "$! opening $filename for reading";

  my $h = {};

  while(my $line = <$fh>) {

    next unless $line =~ /\S/ && $line !~ /^\s*#/;

    my ($key,@vals) = split /\s+/,$line;
    if($key eq 'include') {
      $visited ||= {};
      $visited->{$filename} = 1;
      foreach my $f (@vals) {
	if($visited->{$f}) {
	  croak "error parsing key-value file, circular includes found ($filename,$f)";
	}

	my $i = parse_kv_file(File::Spec->catfile($thisdir,$f),$visited); #< will die if file not found
	%$h = (%$h,%$i); #< add the included values to this one
      }
    } else {
      $h->{$key} = \@vals;
    }
  }

  return $h;
}


=head2 parse_release_dirname

  Usage: my $p = parse_release_dirname($dirname);
  Desc : parse the directory name of an ITAG release directory
  Args : directory name string, either at the end of a path or alone
  Ret  : hashref as:
         {  releasenum => decimal release number, or 0 if a devel release
            pre     => boolean flag of whether this is a prerelease,
            devel   => boolean flag of whether this is a development release,
            dir     => original dir passed in, for your convenience
            basename => basename of dir passed in, for your convenience
         }
  Side Effects: none

=cut

sub parse_release_dirname {
  my ($dir) = @_;

  my $bn = basename($dir);

  my %p;
  @p{'releasenum','devel','pre'} = my @matches = $bn =~ m!^ITAG([^_]*)_(devel_)?(pre_)?release$!i or return;

  return unless $p{devel} || $p{releasenum};

  #convert pre and devel to booleans
  $p{pre}   = !! $p{pre};
  $p{devel} = !! $p{devel};
  $p{releasenum} ||= 0;
  $p{dir} = $dir;
  $p{basename} = $bn;

  return \%p;
}

=head2 assemble_release_dirname

  Usage: my $dirname = assemble_release_dirname({ releasenum => 2, pre => 0, devel => 0});
         #returns 'ITAG2_release'
  Desc : assemble an ITAG release directory name from a hashref
         like that returned by parse_release_dirname
  Args : hashref like that returned by parse_release_dirname() above
  Ret  : string directory name without path
  Side Effects: none

=cut

sub assemble_release_dirname {
  my ($p) = @_;

  my $tag = assemble_release_tag($p)
    or return;

  return $tag.'_release';
}

=head2 assemble_release_tag

  Usage: my $fn = assemble_release_tag({ releasenum=>2})
         #returns 'ITAG2'
  Desc : assemble an ITAG release tag from a hashref
         like that returned by parse_release_dirname
  Args : hashref like that returned by parse_release_dirname() above
  Ret  : release tag string
  Side Effects: none

=cut


sub assemble_release_tag {
  my ($p) = @_;

  my $pre   = $p->{pre}   ? '_pre'    : '';
  if( $p->{devel} ) {
    return "ITAG_devel$pre";
  } else {
    $p->{releasenum} && $p->{releasenum} =~ /^\d+(\.\d+)?$/ && $p->{releasenum} > 0 or croak "releasenum must be a positive integer or decimal, not '$p->{releasenum}'\n";
    return "ITAG$p->{releasenum}${pre}";
  }
}

=head2 find_published_releases

  Usage: my @releases = find_published_releases();
  Desc : finds all ITAG releases that have been published so far,
         returns a list of them, each parsed with parse_release_dirname()
         above
  Args : (optional) root directory to look in.  Defaults
         to tomato_genome/ in the public FTP site
  Ret  : list of parsed directory names, in descending order of modtime
  Side Effects: none

=cut

sub find_published_releases {
  my ($dirname) = @_;
  $dirname ||= do {
    my $ftpsite = CXGN::ITAG::Config->load->{'ftpsite_root'};
    File::Spec->catdir($ftpsite,'tomato_genome');
  };

  return sort { #< sort the results in descending order of modtime, and remember the modtimes
    $b->{modtime} <=> $a->{modtime}
  }
    map { $_->{modtime} = (stat($_->{dir}))[9]; $_ }
    grep $_,
    map parse_release_dirname($_),
      glob File::Spec->catfile($dirname,'*');

}



=head1 AUTHOR(S)

Robert Buels

=cut

###
1;#do not remove
###
