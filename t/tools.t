#!/usr/bin/perl
use strict;
use warnings;
use English;
use Data::Dumper;

use File::Temp qw/tempfile/;
use File::Basename;

BEGIN {
  our %file_tests = (
		     'foo/bar/AC123456.1.repeatmasker.masked.itag002.batch047.v12.fasta' =>
		     { dir     => 'foo/bar/',
		       seqname => 'AC123456.1',
		       desc => 'masked',
		       analysis_tag => 'repeatmasker',
		       pipe_ver => 2,
		       batch_num => 47,
		       file_ver => 12,
		       ext => 'fasta',
		     },
		     'foo/bar/AC123456.1.repeatmasker.itag002.fasta' => undef,
		    );
  our %dir_tests = (
		    'ITAG2.10_release' =>
		    { releasenum => '2.10',
		      pre => '',
		      devel => '',
		      basename => 'ITAG2.10_release',
		    },
		    ITAG3_release =>
		    { releasenum => 3,
		      pre => '',
		      devel => '',
		      basename => 'ITAG3_release',
		    },
		    '/data/shared/ITAG3_release' =>
		    { releasenum => 3,
		      pre => '',
		      devel => '',
		      basename => 'ITAG3_release',
		    },
		    'ITAG2.2_pre_release' =>
		    { releasenum => '2.2',
		      pre => 1,
		      devel => '',
		      basename => 'ITAG2.2_pre_release',
		    },
		    'ITAG3_pre_release' =>
		    { releasenum => 3,
		      pre => 1,
		      devel => '',
		      basename => 'ITAG3_pre_release',
		    },
		    'ITAG_devel_release' =>
		    { releasenum => 0,
		      pre => '',
		      devel => 1,
		      basename => 'ITAG_devel_release',
		    },
		    'ITAG_devel_pre_release' =>
		    { releasenum => 0,
		      pre => 1,
		      devel => 1,
		      basename => 'ITAG_devel_pre_release',
		    },
		   );
}
use Test::More tests => 1 + 2*scalar(keys our %file_tests) + 2 + 2*scalar(keys our %dir_tests);
BEGIN {
  use_ok(  'CXGN::ITAG::Tools', qw/parse_filename assemble_filename parse_kv_file parse_release_dirname assemble_release_dirname find_published_releases/)
    or BAIL_OUT('could not include the module being tested');
}

while(my ($filename,$test) = each our %file_tests) {
  if($test) {
    my $parsed = parse_filename($filename);
    is_deeply( $parsed, $test,
	       "'$filename' parsed correctly" );
    is( assemble_filename($test),
        $filename,
        "'$filename' assembled correctly" );
  } else {
    is(parse_filename($filename),undef,"'$filename' parse fails, as expected")
      or diag Dumper(parse_filename($filename));
    SKIP: { skip 'invalid filename, skip assemble',1 }
  }
}



#test parse_kv_file
my (undef,$t) = tempfile;
open my $tfh, ">$t" or die "$! opening $t";
print $tfh "foo  bar baz bo\n";
close $tfh;

my $h = parse_kv_file($t);
is_deeply($h,{foo => ['bar','baz','bo']}, 'parse_kv_file');

sleep 1; #wait for file caching

open $tfh, ">$t" or die "$! opening $t";
print $tfh "monkey_locs  bushes trees houses\n";
close $tfh;
$h = parse_kv_file($t);

is_deeply($h,{monkey_locs => ['bushes','trees','houses']}, 'parse_kv_file')
  or diag 'actually contained '.Dumper($h);


while(my ($dirname,$test) = each our %dir_tests) {
  $test->{dir} = $dirname;
  my $bn = basename($dirname);
  is( assemble_release_dirname($test), $bn, "assemble_release_dirname works for $dirname" );
  is_deeply( parse_release_dirname($dirname), $test, "parse_release_dirname works for $dirname" );
}

#warn "currently published releases:\n",Dumper find_published_releases();
