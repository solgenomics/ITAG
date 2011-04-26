use strict;
use warnings;

use Test::More;

use CXGN::ITAG::Release;
use CXGN::ITAG::Release::Statistics;

unless( $ENV{ITAG_TEST_RELEASE} ) {
    plan skip_all => 'set ITAG_TEST_RELEASE to "name:dir" to set a release to test against';
}

use_ok 'CXGN::ITAG::Release::View::README';

my ( $num, $dir ) = split /:/, $ENV{ITAG_TEST_RELEASE};

my $rel = CXGN::ITAG::Release->new( releasenum => $num, dir => $dir );

my $readme = CXGN::ITAG::Release::View::README->new( release => $rel );
my $txt = $readme->render;
diag $txt;


done_testing;
