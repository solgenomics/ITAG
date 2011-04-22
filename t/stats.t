use strict;
use warnings;

use Test::More;

use CXGN::ITAG::Release;
use CXGN::ITAG::Release::Statistics;


unless( $ENV{ITAG_TEST_RELEASE} ) {
    plan skip_all => 'set ITAG_TEST_RELEASE to "name:dir" to set a release to test against';
}

my ( $num, $dir ) = split /:/, $ENV{ITAG_TEST_RELEASE};

my $rel = CXGN::ITAG::Release->new( releasenum => $num, dir => $dir );
my $stats = CXGN::ITAG::Release::Statistics->new( release => $rel );

cmp_ok( $stats->gene_model_length->count, '>', 0, 'got some genes in the stats' );


done_testing;
