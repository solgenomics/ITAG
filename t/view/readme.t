use strict;
use warnings;

use Test::More;

use CXGN::ITAG::Release;
use CXGN::ITAG::Release::Statistics;

use Mock::Quick;

use_ok 'CXGN::ITAG::Release::View::README';

{
    my $readme = CXGN::ITAG::Release::View::README->new( release => 'fake' );
    my $bins = $readme->_bins_with_cutoff(scalar( qobj( min => 0, max => 400_000 )), 3, 20_000 );
    is_deeply( $bins, [10_000,20_000,400_000] ) or diag explain $bins;
}

SKIP: {

    skip 'set ITAG_TEST_RELEASE to "name:dir" to set a release to test against',1
        unless $ENV{ITAG_TEST_RELEASE};

    my ( $num, $dir ) = split /:/, $ENV{ITAG_TEST_RELEASE};

    my $rel = CXGN::ITAG::Release->new( releasenum => $num, dir => $dir );
    my $readme = CXGN::ITAG::Release::View::README->new( release => $rel );

    my $txt = $readme->render;
    diag $txt;

}

done_testing;
