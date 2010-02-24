package CXGN::ITAG::Config;
use base 'CXGN::Config';
my $defaults = {

    # location of ITAG releases (relative to ftpsite_root)
    itag_release_base        => 'tomato_genome/annotation',

};
sub defaults { shift->SUPER::defaults( $defaults, @_ )}
1;


