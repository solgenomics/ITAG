package CXGN::ITAG::Config;
use base 'CXGN::Config';
my $defaults = {

    # itag pipeline variables
    itag_pipeline_base       => '/data/shared/tomato_genome/itagpipeline/itag',

    # location of ITAG releases (relative to ftpsite_root)
    itag_release_base        => 'tomato_genome/annotation',

    itag_loading_log_table   => 'itag_loading_log',

};
sub defaults { shift->SUPER::defaults( $defaults, @_ )}
1;


