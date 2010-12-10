package CXGN::ITAG::Pipeline::Analysis::itag1_ref;
use strict;
use warnings;
use base qw/CXGN::ITAG::Pipeline::Analysis::mummer_base/;

sub gff3_source {
    'ITAG_SL1.00_scaffolds'
}

sub query_file_url {
    'cxgn-resource://genome_SL1.00_scaffolds'
}

1;
