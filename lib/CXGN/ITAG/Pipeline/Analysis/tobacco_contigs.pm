package CXGN::ITAG::Pipeline::Analysis::tobacco_contigs;
use strict;
use warnings;

use base qw/CXGN::ITAG::Pipeline::Analysis::blat_base/;

use CXGN::Tools::Wget 'wget_filter';

sub gff3_source {
    'ITAG_tobacco_contigs'
}

sub query_file_url {
    'ftp://ftp.solgenomics.net/genomes/Nicotiana_tabacum/assembly/curr/tobacco_genome_sequences_assembly.fasta'
}

1;
