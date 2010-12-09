package CXGN::ITAG::Pipeline::Analysis::tobacco_contigs;
use strict;
use warnings;

use base qw/CXGN::ITAG::Pipeline::Analysis::blat_base/;

use CXGN::Tools::Wget 'wget_filter';

sub gff3_source {
    'ITAG_tobacco_contigs'
}

sub query_file {
    my $class = shift;
    return wget_filter(
        'ftp://ftp.solgenomics.net/genomes/Nicotiana_tabacum/assembly/curr/tobacco_genome_sequences_assembly.fasta'
            => $class->cluster_temp('tobacco_genome.fasta'),
       );
}

1;
