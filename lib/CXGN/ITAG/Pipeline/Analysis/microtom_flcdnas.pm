package CXGN::ITAG::Pipeline::Analysis::microtom_flcdnas;
use strict;
use warnings;

use CXGN::Tools::Wget qw/ wget_filter /;

use base 'CXGN::ITAG::Pipeline::Analysis::gth_base';

sub locally_runnable { 1 }

sub _source {
    'MicroTom_FLcDNA'
}

sub _cdna_file {
    my $class = shift;
    wget_filter(
        'ftp://ftp.solgenomics.net/genomes/Solanum_lycopersicum/microtom_flcdna/HTC13227.seq.gz'
            => $class->local_temp('microtom_flcdnas.fasta'),
	{ gunzip => 1 },
	);
}

sub _parse_mode {
    'alignments_merged'
}

1;

