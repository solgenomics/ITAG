package CXGN::ITAG::Pipeline::Analysis::potato_bacs;
use strict;
use warnings;

use base qw/CXGN::ITAG::Pipeline::Analysis::tomato_bacs/;

use CXGN::Tools::Wget 'wget_filter';

sub gff3_source {
    'ITAG_potato_bacs'
}

sub query_file {
    my $class = shift;
    return wget_filter(
        'ftp://ftp.solgenomics.net/genomes/Solanum_tuberosum/bacs/curr/bacs.seq',
            => $class->cluster_temp('potato_bacs.fasta'),
       );
}

sub _accessions_file {
    my ( $class ) = @_;
    return wget_filter(
        'ftp://ftp.solgenomics.net/genomes/Solanum_tuberosum/bacs/curr/genbank_accessions.txt',
            => $class->local_temp('accession_list.txt'),
       );
}

1;
