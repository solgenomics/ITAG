package CXGN::ITAG::Pipeline::Analysis::potato_bacs;
use strict;
use warnings;

use base qw/CXGN::ITAG::Pipeline::Analysis::tomato_bacs/;

use CXGN::Tools::Wget 'wget_filter';

sub gff3_source {
    'ITAG_potato_bacs'
}

sub _query_file_url {
    'ftp://ftp.solgenomics.net/genomes/Solanum_tuberosum/bacs/curr/all_clone_sequence.fasta.gz'
}

sub _accessions_file {
    my ( $class ) = @_;
    return wget_filter(
        'ftp://ftp.solgenomics.net/genomes/Solanum_tuberosum/bacs/curr/genbank_accessions.txt',
            => $class->local_temp('accession_list.txt'),
       );
}

1;
