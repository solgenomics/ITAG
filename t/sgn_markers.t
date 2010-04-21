use strict;
use warnings;

use FindBin;
use File::Temp;
use File::Spec;

use IPC::Cmd qw/ can_run /;

use Test::More tests => 13;

use CXGN::ITAG::Pipeline;

my $class = 'CXGN::ITAG::Pipeline::Analysis::sgn_markers';
use_ok $class;

my ( $sgn_markers_gthxml, $un_xed_seqs ) =
    map File::Spec->catfile( $FindBin::RealBin, 'data', $_ ),
    qw/ sgn_markers.gthxml  sgn_markers.un_xed_seqs /;

my $temp_out = File::Temp->new;
$temp_out->close;

my $seqname = 'C00.6_contig21';


# test _gthxml_to_gff3
$class->_gthxml_to_gff3( $sgn_markers_gthxml, $un_xed_seqs, $seqname, $temp_out->filename );
my @correct_features =
    (
     { Name => 'C2_At4g00560',
       Alias => 'SGN-M8140',
     },
    );
my $gff3_in = Bio::FeatureIO->new( -format => 'gff', -version => 3, -file => $temp_out->filename );
while( my $f = $gff3_in->next_feature ) {
    my $cf = shift @correct_features;
    isa_ok( $f, 'Bio::AnnotatableI' );
    while( my ($k,$v) = each %$cf ) {
        is( ($f->annotation->get_Annotations($k))[0]->value, $v );
    }
}

SKIP: {
    skip 'set INTENSIVE_TESTS environment variable to run long-running, cpu-intensive tests', 1
        unless $ENV{INTENSIVE_TESTS};
    skip 'genomethreader not installed', 2
        unless can_run( 'gth' );

    my $xml_out = File::Temp->new;
    $class->run_gth( $seqname,
                     $un_xed_seqs,
                     $xml_out->filename,
                     $temp_out->filename,
                    );

    my $out_gff3 = slurp( $temp_out );
    like( $out_gff3, qr/Name=C2_At4g/, 'report looks good' );

    my $out_xml = slurp( $xml_out );
    like( $out_xml, qr/no longer provided/, 'XML is stubbed out' );

#    print $out_gff3;
}



sub slurp {
    open( my $f, '<', +shift ) or die;
    local $/;
    return <$f>
}
