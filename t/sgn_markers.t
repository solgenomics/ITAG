use strict;
use warnings;

use FindBin;
use File::Temp;
use File::Spec;

use Test::More tests => 11;

my $class = 'CXGN::ITAG::Pipeline::Analysis::sgn_markers';
use_ok $class;

my ( $sgn_markers_gthxml, $un_xed_seqs ) =
    map File::Spec->catfile( $FindBin::RealBin, 'data', $_ ),
    qw/ sgn_markers.gthxml  sgn_markers.un_xed_seqs /;

my $temp_out = File::Temp->new;
$temp_out->close;

$class->_gthxml_to_gff3( $sgn_markers_gthxml, $un_xed_seqs, 'C00.6_contig21', $temp_out->filename );

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
