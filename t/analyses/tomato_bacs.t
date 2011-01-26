use strict;
use warnings;

use Test::More;
use Path::Class;

use aliased 'CXGN::ITAG::Pipeline::Analysis::tomato_bacs';

ok( tomato_bacs->gff3_source, 'gff3_source returns something' );

my $accessions_file = tomato_bacs->_accessions_file;
ok( -f $accessions_file, 'able to download bacs accessions file' );

my %aliases = tomato_bacs->_all_aliases;

is_deeply( tomato_bacs->_get_aliases('nonexistent'),
           [],
           '_get_aliases for something nonexistent returns empty array' );


my ($test_key) = keys %aliases;

is_deeply( tomato_bacs->_get_aliases( $test_key ),
           $aliases{$test_key},
           '_get_aliases gets actual data' );

cmp_ok( scalar( @{tomato_bacs->_get_aliases( $test_key )} ),
        '>=',
        2,
        "alias for $test_key has at least 2 elements in it"
       );
#diag "aliases for $test_key are: ", explain tomato_bacs->_get_aliases( $test_key );

tomato_bacs->_load_deflines( file(qw( t data sgn_markers.subseqs ) ) );

is( tomato_bacs->_get_defline( 'foo' ),
    'perfect match',
    'got a right defline',
   );

done_testing;
