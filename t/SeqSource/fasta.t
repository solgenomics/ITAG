#!/usr/bin/perl
use strict;
use warnings;
use English;

use FindBin;
use File::Spec;

use Data::Dumper;

use Test::More tests => 14;

BEGIN {
  use_ok(  'CXGN::ITAG::SeqSource::Fasta'  )
    or BAIL_OUT('could not include the module being tested');
}

my $source = CXGN::ITAG::SeqSource::Fasta->new( file => File::Spec->catfile( $FindBin::RealBin, 'data', 'foo.seq' ) );

#now get the sequences for all of them
my @names = $source->all_seq_names;
is( scalar(@names), 4, 'correct name count' );


my %seqs = qw[
foo1
ACTGACACGATCGACTAGCTACATCATCAGGCGGATCGA
foo2
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
foo3
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
foo4
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
];

my %contigs = (
               foo1 => [
                        {
                         'objname' => 'foo1',
                         'ostart' => 1,
                         'cend' => 39,
                         'oend' => 39,
                         'ident' => 'foo1',
                         'typedesc' => 'finished',
                         'orient' => '+',
                         'linenum' => 1,
                         'type' => 'F',
                         'cstart' => 1,
                         'partnum' => 1
                        }
                       ],
               foo2 => [
                        {
                         'objname' => 'foo2',
                         'ostart' => 1,
                         'cend' => 50,
                         'oend' => 50,
                         'ident' => 'foo2',
                         'typedesc' => 'finished',
                         'orient' => '+',
                         'linenum' => 1,
                         'type' => 'F',
                         'cstart' => 1,
                         'partnum' => 1
                        }
                       ],
               foo3 => [
                        {
                         'objname' => 'foo3',
                         'ostart' => 1,
                         'cend' => 50,
                         'oend' => 50,
                         'ident' => 'foo3',
                         'typedesc' => 'finished',
                         'orient' => '+',
                         'linenum' => 1,
                         'type' => 'F',
                         'cstart' => 1,
                         'partnum' => 1
                        }
                       ],
               foo4 => [
                        {
                         'objname' => 'foo4',
                         'ostart' => 1,
                         'cend' => 57,
                         'oend' => 57,
                         'ident' => 'foo4',
                         'typedesc' => 'finished',
                         'orient' => '+',
                         'linenum' => 1,
                         'type' => 'F',
                         'cstart' => 1,
                         'partnum' => 1
                        }
                       ],
              );


foreach my $name ($source->all_seq_names) {
    my ($seq, $contig) = $source->get_seq( $name );
    is( $seq->id, $name );
    is( $seq->seq, $seqs{$name} );
    is_deeply( $contig, $contigs{$name}, 'got correct contig contents')
        or diag "got contig:\n".Dumper $contig;
}
