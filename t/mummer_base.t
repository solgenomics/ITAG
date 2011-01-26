use strict;
use warnings;
use Test::Most;

use IO::Scalar;

use File::Temp;
use Path::Class;

use String::Random;

use aliased 'CXGN::ITAG::Pipeline::Analysis::mummer_base';

{
    my $mummer_text = <<'';
> SL2.30ch00  Len = 21839854
  C03HBa0060B10.1     71156  16474710      2763
> SL2.30ch00 Reverse  Len = 21839854
  C03SLe0022L21.1     17275   4332044      1171
  C03SLm0006J07.1     89937   6074107      1108
  C03SLm0006J07.1     38337  11594424      3240
  C03SLm0006J07.1     43100  11598967      1865
  C03SLm0006J07.1     45616  11601560      3324
  C03SLm0006J07.1     51129  11606240      5319
  C03SLm0006J07.1     57850  11612414      1139
  C03SLm0006J07.1     59867  11614703      6960

    my $mummer_results =
        [
            {
                'match_len'      => '2763',
                'q_seq_len'      => '21839854',
                'q_start'        => '71156',
                'query'          => 'SL2.30ch00',
                'query_reversed' => 0,
                's_start'        => '16474710',
                'subject'        => 'C03HBa0060B10.1'
            },
            {
                'match_len'      => '1171',
                'q_seq_len'      => '21839854',
                'q_start'        => '17275',
                'query'          => 'SL2.30ch00',
                'query_reversed' => 1,
                's_start'        => '4332044',
                'subject'        => 'C03SLe0022L21.1'
            },
            {
                'match_len'      => '1108',
                'q_seq_len'      => '21839854',
                'q_start'        => '89937',
                'query'          => 'SL2.30ch00',
                'query_reversed' => 1,
                's_start'        => '6074107',
                'subject'        => 'C03SLm0006J07.1'
            },
            {
                'match_len'      => '3240',
                'q_seq_len'      => '21839854',
                'q_start'        => '38337',
                'query'          => 'SL2.30ch00',
                'query_reversed' => 1,
                's_start'        => '11594424',
                'subject'        => 'C03SLm0006J07.1'
            },
            {
                'match_len'      => '1865',
                'q_seq_len'      => '21839854',
                'q_start'        => '43100',
                'query'          => 'SL2.30ch00',
                'query_reversed' => 1,
                's_start'        => '11598967',
                'subject'        => 'C03SLm0006J07.1'
            },
            {
                'match_len'      => '3324',
                'q_seq_len'      => '21839854',
                'q_start'        => '45616',
                'query'          => 'SL2.30ch00',
                'query_reversed' => 1,
                's_start'        => '11601560',
                'subject'        => 'C03SLm0006J07.1'
            },
            {
                'match_len'      => '5319',
                'q_seq_len'      => '21839854',
                'q_start'        => '51129',
                'query'          => 'SL2.30ch00',
                'query_reversed' => 1,
                's_start'        => '11606240',
                'subject'        => 'C03SLm0006J07.1'
            },
            {
                'match_len'      => '1139',
                'q_seq_len'      => '21839854',
                'q_start'        => '57850',
                'query'          => 'SL2.30ch00',
                'query_reversed' => 1,
                's_start'        => '11612414',
                'subject'        => 'C03SLm0006J07.1'
            },
            {
                'match_len'      => '6960',
                'q_seq_len'      => '21839854',
                'q_start'        => '59867',
                'query'          => 'SL2.30ch00',
                'query_reversed' => 1,
                's_start'        => '11614703',
                'subject'        => 'C03SLm0006J07.1'
            }
        ];

    {
        my $f = File::Temp->new;
        $f->print( $mummer_text );
        $f->close;
        my $r = mummer_base->_parse_mummer( $f );
        is_deeply( $r,
                   $mummer_results,
                   'mummer parse works',
                  )
            or diag explain $r;
    }


    my $test_gff3 = '';
    mummer_base->_mummer_to_gff3(
        $mummer_results,
        IO::Scalar->new( \$test_gff3 ),
       );

    eq_or_diff( $test_gff3, <<"", 'got right test gff3' );
##gff-version 3
##sequence-region SL2.30ch00 1 21839854
SL2.30ch00	ITAG_genome	match	71156	73918	.	+	.	Name=C03HBa0060B10.1;Target=C03HBa0060B10.1 16474710 16477472 +
SL2.30ch00	ITAG_genome	match	21821410	21822580	.	-	.	Name=C03SLe0022L21.1;Target=C03SLe0022L21.1 4332044 4333214 +
SL2.30ch00	ITAG_genome	match	21748811	21749918	.	-	.	Name=C03SLm0006J07.1;Target=C03SLm0006J07.1 6074107 6075214 +
SL2.30ch00	ITAG_genome	match	21798279	21801518	.	-	.	Name=C03SLm0006J07.1;Target=C03SLm0006J07.1 11594424 11597663 +
SL2.30ch00	ITAG_genome	match	21794891	21796755	.	-	.	Name=C03SLm0006J07.1;Target=C03SLm0006J07.1 11598967 11600831 +
SL2.30ch00	ITAG_genome	match	21790916	21794239	.	-	.	Name=C03SLm0006J07.1;Target=C03SLm0006J07.1 11601560 11604883 +
SL2.30ch00	ITAG_genome	match	21783408	21788726	.	-	.	Name=C03SLm0006J07.1;Target=C03SLm0006J07.1 11606240 11611558 +
SL2.30ch00	ITAG_genome	match	21780867	21782005	.	-	.	Name=C03SLm0006J07.1;Target=C03SLm0006J07.1 11612414 11613552 +
SL2.30ch00	ITAG_genome	match	21773029	21779988	.	-	.	Name=C03SLm0006J07.1;Target=C03SLm0006J07.1 11614703 11621662 +


    my $random_dna = random_dna(5000);
    my $genome_fasta = tempfile_containing( ">FooGenome\n$random_dna\n" );

    my $seq_file = tempfile_containing( ">BarSequence\n".substr( $random_dna, 123, 4222 )."\n" );

    SKIP : {
        eval {
            mummer_base->run_mummer( 'BarSequence', $genome_fasta, $seq_file );
        };
        skip "mummer not present", 2 if $@;

        my $job = mummer_base->run_mummer( 'BarSequence', $genome_fasta, $seq_file );
        is( $job->out, <<'', 'right mummer out' );
> BarSequence  Len = 4222
  FooGenome       124         1      4222
> BarSequence Reverse  Len = 4222

        my $gff3;
        mummer_base->_mummer_to_gff3( mummer_base->_parse_mummer( $job->out_file ), IO::Scalar->new( \$gff3 ) );


        is( $gff3, <<'', 'right mummer gff3' );
##gff-version 3
##sequence-region BarSequence 1 4222
BarSequence	ITAG_genome	match	124	4345	.	+	.	Name=FooGenome;Target=FooGenome 1 4222 +

    }

}

done_testing;


sub tempfile_containing {
    my ($contents) = shift;

    my $f = File::Temp->new;
    $f->print($contents);
    $f->close;
    return $f;
}

sub random_dna {
    my ($length) = @_;
    my $r = String::Random->new;
    $r->{D} = [qw[ A C T G ]];
    return $r->randpattern('D'x$length);
}
