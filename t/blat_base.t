use strict;
use warnings;
use Test::More;

use File::Temp;
use Path::Class;

use aliased 'CXGN::ITAG::Pipeline::Analysis::blat_base';

my $t1 = temp_pathclass();
my $t2 = temp_pathclass();

# test blat2gff3
blat_base->blat2gff3(
    raw_out  => file(qw( t data blat.subset )),
    gff3_out => $t1,
   );

ok( -s $t1, 'got some gff3 output' );
my $gff3 = $t1->slurp;
is( $gff3, <<'', 'got correct gff3' );
##gff-version 3
Solyc11g011350.1.1	ITAG_blat_base	match	118	532	0.85	+	.	Name=Potato_SGNlocusID_1634_SGN-U287773;Target=Potato_SGNlocusID_1634_SGN-U287773 89 503 +
Solyc11g011090.1.1	ITAG_blat_base	match	3365	3826	0.87	+	.	Name=Potato_SGNlocusID_1634_SGN-U287773;Target=Potato_SGNlocusID_1634_SGN-U287773 81 557 +
Solyc11g011090.1.1	ITAG_blat_base	match	127	476	0.83	+	.	Name=Potato_SGNlocusID_1634_SGN-U287773;Target=Potato_SGNlocusID_1634_SGN-U287773 170 519 +
Solyc11g011080.1.1	ITAG_blat_base	match	197	672	0.89	+	.	Name=Potato_SGNlocusID_1634_SGN-U287773;Target=Potato_SGNlocusID_1634_SGN-U287773 72 547 +
Solyc08g005510.1.1	ITAG_blat_base	match	42	502	0.81	+	.	Name=Potato_SGNlocusID_1634_SGN-U287773;Target=Potato_SGNlocusID_1634_SGN-U287773 73 557 +
Solyc05g024360.1.1	ITAG_blat_base	match	3	178	0.81	-	.	Name=Potato_SGNlocusID_1634_SGN-U287773;Target=Potato_SGNlocusID_1634_SGN-U287773 278 455 +
Solyc05g007850.1.1	ITAG_blat_base	match	1	478	0.87	-	.	Name=Potato_SGNlocusID_1634_SGN-U287773;Target=Potato_SGNlocusID_1634_SGN-U287773 74 557 +
Solyc01g066020.1.1	ITAG_blat_base	match	5	460	0.88	-	.	Name=Potato_SGNlocusID_1634_SGN-U287773;Target=Potato_SGNlocusID_1634_SGN-U287773 81 542 +
Solyc01g008800.1.1	ITAG_blat_base	match	1	471	0.86	+	.	Name=Potato_SGNlocusID_1634_SGN-U287773;Target=Potato_SGNlocusID_1634_SGN-U287773 74 547 +
Solyc11g011350.1.1	ITAG_blat_base	match	118	401	0.89	+	.	Name=Potato_SGNlocusID_1634_SGN-U287774;Target=Potato_SGNlocusID_1634_SGN-U287774 305 593 +
Solyc11g011090.1.1	ITAG_blat_base	match	3365	3644	0.89	-	.	Name=Potato_SGNlocusID_1634_SGN-U287774;Target=Potato_SGNlocusID_1634_SGN-U287774 297 593 +
Solyc11g011080.1.1	ITAG_blat_base	match	197	500	0.93	+	.	Name=Potato_SGNlocusID_1634_SGN-U287774;Target=Potato_SGNlocusID_1634_SGN-U287774 288 593 +
Solyc08g005510.1.1	ITAG_blat_base	match	42	255	0.85	-	.	Name=Potato_SGNlocusID_1634_SGN-U287774;Target=Potato_SGNlocusID_1634_SGN-U287774 289 528 +
Solyc05g024360.1.1	ITAG_blat_base	match	1	99	0.85	+	.	Name=Potato_SGNlocusID_1634_SGN-U287774;Target=Potato_SGNlocusID_1634_SGN-U287774 494 594 +

ok( scalar( blat_base->blat_params ), 'blat_params returned something' );
ok( scalar( blat_base->gff3_source ), 'gff3_source returned something' );

#test run_blat
blat_base->run_blat(
    'C12.5_contig10',
    file(qw( t data sgn_markers.subseqs )),
    file(qw( t data sgn_markers.un_xed_seqs )),
    $t1,
    $t2,
   );

is( $t1->slurp, <<EOT, 'got right psl from blat' );
psLayout version 3

match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts
     	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count
---------------------------------------------------------------------------------------------------------------------------------------------------------------
334	0	0	0	0	0	0	0	+	foo	334	0	334	C12.5_contig10	134255	1380	1714	1	334,	0,	1380,
2199	1	0	0	1	3	1	1	+	bar	2203	0	2203	C12.5_contig10	134255	3540	5741	2	323,1877,	0,326,	3540,3864,
107	0	0	0	1	77	1	131	+	bar	2203	1270	1454	C12.5_contig10	134255	2766	3004	2	57,50,	1270,1404,	2766,2954,
EOT

is( $t2->slurp, <<'', 'translated into right gff3' );
##gff-version 3
C12.5_contig10	ITAG_blat_base	match	1381	1714	1.00	+	.	Name=foo;Target=foo 1 334 +
C12.5_contig10	ITAG_blat_base	match	3541	5741	1.00	+	.	Name=bar;Target=bar 1 2203 +
C12.5_contig10	ITAG_blat_base	match	2767	3004	1.00	+	.	Name=bar;Target=bar 1271 1454 +

done_testing;


# make a tempfile, shove it into a Path::Class::File so it does not
# get destroyed, and return the Path::Class::File
sub temp_pathclass {
    my $tmp = File::Temp->new;
    $tmp->close;
    my $t = file("$tmp");
    $t->{myhackeduptemp} = $tmp;
    return $t;
}
