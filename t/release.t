#!/usr/bin/perl
use strict;
use warnings;
use English;

use Test::More tests => 31;

use File::Temp qw/tempdir/;

BEGIN {
  use_ok(  'CXGN::ITAG::Release'  )
    or BAIL_OUT('could not include the module being tested');
}

my $tempdir = tempdir(File::Spec->catdir(File::Spec->tmpdir,'cxgn-itag-release-test-XXXXXX'),CLEANUP => 1);

mkdir "$tempdir/ITAG1_release";
mkdir "$tempdir/ITAG1_pre_release";
mkdir "$tempdir/ITAG_devel_release";
symlink './ITAG_devel_release', "$tempdir/ITAG_current";

#system 'ls', '-l', $tempdir;

CXGN::ITAG::Release->releases_root_dir($tempdir);
my @rels = sort {$a->release_tag cmp $b->release_tag} CXGN::ITAG::Release->find();
ok(@rels == 3,'found 2 releases');
my @only_dev = CXGN::ITAG::Release->find(devel => 1);
my @only_curr = CXGN::ITAG::Release->find(current => 1);
#print $_->dir."\n" foreach @only_dev;
is(scalar(@only_dev),1,'devel search found 1');
is(scalar(@only_curr),1,'curr search found 1');
is($only_dev[0]->release_tag,'ITAG_devel','and it is a devel');
ok( $only_dev[0]->is_current_release, 'this sole devel release should read as the current release' );
is( $only_curr[0]->dir, $only_dev[0]->dir, 'is_current_release search works');


ok(! $rels[0]->is_devel_release,'it is not a devel');
ok(! $rels[0]->is_pre_release,'it is not a pre');

# test get_file_info
eval{ $rels[0]->get_file_info('combi_gff3') };
like( $EVAL_ERROR, qr/unknown/, 'get_file_info dies on unknown shortname');
is($rels[0]->get_file_info('combi_genomic_gff3')->{file},"$tempdir/ITAG1_release/ITAG1_all_genomic.gff3",'got correct combi_genomic_gff3 filename');
is($rels[1]->get_file_info('combi_genomic_gff3')->{file},"$tempdir/ITAG1_pre_release/ITAG1_pre_all_genomic.gff3",'got correct combi_genomic_gff3 filename');

# test release_tag
is($rels[0]->release_tag, 'ITAG1');
is($rels[1]->release_tag, 'ITAG1_pre');
is($rels[2]->release_tag, 'ITAG_devel');

#test release_number
is($rels[0]->release_number,1,'correct release number');
is($rels[1]->release_number,1,'correct release number');
is($rels[2]->release_number,0,'correct release number');

#test dir_modtime
cmp_ok($rels[0]->dir_modtime,'<=',time(),'reasonable dir_modtime');
#test text_description
like($rels[2]->text_description,qr/development/,'text_description for development release');
unlike($rels[1]->text_description,qr/development/,'text_description for non-development release');

# now let's make a new release and test it out
my $newrel = CXGN::ITAG::Release->new( releasenum => 2);
is($newrel->dir(),"$tempdir/ITAG2_release",'got correct dirname for new release');
ok(! -d $newrel->dir(), 'and dir does not exist yet');
ok($newrel->mkdir(), 'mkdir returned true');
ok( -w $newrel->dir(), 'and now the dir exists and is writable');
is( $newrel->release_tag, 'ITAG2', 'correct release tag');
is( $newrel->get_file_info('cdna_fasta')->{file},"$tempdir/ITAG2_release/ITAG2_cdna.fasta", 'got correct cDNA fasta file name');

# test tarfile handling
is( $newrel->tarfile_name, "$tempdir/ITAG2_release.tar.gz", 'got correct tarfile name');
ok( ! -f $newrel->tarfile_name, 'tarfile does not exist yet');
ok( $newrel->make_tarfile, 'make_tarfile returns true');
ok( -f $newrel->tarfile_name, 'and now tarfile exists');
