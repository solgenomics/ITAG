#!/usr/bin/perl
use strict;
use warnings;
use English;

use File::Temp qw/tempdir/;
use File::Spec;

use CXGN::DB::Connection;

use CXGN::Tools::Run;

use Test::More tests => 42;
use Test::Warn;

BEGIN {
  use_ok(  'CXGN::ITAG::Pipeline'  )
    or BAIL_OUT('could not include the module being tested');
}

$ENV{CXGNITAGPIPELINEANALYSISTESTING} = 1; #< bypass some file
                                           #permissions checks that
                                           #can't always be satisfied
                                           #in a test script

my $tmp = tempdir(File::Spec->catdir(File::Spec->tmpdir,'pipeline-t-XXXXXXXX'),CLEANUP=>1);
#diag "using tempdir $tmp\n";
my $pipe = CXGN::ITAG::Pipeline->create(basedir => $tmp);

isa_ok($pipe,'CXGN::ITAG::Pipeline');
is($pipe->version,0);
ok(-d File::Spec->catdir($tmp,'pipe000'),"made pipe dir in $tmp");

is_deeply([$pipe->valid_file_exts],
	  [],
	  'with no global defs file, no valid file extensions',
	 );
open my $gdef, '>', "$tmp/pipe000/global.def.txt" or die "$! writing temp def file";
print $gdef <<EOF;
file_exts foosta foogp
feature_ontology_url fake_url
EOF
close $gdef;

is_deeply([$pipe->valid_file_exts],
	  ['foosta','foogp'],
	  'with a global defs file, has valid file extensions',
	 );


is_deeply([$pipe->list_analyses],
	  [],
	  'with no def files, analysis list should be empty');


warning_like {
  is($pipe->analysis('nonexistentanalysis'), undef, 'requesting a nonexistent analysis returns undef');
} qr/no analysis found with tag 'nonexistentanalysis'/, 'correct warning for previous test';
#no def file, trying to get the seq analysis should fail
warning_like {
  is($pipe->analysis('seq'),undef,'get analysis fails with no def file');
} qr/no def file/, 'correct warning for previous test';

#make a bad def file for the 'seq' analysis
{ open my $def,'>',$pipe->analysis_def_file('seq')
    or die "cannot open analysis def file for writing, wtf?";
#  warn "opened def file ".$pipe->analysis_def_file('seq');
  print $def 'owner_user ' => (getpwuid($UID))[0] , "\n";
  print $def 'owner_email ' => 'sgn-bugs@sgn.cornell.edu',"\n";
  print $def "fake_key\n";
  close $def;
  #diag "successfully wrote to def file ".$pipe->analysis_def_file('seq');
}

#diag "reopening pipeline";
#reopen the pipeline, rereading analysis stuff
$pipe = CXGN::ITAG::Pipeline->open(basedir => $tmp, version => 0);

isa_ok($pipe->analysis('seq'),'CXGN::ITAG::Pipeline::Analysis','get analysis succeeds with bad def file');
ok(scalar($pipe->analysis('seq')->errors),'analysis has errors with bad def file');
#warn "1 seq analysis is ".$pipe->analysis('seq');
is_deeply([$pipe->list_analyses],
	  ['seq'],
	  "analyses with bad defs still show up in the list",
	 );

#make a good def file for the 'seq' analysis
{ sleep 1; #makes sure the file modtimes are different
  open my $def,'>',$pipe->analysis_def_file('seq')
    or die "cannot open analysis def file for writing, wtf?";
  print $def 'owner_user '.(getpwuid($UID))[0] , "\n";
  print $def 'owner_email  sgn-bugs@sgn.cornell.edu',"\n";
  print $def "depends_on\n";
  print $def "produces_files fonebone:foosta monkey:foogp\n";
#  warn "2 wrote to def file ".$pipe->analysis_def_file('seq');
}

#reopen the pipeline, rereading analysis stuff
$pipe = CXGN::ITAG::Pipeline->open(basedir => $tmp, version => 0);

my $seq_an = $pipe->analysis('seq');
#warn "2 seq analysis is $seq_an";
isa_ok($seq_an,'CXGN::ITAG::Pipeline::Analysis','get analysis succeeds with good def file');
is(scalar($seq_an->errors),0,'with good def file, analysis has no errors');
#diag map "- $_\n", $seq_an->errors;
is_deeply([$pipe->list_analyses],
	  ['seq'],
	  'new seq analysis shows up in analysis list');

my @batches = $pipe->list_batches;
is_deeply(\@batches,[],'initially pipeline has no batches');

my $batch;
warnings_like {
  $batch = $pipe->create_batch(size => 10, seq_source => 'CXGN::ITAG::SeqSource::Fake');
} qr/new batch size.+only 3 were available/, 'new batch made low seq number warning';
is($batch->batchnum,1,'initial batch created is batch number 1');
#diag "made batch with seqnames ".join(', ',$batch->seqlist),"\n";

is($batch->size,3,'batch size is set correctly');

#make sure all the analysis dirs are created
my @dirs = map {s/\/$//; s!.+/!!; $_} glob($batch->dir.'/*/');
is_deeply([sort @dirs],
	  [sort $pipe->list_analyses],
	  'dirs are created for all enabled analyses',
	 );

#seq-fetching analysis should be ready at this point
is($pipe->analysis('seq')->status(1),'ready','seq analysis is ready');

#now write a 'running 1' to the control file and check that the seq analysis reports running
system "echo running 1 > $tmp/pipe000/batch001/seq/control.txt";
flush_cache();
is($pipe->analysis('seq')->status(1),'running','seq analysis is running');
system "echo running 0 > $tmp/pipe000/batch001/seq/control.txt";
flush_cache();
is($pipe->analysis('seq')->status(1),'ready','seq analysis is ready when control file is on 0');
system "echo running 1 > $tmp/pipe000/batch001/seq/control.txt";
flush_cache();
is($pipe->analysis('seq')->status(1),'running','seq analysis is running');
system "rm $tmp/pipe000/batch001/seq/control.txt";
flush_cache();
is($pipe->analysis('seq')->status(1),'ready','seq analysis is running');

#run the seq-fetching analysis on the batch
$pipe->analysis('seq')->run(1,seq_source => 'CXGN::ITAG::SeqSource::Fake');

#make sure analysis has 'done' status with no md5 sum file
flush_cache();
is( $pipe->analysis('seq')->status(1),'done','seq analysis shows done after running with no md5sum' )
  or do{ my $dir = $pipe->analysis('seq')->work_dir(1);  diag "ls of work dir '$dir':\n".`ls $dir`; };

#make an md5sum of its files
system "cd $tmp/pipe000/batch001/seq; md5sum *.* > md5sums.txt";

#make sure analysis has 'validating' status now
flush_cache();
is( $pipe->analysis('seq')->status(1),'validating','seq analysis shows validating after running' )
  or do{ my $dir = $pipe->analysis('seq')->work_dir(1);  diag "ls of work dir '$dir':\n".`ls $dir`; };

#now test analysis validation
$pipe->analysis('seq')->run_intense_validation(1);

flush_cache();
is( $pipe->analysis('seq')->status(1),'done','seq analysis output validated OK' )
  or do{ my $dir = $pipe->analysis('seq')->work_dir(1);  diag "ls of work dir '$dir':\n".`ls $dir`; };

#reopen the pipeline, rereading analysis stuff
$pipe = CXGN::ITAG::Pipeline->open(basedir => $tmp, version => 0);

my $workdir = $pipe->analysis('seq')->work_dir(1);
my @seqfiles = glob("$workdir/*");
is(scalar(@seqfiles),7,'7 files present in sequence fetch dir')
  or diag('actual files present: '.CXGN::Tools::Run->run("ls ".$pipe->analysis('seq')->work_dir(1).'/*')->out);

#remove one of the output files.  seq should report error after that
my $cmd = "rm ".shift(@seqfiles);
system $cmd;
#diag $cmd;

#reopen the pipeline, rereading analysis stuff
$pipe = CXGN::ITAG::Pipeline->open(basedir => $tmp, version => 0);

is($pipe->analysis('seq')->status(1),'error','seq analysis shows error if you delete one of its output files');
$cmd = "rm ".join(' ',@seqfiles);

system $cmd;
#diag $cmd;

#reopen the pipeline, rereading analysis stuff
$pipe = CXGN::ITAG::Pipeline->open(basedir => $tmp, version => 0);
is($pipe->analysis('seq')->status(1),'ready','seq analysis shows ready if you delete all of its output files')
  or diag `ls $workdir`, map "- $_\n", $pipe->analysis('seq')->errors;

#delete the batch
$batch->delete;

#make sure list_batches doesn't return this  batch
@batches = $pipe->list_batches;
is_deeply(\@batches,[],'after delete, pipeline has no batches');

#make another batch
warnings_like {
  $batch = $pipe->create_batch(size => 10, seq_source => 'CXGN::ITAG::SeqSource::Fake');
} qr/new batch size.+only 3 were available/, 'new batch made low seq number warning';
is($batch->batchnum,2,'next batch created is batch number 2, even after deletion');
#diag "made batch with seqnames ".join(', ',$batch->seqlist),"\n";
$pipe->analysis('seq')->run($batch,seq_source => 'CXGN::ITAG::SeqSource::Fake');
is($pipe->analysis('seq')->status($batch),'done','second seq analysis is done');


#create another pipeline
my $pipe2 = CXGN::ITAG::Pipeline->create(basedir => $tmp);
isa_ok($pipe2,'CXGN::ITAG::Pipeline');
is_deeply([$pipe2->valid_file_exts],
	  [qw/foosta foogp/],
	  'pipeline valid_file_exts propagated to new pipeline version',
	 );
is($pipe2->version,1);
ok(-d File::Spec->catdir($tmp,'pipe001'),"made second pipe dir in $tmp/pipe001");

is_deeply([$pipe2->list_analyses],
	  ['seq'],
	  'analysis def files also got carried over',
	 );


#diag `find $tmp`;
sub flush_cache {
  Memoize::flush_cache('CXGN::ITAG::Pipeline::Analysis::status');
  sleep 1;
}
