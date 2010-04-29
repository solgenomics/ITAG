use strict;
use warnings;

use File::Spec::Functions;
use FindBin;
use IPC::Cmd qw/ can_run /;
use Test::More;
use Path::Class;
use Storable;

use Capture::Tiny qw/ capture /;

use CXGN::ITAG::Pipeline;

$ENV{CXGNITAGPIPELINEANALYSISTESTING} = 1;

my $class = 'CXGN::ITAG::Pipeline::Analysis::OutputValidator::gff3_file_format';

use_ok $class;

# test the job file creation
my $pipe = CXGN::ITAG::Pipeline->open( basedir => dir()->subdir('t')->subdir('data')->subdir('test1') );

is( $pipe->version, 1, 'opened test pipeline' );

my $batch = $pipe->batch(4);
my $analysis = $pipe->analysis('sgn_unigenes');

# test files_to_validate
my @gff3_files = $class->files_to_validate( $analysis, $batch );
ok( scalar(@gff3_files), 'got some gff3 files' );
ok( -f, 'file exists' ) for @gff3_files;

# test _make_jobfiles
my @jobfiles = $class->_make_jobfiles( $analysis, $batch, \@gff3_files );
is( scalar(@jobfiles), scalar(@gff3_files), 'right number of jobs' );
foreach my $f (@jobfiles) {
    my $args = retrieve($f);
    is( $args->{batch}->batchnum, $batch->batchnum, 'got right batch' );
    is( $args->{self}, $class, 'got right class' );
}

# test _read_job_failures
is_deeply( $class->_read_job_failures([ mockjob->new ]),
           { testfile => 'failure text' },
           '_read_job_failures works',
          );

SKIP: {
    skip 'validate_gff3.pl not in path, cannot test validation', 2
        unless can_run('validate_gff3.pl');

    my ( $job_stdout, $job_stderr ) = capture {
        eval { $class->_run( $jobfiles[0] ) };
        warn $@ if $@;
    };

    like( $job_stdout, qr/^{/, 'job stdout looks OK' )
        or diag $job_stdout;
    is( $job_stderr, '', 'nothing on job stderr' )
        or diag $job_stderr;

    my $queue = $ENV{ITAG_TEST_CLUSTER_QUEUE}
        or skip 'set ITAG_TEST_CLUSTER_QUEUE to test cluster job submission', 2;

    local $ENV{PATH} = catdir($FindBin::RealBin, updir(), updir(), 'bin' ).":$ENV{PATH}";
    my ($stdout, $stderr) = capture {
        eval {
            $class->run_offline( $analysis, $batch, { job_queue => $ENV{ITAG_TEST_CLUSTER_QUEUE} } );
        };
        warn $@ if $@;
    };
    is($stderr, "running gff3 validate on sgn_unigenes\n", 'got expected stderr')
        or diag "stderr was:\n$stdout";
    is($stdout, '', 'got no stdout' )
        or diag "stdout was:\n$stdout";

    my @errors = $class->errors( $analysis, $batch );
    is( scalar(@errors), 0, 'no errors reported' );
}

done_testing;

package mockjob;
use Data::Dumper;

sub new {  bless {}, shift }

sub out {
    local $Data::Dumper::Terse = 1;
    Dumper { testfile => 'failure text' };
}
