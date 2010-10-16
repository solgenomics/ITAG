package CXGN::ITAG::Pipeline::Analysis::local_cluster_base;
use strict;
use warnings;

use base qw/CXGN::ITAG::Pipeline::Analysis/;

sub locally_runnable { 1 }

sub run {
  my ( $self, $batch ) = @_;

  # for each seq file, set up and run an analysis job on the
  # cluster, format its output, and queue it up to be copied into
  # place
  my @jobs;
  my @ops;
  for my $seqname ( $batch->seqlist ) {

      my ( $job, @job_ops ) = $self->launch_job( $batch, $seqname );

      push @jobs, $job;
      push @ops,  @job_ops;

      # check for any failed jobs, so that we don't have to submit
      # all the jobs before we die.  submitting all the jobs could
      # take a very long time if there are many jobs, because the
      # job submission will block and wait for the torque server to
      # clear a bit if there are too many jobs in its queue.
      $_->alive for @jobs;
  }

  # wait for all the jobs to finish (also running their
  # on_completion hooks)
  $_->wait for @jobs;

  #atomically move the results into position
  $self->atomic_move(@ops);
}

1;
