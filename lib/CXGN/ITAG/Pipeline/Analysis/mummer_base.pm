package CXGN::ITAG::Pipeline::Analysis::mummer_base;
use strict;
use warnings;

use autodie ':all';


use Bio::SeqIO;

use List::Util qw/shuffle/;
use Scalar::Util qw/ blessed /;

use Path::Class;
use URI::Escape;

use CXGN::Tools::Run;
use CXGN::Tools::List 'balanced_split_sizes';

=head1 NAME

CXGN::ITAG::Pipeline::Analysis::mummer_base - pipeline analysis base class for producing alignments using mummer

=head1 BASE CLASS(ES)

L<CXGN::ITAG::Pipeline::Analysis>

=cut

use base qw/CXGN::ITAG::Pipeline::Analysis::local_cluster_base/;

=head1 METHODS

implements all abstract methods defined in L<CXGN::ITAG::Pipeline::Analysis>,

=cut


##### methods to override in subclasses
sub mummer_params {

    return (
        [ qw( -mum -b -n -L ) ],

        {
            -l => 1000,
        }
       );
}
sub gff3_source {
    'ITAG_genome'
}
sub genome_file {
    die 'must implement genome_file in subclass!';
}


######## implemented methods

sub run {
  my ( $self, $batch ) = @_;

  # try to avoid having all the big stuff on one node
  my @query_files = shuffle $self->split_query;

  # for each seq file, set up and run an analysis job on the
  # cluster, format its output, and queue it up to be copied into
  # place
  my @jobs;
  my @move_files;
  for my $seqname ( $batch->seqlist ) {

      my ( $outfile, $gff3_out_file ) =
          map { $self->cluster_temp( $seqname, $_ ) }
          'txt', 'gff3';

      push @jobs, $self->launch_jobs( $batch, $seqname, \@query_files, $outfile, $gff3_out_file );

      my ( $out_dest, $gff_out_dest ) =
          $self->files_for_seq( $batch, $seqname );

      push @move_files, [ $gff3_out_file => $gff_out_dest ],
                        [ $outfile       => $out_dest     ];

      # check for any failed jobs, so that we don't have to submit
      # all the jobs before we die.  submitting all the jobs could
      # take a very long time if there are many jobs, because the
      # job submission will block and wait for the torque server to
      # clear a bit if there are too many jobs in its queue.
      $_->alive for @jobs;
  }

  # wait for all the jobs to finish (also running their completion hooks)
  $_->wait for @jobs;

  #atomically move the results into position
  $self->atomic_move( @move_files );
}

sub split_query {
    my ( $self ) = @_;

    # count the seqs in the file
    my $count = $self->_seq_count( $self->query_file );

    # split into 20 jobs (or fewer if fewer than 20 seqs)
    my $job_sizes = balanced_split_sizes( 50, $count );

    my $query_seqs = Bio::SeqIO->new( -file => $self->query_file, -format => 'fasta' );

    # make the seq files and cluster jobs for each batch of query seqs
    my @files;
    for( my $job_number = 0; $job_number < @$job_sizes; $job_number++ ) {
        my $job_size = $job_sizes->[$job_number];

        my $job_seq_file = $self->cluster_temp( 'query_seqs', "$job_number.fasta" );
        my $job_seqs = Bio::SeqIO->new( -file => ">$job_seq_file", -format => 'fasta' );

        $job_seqs->write_seq( $query_seqs->next_seq )
            for 1..$job_size;

	push @files, $job_seq_file;
    }

    return @files;
}

sub launch_jobs {
    my ( $self, $batch, $seqname, $query_files, $outfile, $gff3_out_file ) = @_;

    my ( $reference_seq_file ) =
        $batch->pipeline
              ->analysis('seq')
              ->files_for_seq( $batch, $seqname );

    return map {
	my $job_seq_file = $_;
	$self->cluster_run_class_method(
            $batch,
            $seqname,
            'run_mummer',
            $seqname,
            $job_seq_file,
            $reference_seq_file,
            {
                on_completion => sub {
                    my $job = shift;

                    #concatenate the output files
                    my $job_out = file( $job->out_file )->openr;
                    my $out_fh  = file( $outfile )->open('>>');
                    $out_fh->print( $_ ) while <$job_out>;

                    #convert the mummer results to gff3
                    my $mummer_results = $self->_parse_mummer( $job->out_file );
                    my $gff_fh = file( $gff3_out_file )->open('>>');
                    $self->_mummer_to_gff3( $mummer_results, $gff_fh );
                },
            },
           )
    } @$query_files;
}
sub _seq_count {
    my ( $class, $file ) = @_;
    open my $f, '<', $file;
    my $c = 0;
    while( my $l = <$f> ) {
        $c++ if $l =~ /^>/;
    }
    return $c;
}

sub run_mummer {
    my ( $class, $seqname, $genome_fasta_file, $seq_file ) = @_;

    my ($mummer_param_singles, $mummer_param_values) = $class->mummer_params;

    my $mummer_job = CXGN::Tools::Run->run(
        'mummer',

        @$mummer_param_singles,
        %$mummer_param_values,

        -F => $genome_fasta_file,

        $seq_file,

        { die_on_destroy => 1,
          working_dir => File::Spec->tmpdir,
	  vmem => '1500mb',
        },
       )
      or die "failed to run mummer on fasta files '$genome_fasta_file' versus '$seq_file'";

    print $mummer_job->out;
}

sub _mummer_to_gff3 {
    my ( $class, $mummer_results, $gff_fh ) = @_;

    return unless @$mummer_results;

    # print gff version and sequence-region pragmas if this is the
    # first gff3 in the file
    #unless( blessed($gff_fh) && $gff_fh->can('tell') && $gff_fh->tell ) {
    unless( $gff_fh->tell ) {
        $gff_fh->print(<<'');
##gff-version 3

        my $seq_length = $mummer_results->[0]->{q_seq_len}
            or die "must run mummer with -L switch to include sequence lengths in output";
        $gff_fh->print(<<"");
##sequence-region $mummer_results->[0]->{query} 1 $seq_length

    }

    for my $match ( @$mummer_results ) {
        my ( $start, $end ) = ( $match->{q_start}, $match->{q_start} + $match->{match_len} - 1 );
        ### start: $start
        ### end: $end
        ### match: $match

        if ( $match->{query_reversed} ) {
            ( $start, $end ) = ( $end, $start );
            $_ = $match->{q_seq_len} - $_ + 1 for $start, $end;
        }

        my @fields = (
            $match->{query},
            $class->gff3_source,
            'match',
            $start,
            $end,
            '.',
            ( $match->{query_reversed} ? '-' : '+' ),
            '.',
           );

        my %attrs = (
            Name   => $match->{subject},
            Target => join( ' ',
                            $match->{subject},
                            $match->{s_start},
                            ($match->{s_start} + $match->{match_len} - 1 ),
                            '+',
                           ),

           );

        $class->munge_gff3( {}, \@fields, \%attrs );

        $gff_fh->print( join( "\t",
                              @fields,
                              join( ';',
                                    map { "$_=".uri_escape($attrs{$_}, ';=%&,') }
                                    sort keys %attrs
                                   )
                             ),
                        "\n",
                       );
    }

    return;
}

sub munge_gff3 {} #< default does no gff3 munging

# returns arrayref as
#  [
#    { match info }
#    ...
#  ]
sub _parse_mummer {
    my ( $class, $mummer_output_file ) = @_;

    open (my $f, $mummer_output_file)
        or die "Can't open file $mummer_output_file for reading...";

    my @return;

    my $query = "";
    my $reverse;
    my $qlength;
    local $_;
    while ( my $line = <$f>) {
        chomp $line;

        if ( $line =~ / ^> \s* (\S+) \s+ (Reverse)? \s* Len \s* = \s* (\d+) /xi ) {
            ( $query, $reverse, $qlength ) = ( $1, $2, $3 );
            $reverse = $reverse ? 1 : 0;
        } else {
            $line =~ s/^\s+|\s+$//g;
            my ($subject, $query_start, $subject_start, $length) = split /\s+/, $line;

            push @return, {
                query          => $query,
                query_reversed => $reverse,
                subject        => $subject,
                q_start        => $query_start,
                s_start        => $subject_start,
                match_len      => $length,
                q_seq_len      => $qlength,
            };
        }
    }

    return \@return;
}


=head1 AUTHOR(S)

Robert Buels

=cut

###
1;#do not remove
###
