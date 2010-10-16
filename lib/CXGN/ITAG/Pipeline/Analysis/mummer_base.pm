package CXGN::ITAG::Pipeline::Analysis::mummer_base;
use strict;
use warnings;
use English;
use Carp;

use Path::Class;
use File::Basename;

use Bio::SeqIO;

use CXGN::Tools::Run;

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

sub launch_job {
    my ( $self, $batch, $seqname ) = @_;

    my ( $outfile, $gff3_out_file ) =
        map { $self->cluster_temp( $seqname, $_ ) }
       'txt', 'gff3';

    my $seq_file = $batch->pipeline
                         ->analysis('seq')
                         ->files_for_seq( $batch, $seqname );

    my $job = $self->cluster_run_class_method(
        $batch,
        $seqname,
        'run_mummer',
        $seqname,
        $self->genome_file,
        $seq_file,
        $outfile,
        $gff3_out_file,
        #{ err_file => $errfile,
        #},
       );

    my ( $out_dest, $gff_out_dest ) = $self->files_for_seq($batch,$seqname);

    return ( $job,
             [ $gff3_out_file => $gff_out_dest ],
             [ $outfile       => $out_dest     ],
           );
}

sub run_mummer {
    my ( $class, $seqname, $genome_fasta_file, $seq_file, $outfile, $gff3_outfile ) = @_;

    my ($mummer_param_singles, $mummer_param_values) = $class->mummer_params;

    my $mummer_job = CXGN::Tools::Run->run(
        'mummer',

        @$mummer_param_singles,
        %$mummer_param_values,

        -F => $genome_fasta_file,

        $seq_file,

        { die_on_destroy => 1,
          working_dir => File::Spec->tmpdir,
          on_completion => sub {
              my $job = shift;
              my $mummer_results = $class->_parse_mummer( $job->out_file );

              # print our outfile
              my $job_out = file( $job->out_file )->openr;
              my $out_fh  = file( $outfile )->openw;
              $out_fh->print( $_ ) while <$job_out>;

              #convert the mummer results to gff3
              my $gff_fh = file($gff3_outfile)->openw;
              $class->_mummer_to_gff3( $mummer_results, $gff_fh );
          },
        },
       )
      or die "failed to run mummer on fasta files '$genome_fasta_file' versus '$seq_file'";

}

sub _mummer_to_gff3 {
    my ( $class, $mummer_results, $gff_fh ) = @_;

    $gff_fh->print(<<'');
##gff3-version 3

    return unless @$mummer_results;

    my $seq_length = $mummer_results->[0]->{q_seq_len}
        or die "must run mummer with -L switch to include sequence lengths in output";

    $gff_fh->print(<<"");
##sequence-region $mummer_results->[0]->{query} 1 $seq_length

    for my $match ( @$mummer_results ) {
        my ( $start, $end ) = ( $match->{q_start}, $match->{q_start} + $match->{match_len} - 1 );
        ### start: $start
        ### end: $end
        ### match: $match

        if ( $match->{query_reversed} ) {
            ( $start, $end ) = ( $end, $start );
            $_ = $match->{q_seq_len} - $_ + 1 for $start, $end;
        }

        $gff_fh->print( join( "\t",
                              $match->{query},
                              $class->gff3_source,
                              'match',
                              $start,
                              $end,
                              '.',
                              ( $match->{query_reversed} ? '-' : '+' ),
                              '.',
                              join( ';',
                                    'Name='.$match->{subject},
                                    'Target='
                                      .$match->{subject}
                                      .' '.$match->{s_start}
                                      .' '.($match->{s_start} + $match->{match_len} - 1 )
                                      .' +'
                                  ),
                             ),
                        "\n"
                       );
    }

    return;
}


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
