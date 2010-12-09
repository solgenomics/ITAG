package CXGN::ITAG::Pipeline::Analysis::blat_base;
use strict;
use warnings;
use Carp;

use autodie ':all';

=head1 NAME

CXGN::ITAG::Pipeline::Analysis::blat_base - pipeline analysis base class for BLAT analyses

=head1 BASE CLASS(ES)

L<CXGN::ITAG::Pipeline::Analysis>

=cut

use base qw/ CXGN::ITAG::Pipeline::Analysis::local_cluster_base /;

sub launch_job {
    my ( $self, $batch, $seqname ) = @_;

    my ( $outfile, $gff3_out_file ) =
        map { $self->cluster_temp( $seqname, $_ ) }
       'txt', 'gff3';

    my ($seq_file) = $batch->pipeline
                           ->analysis('repeats')
                           ->files_for_seq( $batch, $seqname );

    my $job = $self->cluster_run_class_method(
        $batch,
        $seqname,
        'run_blat',
        $seqname,
        $self->query_file,
        $seq_file,
        $outfile,
        $gff3_out_file,
       );

    my ( $out_dest, $gff_out_dest ) = $self->files_for_seq( $batch, $seqname );

    return ( $job,
             [ $gff3_out_file => $gff_out_dest ],
             [ $outfile       => $out_dest     ],
           );
}

sub blat_params {
    return (
        '-minIdentity=85',
        '-fastMap',
       );
}

sub gff3_source {
    'ITAG_blat_base';
}

sub run_blat {
    my ( $class, $seqname, $query_file, $subject_file, $raw_outfile, $gff3_outfile ) = @_;

    my @cmd = ( 'blat',
                $class->blat_params,
                $subject_file,
                $query_file,
                $raw_outfile,
               );
    system @cmd;
    die "$! running @cmd" if $?;

    #convert the blat output into gff3
    $class->blat2gff3(
        raw_out      => $raw_outfile,
        seqname      => $seqname,
        gff3_out     => $gff3_outfile,
        query_seqs   => $query_file,
        subject_seqs => $subject_file,
       );

}

sub munge_gff3 {
    #my ( $class, $args, $gff3, $attr ) = @_;
}


sub blat2gff3 {
    my ( $class, %args ) = @_;

    #convert the blat output into gff3
    open my $psl_in, '<', $args{raw_out};
    open my $gff3, '>', $args{gff3_out};

    $gff3->print( "##gff-version 3\n" );

    my $curr_region = '';
    while(my $line = <$psl_in> ) {
        if($line =~ /^\d/) {
            my ($match,$mismatch,$repmatch,$ns,$qgapc,$qgapb,$tgapc,$tgapb,$strand,$qname,$qsize,$qstart,$qend,$tname,$tsize,$tstart,$tend,$blockc,$blocksizes)
                = split /\s+/,$line;

            #convert the PSL 0-based half-open coords to 1-based closed coords
            $qstart++;
            $tstart++;

            #calculate percentage match, based on the query sequence
            my $percent_id = sprintf "%.2f",
		($match + $repmatch)/( $match + $mismatch + $repmatch);

            my @gff3 = (
                $tname,
                $class->gff3_source,
                'match',
                $tstart,
                $tend,
                $percent_id,
                $strand,
                '.',
               );

            my %attr = (
                Name   => $qname,
                Target => "$qname $qstart $qend +",
               );

            $class->munge_gff3( \%args, \@gff3, \%attr );

            push @gff3, join( ';', map { "$_=$attr{$_}" } sort keys %attr );

            $gff3->print( join( "\t", @gff3 ), "\n" );
        }
    }
}

=head1 AUTHOR(S)

Robert Buels

=cut

1;

