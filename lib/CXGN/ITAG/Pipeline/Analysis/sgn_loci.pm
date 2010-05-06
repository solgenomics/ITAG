package CXGN::ITAG::Pipeline::Analysis::sgn_loci;
use strict;
use warnings;
use English;
use Carp;

use File::Temp qw/tempdir/;
use File::Copy qw/move/;

use Fatal qw/ open close chdir /;

use Bio::SeqIO;

use CXGN::BlastDB;
use CXGN::Tools::Class qw/parricide/;
use CXGN::Tools::Identifiers qw/parse_identifier/;
use CXGN::Tools::Wget qw/wget_filter/;

=head1 NAME

CXGN::ITAG::Pipeline::Analysis::sgn_loci - pipeline analysis to
annotate probable locations of genetic loci from SGN

=head1 BASE CLASS(ES)

L<CXGN::ITAG::Pipeline::Analysis>

=cut

use base qw/CXGN::ITAG::Pipeline::Analysis/;

=head1 METHODS

implements all abstract methods defined in L<CXGN::ITAG::Pipeline::Analysis>,
also overrides FILL IN OVERRIDES HERE

=cut

#this analysis is locally runnable
sub locally_runnable { 1 }

#TODO:
#  - for each accepted match
#     - look up the gene annotation in the eug annots
#     - make another feature with the same extents establishing the relationship
#        - include data about the match in that feature

#and this is the routine to run it
sub run {
    my ($self,$batch) = @_;

    # WE ARE LOOKING FOR:
    #  - almost-complete matches of locus sequences onto bacs

    my ($loci_blastdb) = CXGN::BlastDB->search( file_base => 'loci/sgn_loci' );
    $loci_blastdb or die "no SGN Loci BLAST db defined in database (looking for file_base loci/sgn_loci)";
    $loci_blastdb->source_url or die "no source_url defined in database for loci/sgn_loci blast db!";

    my @ops = map {
        my $seqname = $_;

        my ($genomic_seqs) = $batch->pipeline->analysis('seq')->files_for_seq($batch,$seqname);

        #get the fasta file for this sequence name
        my ($eug_annots,undef,undef,$eug_cdna_seqs) = $batch->pipeline->analysis('renaming')->files_for_seq($batch,$seqname);
        $eug_cdna_seqs && -f $eug_cdna_seqs
            or die "expected sequence file '$eug_cdna_seqs' not found";


        my ( $blat_temp, $gff3_temp ) =
            map $self->cluster_temp($seqname,$_), "blat", "gff3";


        my $blat_run = $self->cluster_run_class_method( $batch,
							$seqname,
							'run_blat',
							$loci_blastdb->source_url,
							$eug_annots,
							$genomic_seqs,
							$eug_cdna_seqs,
							$blat_temp,
							$gff3_temp,
                                                       );

        #the files we will be publishing
        my ( $blat_out, $gff3_out ) = $self->files_for_seq($batch,$seqname);

        #now make the file moves we'll do if all the analyses succeed
        ( [ $blat_temp => $blat_out, $blat_run ],
          [ $gff3_temp => $gff3_out, $blat_run ],
        )
    } $batch->seqlist;

    # wait for all of the cluster jobs to finish
    sleep 2 while grep $_->[2]->alive, @ops;

    # move all the analyzed files into position, attempting to roll it
    # back if it fails
    $self->atomic_move( @ops );
}


sub run_blat {
    my ( $class, $loci_source_url, $eug_annots, $genomic_seqs, $eug_cdna_seqs, $blat_temp, $gff3_temp ) = @_;

    ######

    my $loci_seqs = wget_filter( $loci_source_url,
                                 => $class->local_temp('loci.seq')
                               );
    -s $loci_seqs or die "could not fetch SGN loci sequences";

    #find all the gene elements in the eugene annotation and index them by name
    my %eug_genes = do {
        open my $ea, '<', $eug_annots;
        my @g;
        while (<$ea>) {
            next unless /ID=gene:/;
            my @f = split;
            my ($name) = $f[8] =~ /ID=gene:([^\s;]+)/;
            $name or die "could not extract gene name from attrs '$f[8]'";
            push @g, $name => \@f
        }
        @g
    };

    #make an index of the sequence lengths
    my %seqlengths;
    { my $i = Bio::SeqIO->new( -format => 'fasta', -file => $genomic_seqs );
      while ( my $s = $i->next_seq ) {
          $seqlengths{ $s->display_id } = $s->length;
      }
    }

    #BLAT the sequence against the SGN loci seqs

    my @cmd = ( 'blat',
                '-minIdentity=85',
                '-maxIntron=70000',
                $eug_cdna_seqs,
                $loci_seqs,
                $blat_temp
              );
    system @cmd;
    die "$! running @cmd" if $CHILD_ERROR;

    #convert the blat output into gff3
    $class->blat2gff3( \%eug_genes,
                       \%seqlengths,
                       $blat_temp =>
                       $gff3_temp
                     );

}

sub blat2gff3 {
    my ($class, $eug_genes, $seqlengths, $blat_temp, $gff3_temp ) = @_;
    #convert the blat output into gff3
    open my $psl_in, '<', $blat_temp;
    open my $gff3, '>', $gff3_temp;
    print $gff3 "##gff-version 3\n";
    my %uniq_ctr;
    my $curr_region = '';
    while(my $line = <$psl_in> ) {
      if($line =~ /^\d/) {
	my ($match,$mismatch,$repmatch,$ns,$qgapc,$qgapb,$tgapc,$tgapb,$strand,$qname,$qsize,$qstart,$qend,$tname,$tsize,$tstart,$tend,$blockc,$blocksizes)
	  = split /\s+/,$line;

	#convert the PSL 0-based half-open coords to 1-based closed coords
	$qstart++;
	$tstart++;

	#calculate percentage match, based on the locus sequence
	my $match_pct = $match/$qsize;

	#check the criteria for this match to be included
	next unless
	  $match_pct >= .85  #< match must cover 90% of locus sequence
	    ;

	#start this line as a copy of the gene: eugene line, just change some things
	my $gene_name = $tname; #< remove the splice variant number, if present, to get the gene name
	$gene_name =~ s/(\.\d+)\.\d+$/$1/;
	my $gene_line = $eug_genes->{$gene_name} or die "gene '$gene_name' not found in eugene annotations";


	unless( $curr_region eq $gene_line->[0] ) {
	  my $seqlength = $seqlengths->{$gene_line->[0]}
              or die "genomic seq $gene_line->[0] appears in BLAT output, but not in eugene output?!?!";

	  print $gff3 "##sequence-region $gene_line->[0] 1 $seqlength\n";
	  $curr_region = $gene_line->[0];
	}

	my $p = parse_identifier($qname,'sgn_locus_sequence')
	  or next;

	$gene_line->[1] = 'SGN_Loci';
	$gene_line->[2] = 'match';
	$gene_line->[5] = $match_pct;
	$gene_line->[8] = join(';',map join('=',@$_),
			       #[ ID => "$tname-$qname-".++$uniq_ctr{"$tname-$qname"} ],
			       [ Name => $qname ],
			       [ Target => join ' ',$qname,$qstart,$qend,'+'],
			       [ sgn_locus_id => $p->{id} ],
			       [ cdna_seq_name => $tname ],
			       [ blat_cdna_match_pct => sprintf('%0.2f',$match_pct*100)],
			       [ blat_cdna_base_mismatches => $mismatch ],
			       [ blat_cdna_q_gaps => $qgapc ],
			       [ blat_cdna_q_gap_bases => $qgapb ],
			       [ blat_cdna_t_gaps => $tgapc ],
			       [ blat_cdna_t_gap_bases => $tgapb ],
			      );
	
	print $gff3 join("\t",@$gene_line)."\n";
      }
    }
}


=head1 AUTHOR(S)

Robert Buels

=cut

sub DESTROY {
  return parricide(shift,our @ISA);
}

###
1;#do not remove
###
