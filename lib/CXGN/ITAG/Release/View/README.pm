package CXGN::ITAG::Release::View::README;
use Moose;
use namespace::autoclean;

use File::Basename;
use Template;

use CXGN::ITAG::Release::Statistics;

has 'release' => (
    is => 'ro',
    required => 1,
    );

has 'statistics' => (
    is         => 'ro',
    lazy_build => 1,
    );
sub _build_statistics {
    CXGN::ITAG::Release::Statistics->new( release => shift->release );
}

#write the readme file
sub render {
  my ( $self ) = @_;

  my $readme_template = <<'EOT';
[%- USE Comma -%]
[%- BLOCK based_count -%]
[% num | comma %] ([% num / base * 100 | format('%2.1f') | format('%5s') %]%)
[%- END -%]
[%- BLOCK stat_sparse -%]
   [% label %]
  ---------------------------------------------------
     Min     [% stat.min | comma %]
     Max     [% stat.max | comma %]
     Range   [% stat.max-stat.min | comma %]
     Mean    [% stat.mean | format('%0.1f') | comma %]
     StdDev  [% stat.standard_deviation | format('%0.1f') | comma %]
[%- END -%]
[%- BLOCK stat_full -%]
[% PROCESS stat_sparse %]
     Median  [% stat.median | comma %]
     Frequency Distribution:
            Bin  Frequency
     [% f=stat.frequency_distribution_ref( bins.0.defined ? bins : 8 );
        FOREACH bin IN f.keys.nsort -%]
[% bin | format('%0.0f') | comma | format('%10s') %]  [% f.$bin | comma | format('%-10s') %]
     [% END -%]
[%- END -%]
[%- BLOCK megabases -%]
[% bp / 1000000 | format('%0.0f') | comma %] Mbp
[%- END -%]
[% release_tag %] [% organism %] Genome Release

Contents:
  1. Introduction
  2. Files in this release
  3. Links and other resources
  4. Release statistics

== 1. Introduction ==

The [% project_name %] ([% project_acronym %]) is pleased to announce the [% release_tag %] release of the official [% organism %] genome annotation ([% release_tag %]), covering approximately [% s.genomic_bases / genome_size * 100 | format('%0.0f') %]% of the genome, with [% s.gene_models | comma %] gene models.  This release file set was generated on [% date_str %].

In this release, [% INCLUDE based_count num=s.protein_coding_with_supporting_cdna_or_protein base=s.gene_models %] of the gene models are supported by homology to either existing ESTs or cDNA sequences, with [% INCLUDE based_count num=s.protein_coding_with_supporting_cdna_and_protein base=s.gene_models %] supported by both.  [% IF s.gene_models_with_human_desc == s.gene_models; THEN %]All[% ELSE; INCLUDE based_count num=s.gene_models_with_human_desc base=s.gene_models; END %] of the gene models are annotated with best-guess text descriptions of their function, and [% INCLUDE based_count num=s.gene_models_with_GO_terms base=s.gene_models %] have associated Gene Ontology terms describing their function.  See section 4 for more statistics describing this release.

Please send comments or questions about these annotations to: itag@sgn.cornell.edu

== 2. Files in this release ==

[% file_descriptions %]

== 3. Links and other resources ==

Sequences and annotations can also be viewed and searched on SGN: http://solgenomics.net/genomes/Solanum_lycopersicum/genome_data.pl

The fully annotated chromosome sequences in GFF version 3 format, along with Fasta files of cDNA, CDS, genomic and protein sequences, and lists of genes are available for bulk download from the SGN site at: http://solgenomics.net/itag/release/[% release_number %]/list_files.

For those who are not familiar with the GFF3 file format, the format specification can be found here: http://www.sequenceontology.org/gff3.shtml

A graphical display of the [% organism %] sequence and annotation can be viewed using SGN's genome browser. Browse the chromosomes, search for names or short sequences and view search hits on the whole genome, in a close-up view or on a nucleotide level: http://solgenomics.net/gbrowse/

SGN's BLAST services have also been updated with this dataset, available at: http://solgenomics.net/tools/blast/

[% project_acronym %] is committed to the continual improvement of the [% organism %] genome annotation and actively encourages the community to contact us with new data, corrections and suggestions.

Announcements of new releases, updates of data, tools, and other developments from [% project_acronym %] can be found on SGN: http://solgenomics.net/

== 4. Release statistics ==

4.1 Proportion of Genome Annotated

  Estimated genome size:      [% INCLUDE megabases bp=genome_size %]
  Size of annotated assembly: [% INCLUDE megabases bp=s.genomic_bases %]
  Est. proportion of genome:  [% s.genomic_bases / genome_size * 100 | format('%0.0f') %]%

4.2 Structural Annotation

  Gene model count: [% s.gene_models | comma %]
  Exon count:       [% s.exon_length.count | comma %]
  Intron count:     [% s.intron_length.count | comma %]

[% INCLUDE stat_full stat=s.gene_model_length    label='Gene model length (bp)'   bins=bins.gene_model_length    %]
[% INCLUDE stat_full stat=s.intergenic_length    label='Intergenic distance (bp)' bins=bins.intergenic_length    %]
[% INCLUDE stat_full stat=s.exons_per_gene_model label='Exons per gene model'     bins=bins.exons_per_gene_model %]
[% INCLUDE stat_full stat=s.exon_length          label='Exon length (bp)'         bins=bins.exon_length          %]
[% INCLUDE stat_full stat=s.intron_length        label='Intron length (bp)'       bins=bins.intron_length        %]

4.3 Functional Annotation

  Gene models with GO terms:   [% INCLUDE based_count num=s.gene_models_with_GO_terms base=s.gene_models %]
  Unique GO terms associated:  [% s.unique_GO_terms | comma %]
  Genes with splice variants:  [% s.genes_with_splice_variants | comma %]
  Gene models with functional
  description text:            [% INCLUDE based_count num=s.gene_models_with_human_desc base=s.gene_models %]

[% INCLUDE stat_full stat=s.GO_terms_per_mrna label='Gene Ontology terms associated, per gene model' -%]

4.4 Gene model supporting evidence

  ESTs/cDNAs aligned to the genome: [% s.mapped_ests | comma %]

  Gene models with cDNA OR protein support:       [% INCLUDE based_count num=s.protein_coding_with_supporting_cdna_or_protein  base=s.gene_models  | format('%17s') %]

  Gene models with cDNA homology support:         [% INCLUDE based_count num=s.protein_coding_with_supporting_cdna             base=s.gene_models  | format('%17s') %]
  Gene models without cDNA homology support:      [% INCLUDE based_count num=s.protein_coding_without_supporting_cdna          base=s.gene_models  | format('%17s') %]
  Gene models with protein homology support:      [% INCLUDE based_count num=s.protein_coding_with_supporting_prot             base=s.gene_models  | format('%17s') %]
  Gene models without protein homology support:   [% INCLUDE based_count num=s.protein_coding_without_supporting_prot          base=s.gene_models  | format('%17s') %]

  Gene models with both cDNA and protein support: [% INCLUDE based_count num=s.protein_coding_with_supporting_cdna_and_protein base=s.gene_models  | format('%17s') %]
  Gene models with only cDNA homology support:    [% INCLUDE based_count num=s.protein_coding_with_supporting_only_cdna        base=s.gene_models  | format('%17s') %]
  Gene models with only protein homology support: [% INCLUDE based_count num=s.protein_coding_with_supporting_only_protein     base=s.gene_models  | format('%17s') %]
  Gene models with no homology support:           [% INCLUDE based_count num=s.protein_coding_with_supporting_none             base=s.gene_models  | format('%17s') %]

EOT

  my $readme_text;
  my $tt = Template->new;
  $tt->process(
      \$readme_template,
      {
        s                 => $self->statistics,
        bins              => {
            gene_model_length =>
                $self->_bins_with_cutoff( $self->statistics->gene_model_length, 5_000, 40_000 ),
            intron_length =>
                $self->_bins_with_cutoff( $self->statistics->intron_length, 2_000, 20_000 ),
            exon_length   =>
                $self->_bins_with_cutoff( $self->statistics->exon_length, 5_000, 30_000 ),
            intergenic_length   =>
                $self->_bins_with_cutoff( $self->statistics->intergenic_length, 20_000, 500_000 ),

        },
        release_tag       => $self->release->release_tag,
        release_num       => $self->release->release_number,
        date_str          => POSIX::strftime( "%B %e, %Y", gmtime() ),
        organism          => 'Tomato',
        genome_size       => 930_000_000,
        project_name      => 'International Tomato Annotation Group',
        project_acronym   => 'ITAG',
        file_descriptions => $self->file_descriptions,
        release_dirname   => $self->release->dir_basename,
       },
      \$readme_text,
     ) || die $tt->error;
  warn $tt->error if $tt->error;

  # do some silliness to correct the word wrapping of the text so it
  # comes out looking nice
  $readme_text = $self->wrap_long_lines( $readme_text );

  return $readme_text;

}
sub _bins_with_cutoff {
    my ( $self, $stat, $binsize, $cutoff ) = @_;

    $cutoff = $stat->max if $stat->max < $cutoff;

    my $bincount = sprintf('%0.0f',($cutoff - $stat->min) / $binsize )+1;
    my @bins = ( ( map { $cutoff - $binsize*$_ } reverse 0..($bincount-2) ), $stat->max );
    return \@bins;
}

#generate a block of text describing each file in the release
sub file_descriptions {
  my ( $self ) = @_;

  my $gen_files = $self->release->get_all_files;

  my $file_descriptions = join '', sort map {
    my $bn = basename($_->{file});
    my $desc = wrap_long_lines($_->{desc},75);
    $desc =~ s/\n/\n  /g;
    "$bn\n  $_->{desc}\n\n"
  } values %$gen_files;

  #remove ending newlines
  $file_descriptions =~ s/\n+$//;

  return $file_descriptions;
}

sub wrap_long_lines {
  my $self = shift;
  my $text = shift;
  my $length = shift || 72;
  $length--;
  $text =~
    s{
      \G              # begin where previous match left off
      ([\d\D]*?)      # consume short lines
      (?:(?<=^) | \G) # and pick up at the beginning of a line,
      # or just after the previous replaced space
      (.{1,$length})       # match as many characters on this line as fit
      \ +             # followed by spaces
      (?=(\S+))       # followed by (unconsumed) nonspace
    }{ (length($2) + length($3) >= $length) ? "$1$2\n" : "$1$2 " }mexg;

  return $text;
}


1;

