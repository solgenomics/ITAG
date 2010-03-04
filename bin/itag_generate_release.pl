#!/usr/bin/env perl

=head1 NAME

itag_generate_release.pl - script to generate a genome annotation release for ITAG

=head1 SYNOPSIS

  itag_generate_release.pl [options] -d basepath release_num target_path
  itag_generate_release.pl [options] -d basepath -D target_path

  Using the current output of the ITAG pipeline, generate a genome
  release directory target_path/<release_tag>_<pre?>release/,
  populated with all the right genome release files.

  Options:

    -d itag_pipeline_path
       REQUIRED path at which to find the itag pipeline directory
       structure

    -p version
       pipeline version number to use
       Default: most recent version

    -b <num>  or -b <num>-<num> or -b <num>,<num>-<num>
       range of batches to include in this release
       Default: all batches

    -P this is a prerelease

    -D this is a development snapshot

    -G if passed, generate and write GBrowse configurations only.
       Attempts to write them to release_dir/<conf_file>, relative to
       the given dir.  Does not overwrite other files;

    -S just collect raw stats (from existing release files) and dump.
       Mostly useful for debugging stat collection.

=head1 MAINTAINER

Robert Buels

=head1 AUTHOR(S)

Robert Buels, E<lt>rmb32@cornell.eduE<gt>

=head1 COPYRIGHT & LICENSE

Copyright 2009 The Boyce Thompson Institute for Plant Research

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;
use Pod::Usage;

use POSIX;

use IO::File;
use Hash::Util qw/ lock_hash lock_keys /;
use File::Copy;
use File::Basename;
use File::Temp qw/tempfile/;
use URI::Escape;

use Data::Dumper;

use CXGN::TomatoGenome::BACPublish qw/ genbank_acc_to_seq_name seq_name_to_genbank_acc /;
use CXGN::ITAG::Pipeline;
use CXGN::ITAG::Pipeline::Analysis;
use CXGN::ITAG::Release;

use CXGN::BioTools::AGP qw/agp_parse/;

use CXGN::Tools::Identifiers qw/ parse_identifier identifier_namespace /;
use CXGN::Tools::List qw/str_in max distinct flatten/;
use CXGN::Tools::Text qw/commify_number/;

### parse and validate command line args
our %opt;
getopts('d:p:b:PDGS',\%opt) or pod2usage(1);
$opt{m} ||= 1_000_000;
$opt{P} && $opt{D} and die "-P and -D are mutually exclusive\n";
$opt{d} or pod2usage('must provide -d option');

my $release_num = $opt{D} ? 0 : shift @ARGV;
my $target_path = shift @ARGV or pod2usage();
($release_num + 0 == $release_num && $release_num >= 0) or die "release_num argument must be numeric\n";
-d $target_path or die "target path $target_path does not exist\n";
-w $target_path or die "target path $target_path is not writable\n" unless $opt{G}; #< do not test dir writability if -G passed, since we won't actually be writing to this dir if we're just writing gbrowse confs.

# set the root dir to our target directory so the object will make it in the right place#
CXGN::ITAG::Release->releases_root_dir($target_path);

# make our new release object
my $release = CXGN::ITAG::Release->new( releasenum => $release_num,
					pre   => $opt{P},
					devel => $opt{D},
				      );

# if -S was passed, just collect stats, write the readme, and exit
if( $opt{S} ) {
    die Dumper( collect_stats( $release ) );
}

# if -G was passed, we just write the gbrowse configs and exit
if( $opt{G} ) {
  print "writing GBrowse configs and exiting.\n";
  write_gbrowse_genomic( $release, 'verbose' );
  write_gbrowse_prot(    $release, 'verbose' );
  exit;
}

###### MAKE OUTPUT DIRECTORY AND OPEN OUTPUT FILES

dump_data( $release );

#now collect statistics about this release to use in the README file
my $stats = collect_stats( $release );

#write the readme file
print "writing readme and GBrowse configuration ...\n";
write_readme( $release, 0, $stats );
write_gbrowse_genomic( $release );
write_gbrowse_prot( $release );

#and make the tarfile
#print "making release tarfile ...\n";
#$release->make_tarfile();

#update the ITAG_current link in the target dir, if any
#print "updating symlink ...\n";
#update_itag_current_link($target_path);

#and exit
exit 0;

######################## SUBROUTINES ##################

# for each contig, write the appropriate lines to the output files
sub dump_data {
    my ( $release ) = @_;
    my $gen_files = $release->get_all_files;

    # now open our ITAG pipeline, at the right version and base directory
    my $pipe = CXGN::ITAG::Pipeline->open( basedir => $opt{d},
                                           ( $opt{p} ? (version => $opt{p}) : () ),
                                          )
        or die "pipeline directory not found or not openable, do you need to specify a different -d or -p option?\n";

    # go through each batch, find the ones that have the sequences we
    # need, make sure the eugene analysis is done for it

    #first, figure out which batch number we're including
    my @batchnums = $pipe->list_batches; #< list all batches
    if ( $opt{b} ) {            #< if we got a -b option, parse it and use it to restrict the list
        @batchnums = restrict_batch_numbers( $opt{b}, @batchnums );
    }

    # convenient hash of all the analysis objects
    my %a =
        map {
            my $a = $pipe->analysis($_)
                or die "'$_' analysis not found in pipeline version ".$pipe->version.". Aborting.\n";
            ( $_ => {name => $_,  obj => $a} )
        } $pipe->list_analyses;

    # flag which analyses are required to be done in order for a batch to
    # be included in the release
    my @required_analyses = (
        #( map $_->{name}, grep {$_->{name} =~ /^blastp_/i} values %a), #< all the blastp analyses are required
        'renaming',
        'trnascanse',
        'sgn_markers',
        'infernal',
        'transcripts_sol',
        'transcripts_tomato',
       );
    foreach (@required_analyses) {
        $a{$_} or die "unknown analysis $_";
        $a{$_}->{required} = 1;
    }

    #get a list of the contigs we'll be using, including their batch numbers
    print "finding sequences ...\n";
    my $latest_seqs = list_latest_seqs( $pipe, \@batchnums, \%a );
    print "releasing annotations for ".@$latest_seqs." sequences.\n";

    foreach my $ctg (@$latest_seqs) {
        #check that all the files we need are readable
        -r or die "$_ file not readable\n" foreach ( $a{renaming}->{obj}->files_for_seq( $ctg->{batch}, $ctg->{name} ),
                                                     $a{seq}->{obj}->files_for_seq(    $ctg->{batch}, $ctg->{name} ),
                                                    );
    }

    # move build directory to <build>.prev if present, and delete any old .prev
    if(-d $release->dir ) {
        my $d = $release->dir;

        if( $opt{D} ) {
            system "rm -rf $d.prev";
            die if $CHILD_ERROR;
            system "mv $d $d.prev";
            die if $CHILD_ERROR;
        } else {
            die "$d already exists.  I won't overwrite it!\n";
        }
    }

    $release->mkdir();
    -d $release->dir or die "$! creating build directory ".$release->dir;
    -w $release->dir or die "build directory ".$release->dir." is not writable\n";

    my $gen_fh    = open_data_files( $gen_files );

    #put the header and feature-ontology in all the gff3 files
    my $sofa_url = $pipe->feature_ontology_url;

    #TODO: check that the sofa URL is still accessible
    $gen_fh->{$_}->print( "##gff-version 3\n##feature-ontology $sofa_url\n" )
        foreach grep $gen_files->{$_}->{type} eq 'gff3', keys %{$gen_fh};

    print "dumping data...\n";
    my $dump_count = 0;
    foreach my $ctg (@$latest_seqs) {
        my $ctg_name = $ctg->{name};
	### contig: $ctg_name
	print '.' unless ++$dump_count % 100;

        #get the names of all the files we need
        my (undef,$eug_prot,$eug_cds,$eug_cdna) = $a{renaming}->{obj}->files_for_seq( $ctg->{batch}, $ctg->{name} );
        my ($seq,$ctg_agp)                      =      $a{seq}->{obj}->files_for_seq( $ctg->{batch}, $ctg->{name} );

        copy_or_die( $seq,      $gen_fh->{genomic_fasta} );
        copy_or_die( $eug_cds,  $gen_fh->{cds_fasta}     );
        copy_or_die( $eug_cdna, $gen_fh->{cdna_fasta}    );

        # using this sequence's AGP file, make the sol -> gb ID mapping
        # file, and write some sequence-region lines
        if ( identifier_namespace( $ctg_name ) eq 'tomato_bac_contig' ) {
            process_agp( $ctg_name, $ctg_agp, $gen_fh );
        }

        # fetch the GO terms for each mrna if available
        my $go_terms = {};
        my ( $go_tabular ) = $a{go}->{obj}->files_for_seq( $ctg->{batch}, $ctg->{name} );
        if( ! -f $go_tabular ) {
            report_file_not_available( $ctg, $a{go}, $go_tabular );
        }
        elsif( -s $go_tabular ) {
            $go_terms = get_go_terms_for_mrnas( $go_tabular );
	    ### go terms: $go_terms
        }

        # fetch our human_readable_description if available
        my $human_readable_descriptions = {};
        my ( $prot_w_desc ) = $a{human_readable_description}->{obj}->files_for_seq( $ctg->{batch}, $ctg->{name} );
        if ( -f $prot_w_desc ) {
            $human_readable_descriptions = extract_human_readable_desc_strings( $prot_w_desc );
            format_deflines( $eug_prot, $human_readable_descriptions, $gen_fh->{protein_fasta} );
        } else {
            copy_or_die( $eug_prot, $gen_fh->{protein_fasta} );
        }

        # gff3-escape all the description strings (if any)
        $_ = gff3_escape( $_ )
            for values %$human_readable_descriptions;

        my @gff3_dump_specs =
            (
                { analyses => 'renaming',
                  gff3_output_spec_index => 0,
                  release_files => [qw| combi_genomic_gff3  models_gff3 |],
                  alter_lines_with => sub {
                      add_functional_annotations( shift, $go_terms, $human_readable_descriptions )
                  },
                  errors_fatal => 1,
                },
                { analyses => [qw[ geneid_tomato
                                   genemark_ath
                                   genemark_tom
                                   glimmerhmm_ath
                                   glimmerhmm_tomato
                                   augustus
                                 ]
                              ],
                  gff3_output_spec_index => 0,
                  release_files => [qw| combi_genomic_gff3 genefinders_gff3 |],
                },
                { analyses => 'trnascanse',
                  gff3_output_spec_index => 1,
                  release_files => [qw| combi_genomic_gff3 genefinders_gff3 |],
                },
                { analyses => [qw[
                                  sgn_markers
                                  sgn_loci
                                 ]
                              ],
                  gff3_output_spec_index => 1,
                  release_files => [qw| combi_genomic_gff3 sgn_data_gff3 |],
              },
                { analyses => 'sgn_unigenes',
                  gff3_output_spec_index => 1,
                  release_files => [qw| combi_genomic_gff3 cdna_algn_gff3 sgn_data_gff3|],
              },
                { analyses => 'infernal',
                  gff3_output_spec_index => 0,
                  release_files => 'combi_genomic_gff3',
              },
                { analyses => qr/^transcripts_/i,
                  gff3_output_spec_index => 0,
                  release_files => [qw| combi_genomic_gff3 cdna_algn_gff3 |],
              },
                { analyses => qr/^blastp_/i,
                  gff3_output_spec_index => 1,
                  release_files => 'functional_prot_gff3',
                  alter_lines_with => \&fix_blastp_gff3,
              },
                { analyses => 'interpro',
                  gff3_output_spec_index => 0,
                  alter_lines_with => \&fix_interpro_gff3,
                  release_files => 'functional_prot_gff3',
              },
                { analyses => 'rpsblast',
                  gff3_output_spec_index => 0,
                  alter_lines_with => \&rpsblast_m8_to_gff3,
                  release_files => 'functional_prot_gff3',
              },
               );

        # validate the above dumpspecs so we aren't surprised later with
        # errors after we have already been dumping for a while
        validate_dumpspec(\%a,$gen_fh,$_) for @gff3_dump_specs;

        # now do all the actual dumping, since we are now reasonably certain
        # that it will succeed
        foreach my $dump (@gff3_dump_specs) {
            foreach my $arecord ( find_analyses( \%a, $dump->{analyses} ) ) {


                my $index = $dump->{gff3_output_spec_index};
                my @filehandles = @{$gen_fh}{flatten $dump->{release_files}};

                my $sub = $dump->{alter_lines_with};

                my $file = ($arecord->{obj}->files_for_seq($ctg->{batch},$ctg->{name}))[$index];
                #print STDERR "dumping $arecord->{name}:$file to ".join(' ',flatten $dump->{release_files})."\n";
                unless ( $dump->{errors_fatal} || -f $file ) {
                    report_file_not_available( $ctg, $arecord, $file );
                } else {
                    copy_gff3_or_die( $arecord->{name}, $sub, $file, @filehandles);
                }
            }
        }
    }
    print "finished dumping data.\n";

    #close all the data files
    close_all( $gen_fh );

    # warn about any files that were not available
    warn_not_available_summary();

    # post-process all the GFF3 files in place to clean them up, sort them, add sync marks, etc.
    postprocess_gff3(  $gen_files->{$_->{seq_type}.'_fasta'}->{file},  $_->{file}  )
        foreach grep $_->{type} eq 'gff3', values %$gen_files ;
}

# make sure a given dumpspec is valid
sub validate_dumpspec {
  my ($a,$gen_fh,$d) = @_;
  ref $d eq 'HASH' or die "dumpspec must be a hashref!\n";

  my @required_keys = qw| analyses gff3_output_spec_index release_files |;
  my @optional_keys = qw| alter_lines_with errors_fatal |;
  my %valid_dump_keys = map {$_ => 1} @required_keys, @optional_keys;

  $valid_dump_keys{$_} or die "invalid key '$_' in dumpspec\ns"
    for keys %$d;

  $a->{$_} or die "invalid analysis name '$_' in dumpspec\n"
    for aname_flatten($a, $d->{analyses});

  my %legal_file_idents = map { $_ => 1 } keys %$gen_fh;
  $legal_file_idents{$_} or die "invalid file identifier '$_': see `perldoc CXGN::ITAG::Release`\n"
    for flatten $d->{release_files};

}

# acts like a super-flatten that also does regexp searches against the
# analysis names
sub aname_flatten {
  my ($a, $spec) = @_;
  return map {
    my $s = $_;
    ref $s eq 'Regexp' ? ( grep {$_ =~ $s} keys %$a )
                       : $s
  } flatten $spec;
}

# just uses the output of aname_flatten to return a hash-slice from %$a
sub find_analyses {
  my ( $a, $spec ) = @_;
  return @{$a}{ aname_flatten( $a, $spec) };
}

#parse the contig's AGP file and use its information
sub process_agp {
  my ( $ctg_name, $ctg_agp, $gen_fh ) = @_;

  my @ctg_tab_spec_sol =
#    my @ctg_tab_spec_gbacc =
      ($ctg_name); #< start assembling the line for the contig_members.tab file

  my $parsed_agp = agp_parse( $ctg_agp );

  #find the length of the contig
  my $ctg_length = max map $_->{oend}, @$parsed_agp;
  $_->print( "##sequence-region $ctg_name 1 $ctg_length\n" )
    foreach @{$gen_fh}{qw/contig_gff3 combi_genomic_gff3/};

  #go through each agp line, and make gff3 for it
  foreach my $line (@$parsed_agp) {
    next if $line->{comment};	#< skip comments and gaps
    $line->{is_gap} and die "$ctg_agp:$INPUT_LINE_NUMBER: there should not be any gaps in contig AGP files!";

    #     if( identifier_namespace( $line->{ident} ) eq 'bac_sequence' ) {
    #       #if it's a sol-style bac sequence ident, look up its genbank name
    #       if(my $gbacc = seq_name_to_genbank_acc( $line->{ident} ) ){
    # 	warn "translating $line->{ident} to $gbacc\n";
    # 	$line->{ident} = $gbacc;
    #       }
    #     }

    #transform identifiers if necessary
    my $ident_namespace = identifier_namespace( $line->{ident} ) || '';
    if ( $ident_namespace eq 'genbank_accession' ) {
      #if it's a sol-style bac sequence ident, look up its genbank name
      if (my $solid = genbank_acc_to_seq_name( $line->{ident}, get_dbh() ) ) {
	warn "$ctg_name spec: translating genbank accession $line->{ident} to $solid\n";
	$gen_fh->{sol_gb_mapping}->print( "$solid\t$line->{ident}\n" );
	$line->{ident} = $solid;
      } else {
	die "could not find SOL ID for genbank accession '$line->{ident}'!\n";
      }
    } elsif ( $ident_namespace eq 'bac_sequence' ) {
      if (my $gbacc = seq_name_to_genbank_acc( $line->{ident}, get_dbh() ) ) {
	$gen_fh->{sol_gb_mapping}->print( "$line->{ident}\t$gbacc\n" );
      } else {
	$gen_fh->{sol_gb_mapping}->print( "$line->{ident} NO_ACCESSION_REGISTERED\n" );
      }
    }

    #push this info onto the contig's member listing
    push @ctg_tab_spec_sol, $line->{ident};
    #push @ctg_tab_spec_gb,

    #and write out a GFF3 line for this component in the assembly
    foreach ( @{$gen_fh}{qw/contig_gff3 combi_genomic_gff3/} ) {
        $_->print( join("\t",
                        $ctg_name,
                        'ITAG',
                        'clone',
                        $line->{ostart},
                        $line->{oend},
                        qw/. + ./,
                        "Name=$line->{ident};Target=$line->{ident} $line->{cstart} $line->{cend} +"
                       )
                   ."\n###\n" )
    }

  }
  #print the member listing for this contig
  $gen_fh->{contig_members_sol}->print( join("\t",@ctg_tab_spec_sol),"\n" );
  #$gen_fh->{contig_members_gbacc}->print( join("\t",@ctg_tab_spec_gbacc),"\n" );
}

# input is a line of RPSBLAST -m 8 format, output is a line of gff3,
# or empty string.  converts rpsblast m8 into gff3
sub rpsblast_m8_to_gff3 {
  my $line = shift;

  my ($qname,$hname, $percent_id, $hsp_len, $mismatches,$gapsm,
      $qstart,$qend,$hstart,$hend,$evalue,$bits) = split /\s+/,$line;
  return '' unless $bits;

  my @gff3_fields = ( $qname,
		      'ITAG_rpsblast',
		      'match',
		      $qstart,
		      $qend,
		      $bits,
		      '+',
		      '.',
		      join ';',
		      (
		       map {join '=',@$_}
		       [ Name   => gff3_escape($hname) ],
		       [ Target => gff3_escape($hname)." $hstart $hend +" ],
		       [ rpsblast_evalue => $evalue],
		       [ rpsblast_pct_id => $percent_id],
		       [ rpsblast_mismatches => $mismatches],
		       [ rpsblast_gaps => $gapsm],
		       [ rpsblast_bits => $bits],
		       [ rpsblast_hsp_len => $hsp_len],
		      )
		    );
  return join("\t",@gff3_fields)."\n";
}

# input and output is a single line of gff3.  just munges ITAG
# blastp_* gff3 slightly.
sub fix_blastp_gff3 {
  $_[0] =~ s/(?<=\t|;)(\w+=)\s+/$1/g;
  $_[0] =~ s/name=(\S+)\s+/Name=$1;Note=/;
  #$_[0] =~ s/Name=[a-z]{2,3}\|([^\|;]+)\|/Name=$1/;
  return $_[0];
}
;

# input and output is a single line of gff3.  just munges ITAG
# interpro gff3 slightly
sub fix_interpro_gff3 {
  #take out the polypeptide features, and the 'Parent' references to them,
  #leaving on the bare polypeptide regions
  return '' unless $_[0]=~ /polypeptide_region/;
  $_[0] =~ s/ID=[^-]+-([^-]+)-\d+;/Name=$1;/;
  $_[0] =~ s/description=/Note=/;
  $_[0] =~ s/Parent=[^;]+;//;
  return $_[0];
}

sub gff3_escape {
  uri_escape(shift,'\n\r\t;=%&,\x00-\x1f\x7f-\xff');
}


# get all the gene models from renaming and dump them into the combi and models files
sub add_functional_annotations {
  my ($gff3_line, $ontology_terms, $descriptions ) = @_;
  # currently only operate on mRNA features
  my ($type,$id) =
      $gff3_line =~ m(  \t
                        (?:mRNA)  # feature type
                        \t        # flanked by tabs
                        .+        # and then whatever
                        ID=(mRNA):([^;\n]+) # and then an ID attribute with an mRNA type
                     )x
     or return $gff3_line;

  chomp $gff3_line;
  ### Adding functional annotations...
  ### type: $type
  ### id: $id

  # add the human-readable desc if present
  if( my $desc_string = $descriptions->{$id}) {
      # add it to the gff3 string
      ### desc: $desc_string
      $gff3_line .= ";functional_description=".uri_escape( $desc_string );
  }

  # add ontology terms if present
  if( my $terms = $ontology_terms->{$id} ) {
      ### got terms: $terms
      $gff3_line .= ";Ontology_term=$_" for @$terms;
  }

  $gff3_line .= "\n";

  return $gff3_line;
}

# given a fasta file containing the functional description in the
# description line, return a hashref of { gene_name => 'desc string', ... }
sub extract_human_readable_desc_strings {
  my ( $prot_w_desc ) = @_;

  return {} unless -s $prot_w_desc > 3;

  my %descriptions;
  open my $fa, '<', $prot_w_desc or die "$! reading file '$prot_w_desc'";
  while( my $l = <$fa> ) {
    next unless $l =~ /^\s*>/;
    chomp $l;
    my ($genemodel_name, $functional_desc) = $l =~ /^\s*>\s*(\S+)\s+(.+)/
      or die "$prot_w_desc:$INPUT_LINE_NUMBER: cannot parse fasta file $prot_w_desc.  Is this a valid fasta file?";
    length($functional_desc) > 10 or die "functional description '$functional_desc' too short";
    $descriptions{ $genemodel_name } = $functional_desc;
  }
  %descriptions or die "cannot extract human readable description string from file '$prot_w_desc'.  Is this a valid fasta file?";
  return \%descriptions;
}


sub update_itag_current_link {
  my ($target_path) = @_;

  # set the root dir to our target directory so the object will make it in the right place
  CXGN::ITAG::Release->releases_root_dir($target_path);

  my @releases = CXGN::ITAG::Release->find;

  #if there is an official release, always use that as current
  my ($curr_rel) = grep !$_->is_devel_release && !$_->is_pre_release, @releases;
  #otherwise, just use the most recent devel or pre release
  $curr_rel ||= $releases[0];

  my $curr_link = File::Spec->catfile($target_path,'ITAG_current');
  unlink($curr_link);
  symlink($curr_rel->dir,$curr_link) or die "$! linking ".$curr_rel->dir." -> $curr_link";
}

#take a hashref of name => filename, return a new hashref of  name => opened filehandle
#skips files tagged with readme*
sub open_data_files {
  my $gen_files = shift;
  my %new;
  while( my ($tag,$rec) = each %$gen_files) {
    next if $tag =~ /^readme/i;
    $new{$tag} = IO::File->new("> $rec->{file}")
      or die "$! opening $rec->{file} for writing\n";
  }
  lock_hash(%new);
  return \%new;
}
sub close_all {
  foreach my $fh (values %{+shift}) {
    $fh->close;
  }
}

#go through all the files in the release, collect statistics about them,
#return it as a locked hashref
sub collect_stats {
  my ($release) = @_;
  my $gen_files = $release->get_all_files;

  print "collecting statistics...\n";

  my %stats = map {$_=>''} qw(
			      gene_cnt
			      gene_model_cnt
			      gene_mrna_counts
                              gene_length_avg
			      genes_with_splice_variants_pct
			      genes_with_splice_variants_cnt
                              genes_with_ontology_terms
			      genome_coverage_pct
			      mapped_ests_cnt
			      protein_coding_with_cdna_or_est_cnt
			      protein_coding_without_cdna_or_est_cnt
			      protein_coding_with_prot_cnt
			      protein_coding_without_prot_cnt
			      gene_model_length_classes
			      loci_similar_to_sgn_loci
			      loci_similar_to_unch_prots_cnt
			      loci_no_prot_hits_cnt
			      loci_no_cdna_est_evidence_cnt
                              unique_ontology_terms
			     );
  lock_keys(%stats);

  #open the aggregated GFF3 file
  open my $combi_in, '<', $gen_files->{combi_genomic_gff3}->{file} or die "$! reading combi gff3 file\n";
  my $gene_length_accum = 0;
  my %ontology_terms_seen;
  while( my $line = <$combi_in> ) {
    chomp $line;
    ### line: $line
    my @f = split /\t/, $line, 9;

    my $src = lc $f[1];
    my $type = $f[2];
    ### src: $src
    ### type: $type
    if( $src eq 'itag_renaming' ) {
      if( $type eq 'gene' ) {
	$stats{gene_cnt}++;
        $gene_length_accum +=  $f[4]-$f[3]+1;
	### length: $f[4]-$f[3]+1
      } elsif( $type eq 'mRNA' ) {
	$stats{gene_model_cnt}++;
        $stats{gene_mrna_counts} ||= {};
	my ($parent) = $line =~ /Parent=gene:([^;\n]+)/ or die "cannot parse parent from gff3 line:\n$line\n";
	### parent: $parent
        $stats{gene_mrna_counts}{$parent}++;
        if( my @terms = $line =~ /Ontology_term=([^;\n]+)/g ) {
            $stats{genes_with_ontology_terms}++;
            $ontology_terms_seen{$_} = 1 for @terms;
	    ### terms: @terms
        }
      }
    }
    elsif( $src =~ /^itag_transcripts_/i ) {
      if( $type eq 'match' ) {
        my @tgts = map [split], $line =~ /Target=([^;\n]+)/g;
	### targets: @tgts
	our %stats_seen_est;
	unless($stats_seen_est{$tgts[0][0]}++) {
	  $stats{mapped_ests_cnt}++;
	}
      }
    }
  }
  close $combi_in;

  # calculate average gene length
  $stats{gene_length_avg} = sprintf( '%0.0f', $gene_length_accum / $stats{gene_cnt} );

  # calculate number of unique ontology terms
  $stats{unique_ontology_terms} = scalar keys %ontology_terms_seen;

  ## aggregate the splice variant statistics
  my $variants = delete $stats{gene_mrna_counts};
  # currently just counting how many genes have been annotated with splice variants
  $stats{genes_with_splice_variants_cnt} = scalar grep $_ > 1, values %$variants;
  $stats{genes_with_splice_variants_pct} = sprintf( '%0.1f', 100 * $stats{genes_with_splice_variants_cnt}/$stats{gene_cnt} );

  #estimate the coverage percentage by adding up the lengths of the
  #sequence-regions in the gene models file and dividing that by the
  #estimated genome size
  my $estimated_genome_size = 240_000_000; #< in bases
  my $covered_bases = 0;
  open my $gm, $gen_files->{models_gff3}->{file}
    or die "$! reading $gen_files->{models_gff3}->{file}";
  while( <$gm> ) {
    if( /##\s*sequence-region\s+\S+\s+(\d+)\s+(\d+)/ ) {
      $covered_bases += $2 - $1 + 1;
    }
  }
  close $gm;

  $stats{genome_coverage_pct} = sprintf('%0.0f',$covered_bases/$estimated_genome_size*100);

  # get statistics about gene models from the gene description codes
  # generated by EuGene in the deflines of the protein and CDS fasta
  # files

  my %gene_descriptions;
  open my $deflines, "grep '>' $gen_files->{protein_fasta}->{file} |"
    or die "$! running grep on $gen_files->{protein_fasta}->{file}";
  while(my $line = <$deflines>) {
    $line =~ s/^\s*>//; #<trim off the beginning symbol
    my ($ident,$def) = split /\s+/,$line,2;
    if(my $desc = parse_gene_description($def)) {
      $gene_descriptions{$ident} = $desc;
      if( $desc->{cdna_complete_coverage} || $desc->{cds_from_cdna_aln} ) {
	$stats{protein_coding_with_cdna_or_est_cnt}++;
      } else {
	$stats{protein_coding_without_cdna_or_est_cnt}++;
      }

      if( $desc->{gene_model_from_prot_aln} ) {
	$stats{protein_coding_with_prot_cnt}++;
      } else {
	$stats{protein_coding_without_prot_cnt}++;
      }

      $stats{'gene_model_length_classes'} ||= {};
      my $lc = $desc->{length_class};
      $lc = 'none' unless defined $lc;
      $stats{'gene_model_length_classes'}{$lc}++;
    } else {
      chomp $line;
      die "ERROR: no parsable gene description found in defline $line\n";
      next;
    }
  }
  close $deflines;

  # now go through the functional annotation and gather statistics on
  # the similarity of loci to various things
#   my $functional_in = Bio::FeatureIO->new( -format => 'gff', -version => 3, -file => $gen_files->{functional_prot_gff3}->{file} );
#   while( my $feat = $functional_in->next_feature ) {
#   }
  lock_hash(%stats);

  return \%stats;
}

#generate a block of text describing each file in the release
sub file_descriptions {
  my ($release_info,$gen_files) = @_;

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

#given a gff3 file, sort it and clean up its pragmas
sub postprocess_gff3 {
  my ($seqs_file,$gff3_file) = @_;

  print "postprocessing gff3 $gff3_file with reference sequences $seqs_file\n";

  -f $seqs_file or die "seqs file '$seqs_file' not found!";
  -f $gff3_file or die "gff3 file '$gff3_file' not found!";

  clear_uniq_ids();

  #open a temp file for our output
  my ($temp_out_fh,$temp_out_file) = tempfile(File::Spec->catfile(File::Spec->tmpdir,"$FindBin::Script-gff3-postproc-XXXXXX"), UNLINK => 1);

  #sort the gff3 file by ref sequence and start coordinate into the temp file,
  #and remove duplicates
  open my $sg, "sort -k 1,1 -k 4,4g -s $gff3_file | gff3_reformat.pl -i -S $seqs_file | "
    or die "$! running sort on file '$gff3_file'";

  #print the sorted output into the temp file, massaging the
  #directives: buffering some, dropping others
  my @sequence_regions;
  my @other_pragmas;
  while ( my $line = <$sg> ) {
    chomp $line;
     if ( $line =~ /##\s*(.+)$/ ) {
      #it's a directive
      my $directive = $1;
      if ( $directive =~ /^gff-version/i  || $directive =~ /^#/  ) {
	#ignore gff-version and synchronization pragmas
      } elsif ( $directive =~ /^sequence-region\s+(.+)/i ) {
	print $temp_out_fh "$line\n";
      } elsif ( $directive =~ /^FASTA/ ) {
	die "ERROR: gff3 postprocessing does not currently support gff3 FASTA sections, and '$gff3_file' contains one\n";
      } elsif( grep {index($directive,$_) == 0} qw/feature-ontology attribute-ontology source-ontology/ ) {
	push @other_pragmas,$directive;
      }
    } elsif($line =~ /^\s*#/) {
      #get rid of comments
    } elsif( $line =~ /\S/ ) { #< get rid of blank lines
      my @fields = split /\t/,$line;
      @fields == 9 or die "$gff3_file:$INPUT_LINE_NUMBER: invalid gff3 line:\n$line";
      #print $temp_out_fh "weird line:\n$line\n" unless @fields == 9;
      $fields[1] = 'ITAG' if $fields[1] eq 'ITAG_ITAG';

      #uniqify the identifiers
      $fields[8] =~ s/(?<=ID=)([^;]+)/uniq_id($1)/e;
      $fields[8] =~ s/(?<=Parent=)([^;]+)/last_uniq_id($1)/e;

      print $temp_out_fh join "\t",@fields;
      print $temp_out_fh "\n";
    }
  }
  close $temp_out_fh;

  #now open the original file, print the pragmas we found, and copy the rest into it
  open my $orig, ">$gff3_file" or die "$! writing $gff3_file";
  print $orig "##gff-version 3\n";
  #print $orig "##sequence-region $_\n" foreach distinct sort @sequence_regions;
  print $orig "##$_\n" foreach @other_pragmas;
  $orig->flush();
  copy_or_die( $temp_out_file, $orig );
  close $orig;
}

{
    my %id_ctr;
    sub clear_uniq_ids() {
        %id_ctr = ();
    }
    sub uniq_id {
        my $id = shift;
        my $cnt = $id_ctr{$id}++;
        if ($cnt) {
            return "$id-i$cnt";
        }
        return $id;
    }
    sub last_uniq_id {
        my $id = shift;
        no warnings 'uninitialized';
        my $cnt = $id_ctr{$id}-1;
        if ($cnt > 0) {
            return "$id-i$cnt";
        }
        return $id;
    }
}

#given a line of text containing a gene description somewhere, parse
#it and return a hashref of its contents, or nothing if no description
#was found
sub parse_gene_description {
  my ($text) = @_;

  my %stuff;

  @stuff{qw/year cdna_complete_coverage gene_model_from_prot_aln cds_from_cdna_aln program length_class/}
    = $text =~ / \b (\d\d) F (\d) H (\d) E (\d) I ([A-Z]{2}) (L (\d))? \b /x
      or return;

  my %program_map = ( EG => 'EuGene' );
  $stuff{program} = $program_map{$stuff{program}} or die "unknown program code '$stuff{program}' found in gene description";

  lock_hash(%stuff);

  return \%stuff;
}

sub format_deflines {
  my ($prot_seq_file, $descriptions, $out_fh) = @_;

  my $prot_seqs = Bio::SeqIO->new(-file => $prot_seq_file, -format => 'fasta');

  while( my $s = $prot_seqs->next_seq ) {
    my $desc = $descriptions->{ $s->display_id };
    if( $desc ) {
      $desc =~ s/([\\"])/\$1/g;
      $desc = qq|"$desc"|;
    }

    my ($contigname, $evidence_code, $coords, $timestamp ) = split /\s+/, $s->desc
      or die "could not parse description for ".$s->display_id.": (".$s->desc.")";

    $s->desc( "genomic_reference:$contigname evidence_code:$evidence_code gene_region:$coords"
	      .($desc ? " functional_description:$desc" : '')
	    );

    my $line = '>'.$s->display_id.' '.$s->desc."\n".$s->seq."\n";
    $out_fh->print( $line );
  }
}

#write the readme file
sub write_readme {
  my ( $release, $ncbi_tax_id, $stats ) = @_;

  my $gen_files = $release->get_all_files;

  $ncbi_tax_id ||= 0; #or die "must pass ncbi tax id as second argument to write_readme()\n";
  my $organism = 'Tomato';
  my $project_name = 'International Tomato Genome Annotation';
  my $project_acronym = 'ITAG';

  $stats ||= {};
  my %fmt_stats = %$stats;
  $_ = commify_number($_ || 0) foreach values %fmt_stats;

  # format the splice variants text
  if( $fmt_stats{genes_with_splice_variants_cnt} ) {
      $fmt_stats{genes_with_splice_variants_pct} = " ($fmt_stats{genes_with_splice_variants_pct}%)";
  } else {
      $fmt_stats{genes_with_splice_variants_cnt} = 'no';
      $fmt_stats{genes_with_splice_variants_pct} = '';
  }

  lock_hash(%fmt_stats);

  my $file_descriptions = file_descriptions($release,$gen_files);

  #use Data::Dumper;
  #warn 'stats are '.Dumper(\%fmt_stats);

# The $release_tag release is also available at NCBI:

# http://www.ncbi.nlm.nih.gov/mapview/map_search.cgi?taxid=$ncbi_tax_id

# Datasets are also available from SGN's bulk download tool (paste in or
# upload a list of identifiers and download the corresponding data):

# http://www.solgenomics.net/bulk/input.pl

  #PARAGRAPH ABOUT EST COVERAGE AND PROTEIN SIMILARITY
# $fmt_stats{loci_similar_to_unch_prots_cnt} loci have similarity only to uncharacterised proteins (i.e. hypothetical, predicted, unknown etc), $fmt_stats{loci_no_prot_hits_cnt} have no significant protein similarity to GenBank proteins, and of these $fmt_stats{loci_no_cdna_est_evidence_cnt} also have no supporting EST/cDNA evidence and may represent erroneous gene predictions.

  if( $fmt_stats{genome_coverage_pct} > 100 ) {
      $fmt_stats{genome_coverage_pct} = sprintf('%0.1f-fold',$fmt_stats{genome_coverage_pct}/100);
  } else {
      $fmt_stats{genome_coverage_pct} .= '%';
  }

  my $release_dirname = $release->dir_basename;
  my $release_tag     = $release->release_tag;

  my $date_str = POSIX::strftime( "%A %B %e, %Y", gmtime() );

  my $readme_text = <<EOT;
$release_tag $organism Genome release

The $project_name project ($project_acronym) is pleased to announce the release of the latest version of the official $organism genome annotation ($release_tag).  This set of release files was generated on $date_str.

The $release_tag release contains $fmt_stats{gene_cnt} genes in all, with $fmt_stats{gene_model_cnt} gene models. Average gene length is $fmt_stats{gene_length_avg} base pairs.  Currently, $fmt_stats{genes_with_splice_variants_cnt} genes$fmt_stats{genes_with_splice_variants_pct} have annotated splice variants.  The current $project_acronym annotation provides approximately $fmt_stats{genome_coverage_pct} of the $organism genome.  $fmt_stats{genes_with_ontology_terms} genes have ontology terms associated, with a total of $fmt_stats{unique_ontology_terms} different ontology terms represented.

This $project_acronym release has $fmt_stats{mapped_ests_cnt} cDNA and EST sequences mapped to the genome, resulting in $fmt_stats{protein_coding_with_cdna_or_est_cnt} protein coding genes derived at least partly from supporting cDNA and/or EST alignments and thus $fmt_stats{protein_coding_without_cdna_or_est_cnt} protein coding genes not utilizing transcript support.  With respect to protein homology, $fmt_stats{protein_coding_with_prot_cnt} gene models used homology to known proteins in their construction, while $fmt_stats{protein_coding_without_prot_cnt} did not.

Files included in this release:

$file_descriptions

Sequences and annotations can also be viewed and searched on SGN:

http://www.solgenomics.net/gbrowse/

The fully annotated chromosome sequences in GFF version 3 format, along with Fasta files of cDNA, CDS, genomic and protein sequences, and lists of genes are available from the SGN ftp site at:

ftp://ftp.solgenomics.net/tomato_genome/annotation/${release_dirname}/

For those who are not familiar with the relatively new GFF3 format,
the format specification can be found here:

http://www.sequenceontology.org/gff3.shtml

A graphic display of the $organism sequence and annotation can be viewed using SGN's genome browser. Browse the chromosomes, search for names or short sequences and view search hits on the whole genome, in a close-up view or on a nucleotide level:

http://www.solgenomics.net/gbrowse/

SGN's BLAST services have also been updated with the new datasets and are available from:

http://www.solgenomics.net/tools/blast/

$project_acronym is committed to the continual improvement of the $organism genome annotation and actively encourages the plant community to contact us with new data, corrections and suggestions.

Announcements of new releases, updates of data, tools, and other developments from $project_acronym can be found

http://www.solgenomics.net/

Please send comments or questions to:

itag\@sgn.cornell.edu

The $project_acronym Team

EOT

  #do some silliness to correct the word wrapping of the text so it
  #comes out looking nice no matter what
#  $readme_text =~ s/(?<=[^\n])\n(?=[^\n])/ /g;
  $readme_text = wrap_long_lines($readme_text);


  #and print it to the readme file
  open my $readme, '>',$gen_files->{readme}->{file}
    or die "$! opening $gen_files->{readme}->{file}";
  print $readme $readme_text;
  close $readme;

  return 1;
}

sub wrap_long_lines {
  my $text = shift;
  my $length = shift || 80;
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

# return a arrayref like
# [ { seq name => {num => ctg num, batchnum => batch number, batch => batch obj},
#    ...
# ]
sub list_latest_seqs {
  my ($pipe,$batchnums,$an) = @_;

  my @seqs;

 BATCH:
  foreach my $bn (@$batchnums) {

    my $batch = $pipe->batch($bn);

    my $skip_flag = 0;
    while( my ($aname,$a) = each %$an ) {
      #warn "checking done for $aname ".Dumper($aname,$a);
      my $status = $a->{obj}->status($batch);
      unless( !$a->{required} || $status eq 'done' ) {
	warn "Analysis '$aname' is not done in batch $bn.\n";
	$skip_flag = 1;
      }
    }
    if( $skip_flag ) {
      warn "Skipping batch $bn.\n";
      next BATCH;
    }

    foreach ( $batch->seqlist ) {
        push @seqs, { batch => $batch, name => $_ };
    }
  }

  return [ sort { $a->{name} cmp $b->{name} }
           @seqs
         ];
}


sub gbrowse_conf_head {
  my ($release,$type) = @_;

  my $desc = $release->text_description.": $type Annotations";
  $desc =~ s/\b([a-z]+)\b/ucfirst($1)/ge;
  my $relname = $release->release_name;

  return
q|[GENERAL]
description = |.$desc.q|
adaptor  = dbi::cxgn_pg
db_args  = -db |.$relname.'.'.lc($type).q|

plugins = AttributeHiliter

reference class = Sequence

# Web site configuration info
stylesheet  = /documents/gbrowse/gbrowse.css
buttons     = /documents/gbrowse/images/buttons
tmpimages   = /documents/tempfiles/gbrowse
js          = /documents/gbrowse/js
cachedir    = /documents/tempfiles/gbrowse
help        = /documents/gbrowse

# Default glyph settings
glyph       = generic
height      = 8
bgcolor     = cyan
fgcolor     = black
label density = 25
bump density  = 10000

# where to link to when user clicks in detailed view
link        = AUTO

# what image widths to offer
image widths  = 450 680 980 1550 2400

# default width of detailed view (pixels)
default width = 680
upload_tracks section = closed
display_settings section = closed

# max and default segment sizes for detailed view
max segment     = 500000
default segment = 50000

# zoom levels
zoom levels    = 100 200 1000 2000 5000 10000 20000 40000 100000 200000 500000 1000000

# colors of the overview, detailed map and key
overview bgcolor = white
detailed bgcolor = white
key bgcolor      = beige

language = en

# a footer
# Various places where you can insert your own HTML -- see configuration docs
html1 =
html2 =
html3 =
html4 =
#html4 = <div style="border: 1px solid red; padding: 4px"><b>Troubleshooting tip:</b> if the browser seems to be acting strangely, try hitting the <a class="reset_button" href="?reset=1"><b>[Reset]</b></a> link.</div>
html5 =
html6 =

### defaults examples sub, works for any of them
examples=sub
    {
        my ($dbgff) = @_;
        if($dbgff && $dbgff->can('features_db') ) {
          my $dbh = $dbgff->features_db;
          my $examples = $dbh->selectall_arrayref('select fref,seqlen from fdata join (select fref,max(foffset) as seqlen from fdna group by fref having max(foffset) <=500000 order by seqlen desc limit 70) as lens using(fref) group by fref,seqlen order by seqlen desc limit 4 offset 2');
          return map $_->[0],@$examples
        } else {
          return ();
        }
    }


### DAS CONFIGURATION ####
das mapmaster = SELF
|

}

sub write_gbrowse_genomic {
  my ($release, $verbose) = @_;

  my $type = 'genomic';

  my $file = $release->get_file_info("gbrowse_${type}_conf")->{file};
  my $protconf = $release->get_file_info("gbrowse_protein_conf")->{file};
  my $prot_link_base = "/gbrowse/gbrowse/".basename($protconf,'.conf');

  print "writing $type conf to ' $file '\n" if $verbose;

  my $releasenum = $release->release_number;
  my $is_devel = $release->is_devel_release ? 1 : 0;
  my $is_pre = $release->is_pre_release ? 1 : 0;
  my $releasetag = $release->release_tag;

  open my $f,'>', $file or die "$! writing $file";

  print $f gbrowse_conf_head($release,$type);

  print $f q|
default features = genespan mrna cds tilingpath

aggregators =
	gappedmatch{match_part/match}
	eugtranscript{five_prime_UTR,CDS,three_prime_UTR}
	eugcds{CDS}
        italymatch{match}
# "automatic" classes to try when an unqualified identifier is given
automatic classes = Sequence clone match nucleotide_motif

### TRACK CONFIGURATION ####
# the remainder of the sections configure individual tracks

#track to show DNA and GC content
[DNA]
glyph          = dna
global feature = 1
height         = 40
do_gc          = 1
fgcolor        = red
category       = General
axis_color     = blue
strand         = both
key            = DNA/GC Content
citation       = This track displays a GC content graph of the reference sequence at low magnifications and the DNA sequence itself at higher magnifications.

#track to show a 6-frame translation of the sequence
[Translation]
glyph          = translation
global feature = 1
height         = 40
fgcolor        = purple
category       = General
start_codons   = 0
stop_codons    = 1
translation    = 6frame
key            = 6-frame translation
citation       = This track displays a six-frame translation of the reference DNA sequence.

[genespan]
feature      = gene:ITAG_renaming
key          = Gene span
glyph        = segments
fgcolor      = black
bgcolor      = darkorange
stranded     = 1
category     = Gene Model features
das category = transcription
strand_arrow = 1
height       = 10
citation     = This track shows the spans of gene models annotated by EuGene, the integrative gene predictor for ITAG. (see <a href="http://www.ab.wur.nl/TomatoWiki/AnEugene000">EuGene on ITAG Wiki</a>, <a href="http://www.inra.fr/internet/Departements/MIA/T//EuGene/index.html">EuGene main page</a>)

[mrna]
feature      = eugtranscript:ITAG_renaming
key          = mRNA
glyph        = processed_transcript
fgcolor      = black
bgcolor      = goldenrod
stranded     = 1
description  = sub { use CXGN::Page::FormattingHelpers;  CXGN::Page::FormattingHelpers::truncate_string((shift->attributes('functional_description'))[0], 20, '...') }
font2color   = blue
title        = sub { (shift->attributes('functional_description'))[0] }
category     = Gene Model features
das category = transcription
strand_arrow = 1
height       = 10
citation     = This track shows the mRNAs for gene models annotated by EuGene, the integrative gene predictor for ITAG. (see <a href="http://www.ab.wur.nl/TomatoWiki/AnEugene000">EuGene on ITAG Wiki</a>, <a href="http://www.inra.fr/internet/Departements/MIA/T//EuGene/index.html">EuGene main page</a>)

[cds]
feature      = eugcds:ITAG_renaming
key          = CDS - click models to browse protein domains
title        = sub { my $n = shift->display_name; $n =~ s/(mRNA\|CDS)://;  $n =~ s/(\.\d+)\.\d+$/$1/; "click to browse annotations on $n protein product"}
glyph        = segments
fgcolor      = black
bgcolor      = yellow
label        = sub { my $n = shift->display_name; $n =~ s/mRNA:/CDS:/; $n }
stranded     = 1
category     = Gene Model features
das category = transcription
strand_arrow = 1
link         = sub { my $n = shift->display_name; $n =~ s/(mRNA\|CDS)://; $n =~ s/(\.\d+)\.\d+$/$1/; "|.$prot_link_base.q|/?name=$n" }
height       = 10
citation     = This track shows the CDS sequences for gene models annotated by EuGene, the integrative gene predictor for ITAG.  Provided by Stephane Rombauts and Jeffrey Fawcett, <a class="http" href="http://bioinformatics.psb.ugent.be/">Bioinformatics and Evolutionary Genomics</a>, <a class="http" href="http://www.psb.ugent.be">Plant Systems Biology, VIB, Ghent University</a>.  (see also <a href="http://www.ab.wur.nl/TomatoWiki/AnEugene000">EuGene on ITAG Wiki</a>, <a href="http://www.inra.fr/internet/Departements/MIA/T//EuGene/index.html">EuGene main page</a>)

[cdna_tom]
feature      = italymatch:ITAG_transcripts_tomato
key          = ESTs and cDNAs - Tomato
category     = Genome data and reagents
glyph        = segments
stranded     = 1
link         = http://www.ebi.ac.uk/ebisearch/search.ebi?db=nucleotideSequences&t=$name
citation     = This track shows regions of similarity with EST and other cDNA sequences from Tomato. Provided by <a href="http://cab.unina.it/index2.php">CAB group</a> at <a href="http://www.unina.it">UNINA</a>. (see <a href="http://www.ab.wur.nl/TomatoWiki/AnEST000">ITAG Wiki</a>)

[cdna_sol]
feature      = italymatch:ITAG_transcripts_sol
key          = ESTs and cDNAs - Other Solanaceae
category     = Genome data and reagents
glyph        = segments
stranded     = 1
link         = http://www.ebi.ac.uk/ebisearch/search.ebi?db=nucleotideSequences&t=$name
citation     = This track shows regions of similarity with EST and other cDNA sequences taken from species in the Solanaceae other than Tomato. Provided by <a href="http://cab.unina.it/index2.php">CAB group</a> at <a href="http://www.unina.it">UNINA</a>. (see <a href="http://www.ab.wur.nl/TomatoWiki/AnEST000">ITAG Wiki</a>)

[sgn_unigenes]
feature      = match:ITAG_sgn_unigenes
key          = SGN Unigenes
category     = Genome data and reagents
glyph        = segments
stranded     = 1
link         = http://solgenomics.net/search/quick_search.pl?term=$name
citation     = This track shows regions of similarity with SGN Unigene (SGN-U) sequences.  (see <a href="http://www.ab.wur.nl/TomatoWiki/AnSgnUnigenes000">ITAG Wiki</a>)

[augustus]
feature      = processed_transcript:
key          = AUGUSTUS (de novo, Tomato trained)
category     = Prediction Features
glyph        = processed_transcript
stranded     = 1
citation      = <i>De novo</i> gene predictions from the AUGUSTUS gene predictor, trained on Tomato.  Provided by SGN.  (see <a href="http://www.ab.wur.nl/TomatoWiki/AnAUGUSTUS000">ITAG Wiki</a>)

[geneid_tomato]
feature      = processed_transcript:ITAG_geneid_tomato
key          = GeneID (de novo, Tomato trained)
category     = Prediction Features
glyph        = processed_transcript
stranded     = 1
citation     = <i>De novo</i> predictions from <a href="http://genome.imim.es/software/geneid/">GeneID</a>, trained on Tomato.  Provided by Francisco Camara, <a href="http://genome.crg.es">Genome Bioinformatics Research Lab - Gene Prediction Group, Center for Genomic Regulation (CRG) - Spain</a>. (see also <a href="http://www.ab.wur.nl/TomatoWiki/AnGeneID000">ITAG Wiki</a>)

[genemark_ath]
feature      = processed_transcript:ITAG_genemark_ath
key          = GeneMark (de novo, Arabidopsis trained)
category     = Prediction Features
glyph        = processed_transcript
stranded     = 1
citation     = <i>De novo</i> predictions from <a href="http://exon.gatech.edu/GeneMark/">GeneMark</a>, trained on Arabidopsis.  Provided by the <a href="http://mips.gsf.de/proj/plant/jsf/tomato/index.jsp">MIPS Tomato genome database</a>. (see also <a href="http://www.ab.wur.nl/TomatoWiki/AnGeneMark000">ITAG Wiki</a>)

[glimmerhmm_ath]
feature      = processed_transcript:ITAG_glimmerhmm_ath
key          = GlimmerHMM (de novo, Arabidopsis trained)
category     = Prediction Features
glyph        = processed_transcript
stranded     = 1
citation     = <i>De novo</i> predictions from <a href="http://www.genomics.jhu.edu/GlimmerHMM/">GlimmerHMM</a>, trained on Arabidopsis.  Provided by Erwin Datema at <a href="http://appliedbioinformatics.wur.nl/">WUR</a>. (see also <a href="http://www.ab.wur.nl/TomatoWiki/AnGlimmerHMM000">ITAG Wiki</a>)

[glimmerhmm_tomato]
feature      = processed_transcript:ITAG_glimmerhmm_tomato
key          = GlimmerHMM (de novo, Tomato trained)
category     = Prediction Features
glyph        = processed_transcript
stranded     = 1
citation     = <i>De novo</i> predictions from <a href="http://www.genomics.jhu.edu/GlimmerHMM/">GlimmerHMM</a>, trained on Tomato.  Provided by Erwin Datema at <a href="http://appliedbioinformatics.wur.nl/">WUR</a>. (see also <a href="http://www.ab.wur.nl/TomatoWiki/AnGlimmerHMM000">ITAG Wiki</a>)

[infernal]
feature      = snoRNA:ITAG_infernal U6atac_snRNA:ITAG_infernal miRNA:ITAG_infernal antisense_RNA:ITAG_infernal snRNA:ITAG_infernal 5_snRNA:ITAG_infernal 13_snRNA:ITAG_infernal 10_snRNA:ITAG_infernal snRNA:ITAG_infernal 11_snRNA:ITAG_infernal 6_snRNA:ITAG_infernal snoRNA:ITAG_infernal 12_snRNA:ITAG_infernal U1_snRNA:ITAG_infernal group_I_intron:ITAG_infernal 7_snRNA:ITAG_infernal UTR_region:ITAG_infernal group_II_intron:ITAG_infernal miRNA:ITAG_infernal SRP_RNA:ITAG_infernal snRNA:ITAG_infernal snoRNA:ITAG_infernal 14_snRNA:ITAG_infernal 4_snRNA:ITAG_infernal miRNA:ITAG_infernal U2_snRNA:ITAG_infernal 9_snRNA:ITAG_infernal rRNA_5.8S:ITAG_infernal 8_snRNA:ITAG_infernal RNAi_reagent:ITAG_infernal
key          = Infernal
category     = Prediction Features
glyph        = segments
stranded     = 1
citation     = This track shows RNA regions inferred by <a href="http://infernal.janelia.org/">Infernal</a>.  Provided by Daniel Buchan, <a href="http://www.srcuk.org/">SRCUK</a>.

[trnascanse]
feature      = tRNA:ITAG_trnascanse
key          = tRNAscanSE
category     = Prediction Features
glyph        = segments
stranded     = 1
citation     = This track shows transfer RNAs predicted by the tRNAscan-SE program.  For more on tRNAscan-SE, see <a href="http://selab.janelia.org/software#trnascanse">http://selab.janelia.org/software#trnascanse</a>.  (see also <a href="http://www.ab.wur.nl/TomatoWiki/AnTRNAScanSE000">ITAG Wiki</a>)

[sgn_loci]
feature      = match:ITAG_sgn_loci
key          = SGN Locus Sequences
category     = Genetic Loci
glyph        = segments
link         = /search/quick_search.pl?term=$name
stranded     = 1
citation     = This track shows regions of similarity of Eugene-predicted cDNA sequences with known sequences associated with SGN genetic loci, as detected by BLAT (see also <a href="http://www.ab.wur.nl/TomatoWiki/AnSGNLoci000">ITAG Wiki</a>)

#[markers:overview]
#feature      = match:ITAG_sgn_markers
#fgcolor      = black
#bgcolor      = yellow
#key          = Genetic Markers
#glyph        = generic

[sgn_markers]
feature      = match:ITAG_sgn_markers
key          = SGN Marker Sequences
category     = Genetic Loci
stranded     = 1
fgcolor      = black
link         = http://solgenomics.net/search/quick_search.pl?term=$name
bgcolor      = yellow
stranded     = 1
citation     = This track shows regions of similarity of the genomic sequence with known sequences associated with SGN genetic markers, as detected by GenomeThreader (see also <a href="http://www.ab.wur.nl/TomatoWiki/AnSgnMarkers000">ITAG Wiki</a>)

[tilingpath]
feature      = clone:ITAG
glyph        = anchored_arrow
linewidth    = 2
stranded     = 0
category     = Genome data and reagents
das category = structural
fgcolor      = black
bgcolor      = black
key          = Tiling BACs and Fosmids
link         = /search/quick_search.pl?term=$name
citation     = Shows the genome tiling path of BACs and Fosmids, provided by each sequencing center.  For an AGP view of tiling paths, see <a href="http://solgenomics.net/sequencing/agp.pl">Tomato Assembly Display</a> on SGN.
|
  ;
  close $f;
}

sub write_gbrowse_prot {
  my ($release, $verbose) = @_;

  my $type = 'protein';

  my $file = $release->get_file_info(lc "gbrowse_${type}_conf")->{file};

  print "writing $type conf to ' $file '\n" if $verbose;

  my $releasenum = $release->release_number;
  my $is_devel = $release->is_devel_release ? 1 : 0;
  my $is_pre = $release->is_pre_release ? 1 : 0;

  open my $f, '>', $file or die "$! writing $file";

  print $f gbrowse_conf_head($release,$type);

  print $f q|
default features =

aggregators =

# "automatic" classes to try when an unqualified identifier is given
automatic classes = Sequence clone match nucleotide_motif

### TRACK CONFIGURATION ####
# the remainder of the sections configure individual tracks

#track to show DNA and GC content
[peptides]
glyph          = protein
global feature = 1
height         = 40
do_gc          = 1
fgcolor        = red
category       = General
axis_color     = blue
strand         = both
key            = Peptides/KD Plot
citation       = This track displays a Kyte-Doolittle hydropathy plot of the protein sequence at low magnifications and the peptide sequence itself at higher magnifications.


[blastp_trembl]
feature      = match:ITAG_blastp_trembl
key          = TrEMBL
category     = Similarity
glyph        = segments
stranded     = 1
citation     = This track shows similarities with sequences in TrEMBL (see <a href="http://www.ebi.ac.uk/trembl/">http://www.ebi.ac.uk/trembl/</a>), as detected by BLAST (see also <a href="http://www.ab.wur.nl/TomatoWiki/AnBlastP000">ITAG Wiki</a>)
bgcolor      = lightgreen

[blastp_refseq_pep]
feature      = match:ITAG_blastp_refseq_pep
key          = RefSeq Peptides
category     = Similarity
bgcolor      = lightgreen
glyph        = segments
stranded     = 1
citation     = This track shows similarities with peptide sequences in RefSeq (see <a href="http://www.ncbi.nlm.nih.gov/RefSeq/">http://www.ncbi.nlm.nih.gov/RefSeq/</a>), as detected by BLAST (see also <a href="http://www.ab.wur.nl/TomatoWiki/AnBlastP000">ITAG Wiki</a>)

[blastp_rice_pep]
feature      = match:ITAG_blastp_rice_pep
key          = Rice Peptides
category     = Similarity
bgcolor      = cyan
glyph        = segments
stranded     = 1
citation     = This track shows similarities with Rice peptide sequences from RAP-DB (see <a href="http://rapdb.dna.affrc.go.jp/rapdownload/">http://rapdb.dna.affrc.go.jp/rapdownload/</a>), as detected by BLAST (see also <a href="http://www.ab.wur.nl/TomatoWiki/AnBlastP000">ITAG Wiki</a>)

[blastp_swissprot]
feature      = match:ITAG_blastp_swissprot
key          = SwissProt
category     = Similarity
bgcolor      = lightgreen
glyph        = segments
stranded     = 1
citation     = This track shows similarities with peptide sequences in SwissProt (see <a href="http://www.expasy.ch/sprot/">http://www.expasy.ch/sprot/</a>), as detected by BLAST (see also <a href="http://www.ab.wur.nl/TomatoWiki/AnBlastP000">ITAG Wiki</a>)

[blastp_ath_pep]
feature      = match:ITAG_blastp_ath_pep
key          = Arabidopsis Peptides
category     = Similarity
glyph        = segments
bgcolor      = cyan
stranded     = 1
citation     = This track shows similarities with peptide sequences from <i>Arabidopsis thaliana</i> (see <a href="ftp://ftp.arabidopsis.org/home/tair/Sequences/blast_datasets/TAIR8_blastsets">ftp://ftp.arabidopsis.org/home/tair/Sequences/blast_datasets/TAIR8_blastsets</a>), as detected by BLAST (see also <a href="http://www.ab.wur.nl/TomatoWiki/AnBlastP000">ITAG Wiki</a>)

[interpro]
feature      = polypeptide_region:ITAG_interpro
key          = Interpro
bgcolor      = blue
category     = Protein Domain Similarities
glyph        = segments
stranded     = 1
citation     = Protein domains found by Interproscan (see also <a href="http://www.ab.wur.nl/TomatoWiki/AnInterpro000">ITAG Wiki</a>)

[rpsblast]
feature      = match:ITAG_rpsblast
key          = RPS-BLAST
bgcolor      = lightblue
category     = Protein Domain Similarities
glyph        = segments
stranded     = 1
citation     = Protein domain matches found by RPS-BLAST (see <a href="http://www.ab.wur.nl/TomatoWiki/AnRPSBlast000">ITAG Wiki</a>)

|;
  close $f;
}


sub restrict_batch_numbers {
  my ($bspec,@batchnums) = @_;
  my $max_batch = 0+(sort {$a <=> $b} @batchnums)[-1];
  my $san_bspec = $bspec;
  $san_bspec =~ s/[^\d\-,]//g; #< sanitize the batches spec so that people can't pass in arbitrary code
  $san_bspec =~ s/-/../g; #< convert ranges to perl array slice syntax

  my @expanded_bspec = eval "(0..$max_batch)[$bspec]";
  die "invalid batch range '$bspec'" if $EVAL_ERROR;
  my %bn_lookup = map {$_+0 => 1} @batchnums;
  my @new_batches = sort {$a <=> $b} grep $bn_lookup{$_}, @expanded_bspec;
  print "restricting to batch list: ".join(',',@new_batches)."\n";
  return @new_batches;
}


# this function is self-explanatory
sub copy_or_die {
  my ($file,$fh) = @_;
  open my $in, $file or die "$! reading $file";
  eval {
    while(my $line = <$in> ){
        $line =~ s/[[:cntrl:]]//g; #< remove any control chars
        $fh->print($line);
    }
  }; if($EVAL_ERROR) {
    confess "($file,$fh): $EVAL_ERROR";
  }
  close $in;
}

#same as above, except do a little munging on the gff3
sub copy_gff3_or_die {
  my ($aname,$sub,$file,@fh) = @_;
  open my $in, $file or die "$! reading $file";
  #sub check { my $l = shift; return unless $l =~ /\S/; chomp $l;  my @f = split /\s+/,$l,9; @f == 9 or $l =~ /^#/ or confess "'$l' not valid gff3";};
  if( $sub ) {
    while(<$in>) {
      #change the source column to be ITAG_ plus the analysis name
      my $result = $sub->($_);
      $result =~ s/[[:cntrl:]]//g; #< remove any control chars
      $result =~ s/^\s*(\S+)\t\S+/$1\tITAG_$aname/;
      #check($result);
      for my $fh (@fh) {
	$fh->print( $result );
      }
    }
  } else {
    while(my $line = <$in>) {
      #change the source column to be ITAG_ plus the analysis name
      $line =~ s/^\s*(\S+)\t\S+/$1\tITAG_$aname/;
      $line =~ s/[[:cntrl:]]//g; #< remove any control chars
      #check($line);
      for my $fh (@fh) {
	$fh->print( $line );
      }
    }
  }
  close $in;
}

# couple of functions for aggregating and reporting concisely what files are not available
sub warn_not_available_summary {
  our %not_available;
  my @alist = grep !/^_/, keys %not_available
    or return;
  my $alist = do {
    my $l = pop @alist;
    @alist ? 'analyses '.join(', ', @alist, "and $l") : "analysis $l";
  };
  warn "WARNING: ".$not_available{_global}." files not available from $alist.  Results will be incomplete\n";
}
sub report_file_not_available {
  my ($seq_rec,$a_rec,$file) = @_;
  our %not_available;
  push @{ $not_available{ $a_rec->{name} }{ $seq_rec->{batch} } },$file;
  $not_available{_global}++;
}

# parses a tabular file of GO terms in the format:
# scaffold00010_11.1.1    0003824,0008643,0000166
# returns a hashref like:
# { 'scaffold00010_11.1.1' => ['GO:0003824', 'GO:0008643', 'GO:0000166'],
#   ...
# }
sub get_go_terms_for_mrnas {
    my ( $go_tabular_file ) = @_;

    my %mrna_terms;

    open my $t, '<', $go_tabular_file or die "$! reading '$go_tabular_file'";
    while ( my $line = <$t> ) {
        # parse the line
        chomp $line;

        my ( $mrna_name, $go_nums ) = split /\s+/, $line, 2;
	next unless $go_nums =~ /\d{4,}/;
        my @go_terms = map "GO:$_", split /[\s,]+/, $go_nums;

        # check the go nums
        /^GO:\d+$/ or die "invalid go term $_ in list '$go_nums'" for @go_terms;

        # record them in the hash
        $mrna_terms{$mrna_name} = \@go_terms;
    }

    return \%mrna_terms;
}

my $dbh;
sub get_dbh {
    require CXGN::DB::Connection;
    $dbh ||= CXGN::DB::Connection->new;
}
