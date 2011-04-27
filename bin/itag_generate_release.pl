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

    -S use saved statistics file if available

    -R just collect stats and write a new README for the release

    -f force overwrite of existing release dir, will overwrite all files

=head1 MAINTAINER

Robert Buels

=head1 AUTHOR(S)

Robert Buels, E<lt>rmb32@cornell.eduE<gt>

=head1 COPYRIGHT & LICENSE

Copyright 2009-2011 The Boyce Thompson Institute for Plant Research

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
use Storable qw/ nstore retrieve dclone /;

use IO::File;
use IO::Pipe;
use Hash::Util qw/ lock_hash lock_keys /;
use File::Copy;
use File::Basename;
use File::Temp qw/tempfile/;
use List::Util qw/ min max /;
use Statistics::Descriptive;
use Tie::Function;
use Try::Tiny;
use URI::Escape;

use Data::Dumper;

use CXGN::TomatoGenome::BACPublish qw/ genbank_acc_to_seq_name seq_name_to_genbank_acc /;
use CXGN::ITAG::Pipeline;
use CXGN::ITAG::Pipeline::Analysis;
use CXGN::ITAG::Release;
use CXGN::ITAG::Release::Statistics;
use CXGN::ITAG::Release::View::README;

use CXGN::BioTools::AGP qw/agp_parse/;

use CXGN::Tools::Identifiers qw/ parse_identifier identifier_namespace /;
use CXGN::Tools::List qw/str_in flatten/;
use CXGN::Tools::Text qw/commify_number/;

### parse and validate command line args
our %opt;
getopts('d:p:b:PDSRf',\%opt) or pod2usage(1);
$opt{m} ||= 1_000_000;
$opt{P} && $opt{D} and die "-P and -D are mutually exclusive\n";
$opt{d} or pod2usage('must provide -d option');

my $release_num = $opt{D} ? 0 : shift @ARGV;
my $target_path = shift @ARGV or pod2usage();
($release_num + 0 == $release_num && $release_num >= 0) or die "release_num argument must be numeric\n";
-d $target_path or die "target path $target_path does not exist\n";
-w $target_path or die "target path $target_path is not writable\n"
    unless $opt{S} || $opt{R};

# set the root dir to our target directory so the object will make it in the right place#
CXGN::ITAG::Release->releases_root_dir($target_path);

# make our new release object
my $release = CXGN::ITAG::Release->new( releasenum => $release_num,
					pre   => $opt{P},
					devel => $opt{D},
				      );


# if -R, collect stats and write a new readme
if( $opt{R} || $opt{S} ) {
    #now collect statistics about this release to use in the README file
    my $saved_stat_file = File::Spec->catfile( $release->dir, 'statistics.dat' );
    my $stats = $opt{S} && -f $saved_stat_file
        ? retrieve( $saved_stat_file )    :
          $release->calculate_statistics;
    nstore( $stats, $saved_stat_file );

    write_readme( $release, $stats );
    exit;
}

###### MAKE OUTPUT DIRECTORY AND OPEN OUTPUT FILES

dump_data( $release );

#now collect statistics about this release to use in the README file
my $stats = $release->calculate_statistics;

#write the readme file
print "writing readme ...\n";
write_readme( $release, $stats );

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
        }
	elsif( $opt{f} ) {
            system "rm -rf $d/*";
            die if $CHILD_ERROR;
	}
	else {
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
	my $ns = identifier_namespace( $ctg_name );
        if ( $ns && $ns eq 'tomato_bac_contig' ) {
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
                      my $line = shift;
                      $line = add_name_attr( $line );
                      $line = add_functional_annotations( $line, $go_terms, $human_readable_descriptions );
                      $line =~ s/\tITAG_renaming\t/\tITAG_eugene\t/; #< change the source back to the original
                      return $line;
                  },
                  errors_fatal => 1,
                },
                { analyses => 'augustus',
                  gff3_output_spec_index => 0,
                  alter_lines_with => [ \&add_name_attr, \&filter_augustus ],
                  release_files => [qw| combi_genomic_gff3 genefinders_gff3 |],
                },
                { analyses => [qw[ geneid_tomato
                                   genemark_ath
                                   genemark_tom
                                   glimmerhmm_ath
                                   glimmerhmm_tomato
                                 ]
                              ],
                  gff3_output_spec_index => 0,
                  alter_lines_with => \&add_name_attr,
                  release_files => [qw| combi_genomic_gff3 genefinders_gff3 |],
                },
                { analyses => 'trnascanse',
                  gff3_output_spec_index => 1,
                  release_files => [qw| combi_genomic_gff3 genefinders_gff3 |],
                },
                { analyses => [qw[ tobacco_contigs potato_bacs itag1_ref ]],
                  gff3_output_spec_index => 1,
                  release_files => [qw| combi_genomic_gff3 other_genomes_gff3 |],
                },
                { analyses => 'tomato_bacs',
                  gff3_output_spec_index => 1,
                  release_files => [qw| combi_genomic_gff3 genomic_reagents_gff3 |],
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
                { analyses => 'renaming',
                  gff3_output_spec_index => 4,
                  release_files => ['combi_genomic_gff3'],
                  alter_lines_with => sub {
                    my ($line) = @_;
                    # change the source back to the original
                    $line =~ s/\tITAG_renaming\t/\tITAG_infernal\t/;
                    $line = generalize_infernal_types( $line );
                    return $line;
                  },
                },
                { analyses => 'microtom_flcdnas',
                  gff3_output_spec_index => 1,
                  release_files => [qw| combi_genomic_gff3 cdna_algn_gff3 |],
                },
                { analyses => qr/^transcripts_/i,
                  gff3_output_spec_index => 0,
                  filter_fh => \&transcript_aggregating_pipe,
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

                my $subs  = $dump->{alter_lines_with};
                $subs = [$subs] if $subs && !( ref($subs) eq 'ARRAY' );

                my $file = ( $arecord->{obj}->files_for_seq($ctg->{batch},$ctg->{name}) )[$index];
                #print STDERR "dumping $arecord->{name}:$file to ".join(' ',flatten $dump->{release_files})."\n";
                unless ( $dump->{errors_fatal} || -f $file ) {
                    report_file_not_available( $ctg, $arecord, $file );
                } else {
                    open my $in_fh, $file or die "$! reading $file";

                    # apply any filehandle filters to the gff3 input fh
                    if( my $filters = $dump->{filter_fh} ) {
                        $filters = [$filters] unless ref $filters eq 'ARRAY';
                        $in_fh = $_->( $in_fh ) for @$filters;
                    }

                    try {
                        copy_gff3_or_die( $arecord->{name}, $subs, $in_fh, @filehandles);
                    } catch {
                        use Data::Dump;
                        warn Data::Dump::dump( $subs );
                        warn Data::Dump::dump( $dump );
                        die $_;
                    };
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
  my @optional_keys = qw| alter_lines_with errors_fatal filter_fh |;
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
  $_[0] =~ s/name=(\S+)\s+/'Name='.gff3_escape($1).';Note='/e;
  $_[0] =~ s/Note=([^;]+)/'Note='.gff3_escape($1)/e;
  $_[0] =~ s/Name="/Name=/;
  $_[0] =~ s/Note=([^;]+)";/Note=$1;/;

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
      $gff3_line .= ";Name=$id" unless $gff3_line =~ /Name=/; # line must have a display name in order to be indexed
      $gff3_line .= ";Note=$desc_string"; #< descriptions are already gff3-escaped
  }

  # add ontology terms if present
  if( my $terms = $ontology_terms->{$id} ) {
      ### got terms: $terms
      $gff3_line .= ";Ontology_term=$_" for @$terms;
  }

  $gff3_line .= "\n";

  return $gff3_line;
}

sub add_name_attr {
    my $line = shift;
    return $line if $line =~ /Name=/;
    return $line unless $line =~ /ID=([^;\n]+)/;

    my $name = $1;
    $name =~ s/^[a-z]://i;
    chomp $line;
    $line .= ";Name=$name\n";

    return $line;
}

# skips gff3 lines of type transcription_*_site or intron
sub filter_augustus {
    my $line = shift;
    return '' if/ \t ( transcription_(start|end)_site | intron ) \t /x;
    return $line;
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

#given a gff3 file, sort it and clean up its pragmas
sub postprocess_gff3 {
  my ($seqs_file,$gff3_file) = @_;

  print "postprocessing gff3 $gff3_file with reference sequences $seqs_file\n";

  -f $seqs_file or die "seqs file '$seqs_file' not found!";
  -f $gff3_file or die "gff3 file '$gff3_file' not found!";

  #open a temp file for our output
  my ($temp_out_fh,$temp_out_file) = tempfile(File::Spec->catfile(File::Spec->tmpdir,"$FindBin::Script-gff3-postproc-XXXXXX"), UNLINK => 1);

  #sort the gff3 file by ref sequence and start coordinate into the temp file,
  #and remove duplicates
  open my $sg, "gff3_reformat.pl -l 5000 -s -i -S $seqs_file -U -u '-i' $gff3_file |"
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

      print $temp_out_fh join "\t",@fields;
      print $temp_out_fh "\n";
    }
  }
  close $temp_out_fh;

  #now open the original file, print the pragmas we found, and copy the rest into it
  open my $orig, ">$gff3_file" or die "$! writing $gff3_file";
  print $orig "##gff-version 3\n";
  print $orig "##$_\n" foreach @other_pragmas;
  $orig->flush();

  unless( -s $temp_out_file >= (-s $gff3_file)*0.95 ) {
      die "something wrong in postproc step, $temp_out_file is less than 95% of the size of $gff3_file";
  }
  copy_or_die( $temp_out_file, $orig );
  close $orig;
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
        $fh->print($line);
    }
  }; if($EVAL_ERROR) {
    confess "($file,$fh): $EVAL_ERROR";
  }
  close $in;
}

#same as above, except do a little munging on the gff3
sub copy_gff3_or_die {
  my ( $aname, $subs, $in_fh, @out_fhs ) = @_;
  $subs ||= [];
  #sub check { my $l = shift; return unless $l =~ /\S/; chomp $l;  my @f = split /\s+/,$l,9; @f == 9 or $l =~ /^#/ or confess "'$l' not valid gff3";};
  while( my $result = <$in_fh> ) {
      #change the source column to be ITAG_ plus the analysis name
      $result =~ s/^\s*(\S+)\t\S+/$1\tITAG_$aname/;
      $result = $_->($result) for @$subs;
      $result .= "\n" unless $result =~ /\n$/;
      #check($result);
      for my $fh (@out_fhs) {
          $fh->print( $result );
      }
  }
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



sub generalize_infernal_types {
    my ($line) = @_;
    # for easier config and loading, change the type of infernal lines
    # to 'transcript' and put the real type in rna_type
    if( $line =~ s/ (?<= INFERNAL \t ) ([^\t]+) (?= \t )/transcript/xi ) {
	my $type = $1;
        $line =~ s/[;\s]+$//;
        $line .= ";rna_type=$type\n";
    }
    return $line;
}


sub transcript_aggregating_pipe {
    my $orig_fh = shift;

    my $pipe = IO::Pipe->new;
    if( fork ) {
        # parent process will read the aggregated transcript gff3 from
        # the pipe

        $pipe->reader;
        close $orig_fh;
        return $pipe;

    } else {
        # forked process will aggregate the transcript gff3 and print
        # into the pipe, then hard-exit

        $pipe->writer;
        _aggregate_transcript_gff3( $orig_fh, $pipe );
        $pipe->flush;
        POSIX::_exit(0);

    }
}

sub _aggregate_transcript_gff3 {
    my ( $in_fh, $out_fh ) = @_;

    # hash the lines by Target
    my %seen;
    my %lines;
    while ( <$in_fh> ) {
        unless( /^(\S+)\t.+Target=(\S+)/ ) {
            next;
        }
        next if $seen{$_}++;

        my @f = split /\t/, $_, 9;
        push @{ $lines{"$1:$2"} }, \@f;
    }

    my @groups = map {
        my ($pts,$pte);
        my $pfe;
        my $curr_group = [];
        my @groups = ($curr_group);
        for my $l ( @$_ ) {
            my ($fs,$fe) = @{$l}[3,4];
            my $strand = $l->[6];
            my ($ts,$te) = $l->[8] =~ /Target=\S+ (\d+) (\d+)/ or die;

            if (    $pte
                        && $pfe
                            && (
                                ( $strand eq '+'
                                      ? abs( $pte - $ts ) > 3
                                          : abs( $pts - $te ) > 3
                                         )
                                    || $fs - $pfe > 30_000
                                   )
                           ) {
                $curr_group = [];
                push @groups, $curr_group;
            }
            push @$curr_group, $l;

            $pts = $ts;
            $pte = $te;
            $pfe = $fe;
        }
        @groups;
    } values %lines;

    # make feature groups out of them
    for my $g ( @groups ) {
        if ( @$g > 1 ) {
            my $superfeature = dclone($g->[0]);
            $superfeature->[5] = '.';
            my ( $sf_id ) = $superfeature->[8] =~ /ID=([^;]+)/
                or die "cannot parse line $superfeature";
            my @tcoords = map { $_->[8] =~ /Target=\S+ (\d+) (\d+)/ } @$g;
            @{$superfeature}[3,4] = ( min( map $_->[3], @$g ), max( map $_->[4], @$g ) );
            $superfeature->[8] =~ s/Target=(\S+) \d+ \d+/"Target=$1 ".min(@tcoords).' '.max(@tcoords)/e
                or die;

            for ( @$g ) {
                $_->[8] =~ s/ID=[^;]+;/Parent=$sf_id;/;
                $_->[2] = 'match_part';
            }
            $out_fh->print( join "\t", @$_ ) for $superfeature, @$g;
        } else {
            $out_fh->print( join "\t", @$_ ) for  @$g;
        }
    }
}

sub write_readme {
    my ( $release, $stats ) = @_;

    my $readme_text =
        CXGN::ITAG::Release::View::README
            ->new( release => $release, statistics => $stats)
            ->render;

    #and print it to the readme file
    my $gen_files = $release->get_all_files;
    open my $readme, '>',$gen_files->{readme}->{file}
    or die "$! opening $gen_files->{readme}->{file}";
    print $readme $readme_text;
    close $readme;

    return 1;
}
