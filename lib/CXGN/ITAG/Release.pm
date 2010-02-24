package CXGN::ITAG::Release;
use strict;
use warnings;
use English;
use Carp;

use Hash::Util qw/lock_keys lock_hash/;

use File::Spec;
use File::Basename;

use CXGN::DB::GFF::Versioned;
use CXGN::Tools::Class qw/parricide/;
use CXGN::Tools::Run;
use CXGN::ITAG::Config;
use CXGN::ITAG::Tools qw/ parse_release_dirname /;

=head1 NAME

CXGN::ITAG::Release - object representing an ITAG release, which
consists of a directory full of various data and documentation files,
plus a matching tar.gz file containing the whole release.

Mostly this object just provides a handle for a release directory,
providing services for making the proper names for files, finding
them, etc.

=head1 SYNOPSIS

  #couple of different ways of opening
  my $rel = CXGN::ITAG::Release->open( '/data/prod/ftpsite/tomato_genome/ITAG_devel_release' );
  my @matching_rels = CXGN::ITAG::Release->find( releasenum => 4 );

  #make a new one
  my $new_rel = CXGN::ITAG::Release->new( releasenum => 2 );

  #find the combined gff3 file
  my $combi_gff3 = $rel->get_file('combi_gff3');

=head1 DESCRIPTION

=head1 BASE CLASS(ES)

Class::Data::Inheritable (used for storing releases_root_dir and other class data)

=cut

use base qw/ Class::Data::Inheritable /;

=head1 SUBCLASSES

none

=head1 CLASS METHODS


=head2 new

  Usage: my $rel = CXGN::ITAG::Release->new( releasenum => 4 );
  Desc : make a new release object
  Args : releasenum => release number,
         devel => whether this is to be a development release (true or false),
         pre => whether this is to be a prerelease (true or false),
         dir => force directory to use for this release (defaults to <releaseses_root_dir>/<release_name>),
         basename => force dir basename to use (defaults to <release_name>)
  Ret  : a new CXGN::ITAG::Release object you can use to work with a release directory
  Side Effects: none

=cut

sub new {
  my ($class,%args) = @_;

  #check the keys passed in
  my %allowed_keys = map {$_=>1} qw/releasenum devel pre dir basename/;
  $allowed_keys{$_} or croak "$_ not allowed as argument to new()" foreach keys %args;

  #check values of args
  defined $args{releasenum} || $args{devel}
    or croak "must specify either 'releasenum' or 'devel' for new()";
  $args{devel} && $args{releasenum}
      and croak "cannot specify 'releasenum' if 'devel' is true";

  #normalize boolean args
  $args{$_} = !!$args{$_} foreach qw/devel pre/;

  return bless \%args,$class;
}

=head2 find

  Usage: my @rels = CXGN::ITAG::Release->find( releasenum => 4 );
  Desc : find published releases matching the given arguments
  Args : hash-style list like new() above,
         plus accepts the additional current => boolean
         search condition
  Ret  : list of matching release objects, possibly empty,
         in descending order of modification time
  Side Effects: looks in filesystem

=cut

sub find {
  my ($class,%args) = @_;

  my $dirname = $class->releases_root_dir;

  #subroutine to check if a given release matches the search conditions
  my $cond_match = sub {
    my ($rel) = @_;
    my %ac_map = qw/ pre is_pre_release devel is_devel_release releasenum release_number current is_current_release/;
        foreach my $k (keys %args) {
      my $v = $args{$k};
      my $ac_name = $ac_map{$k};
      $ac_name or croak "invalid condition '$k'";
      return 0 unless $rel->$ac_name eq $v;
    }
    return 1;
  };

  return
    sort { $b->{modtime} <=> $a->{modtime} } #< sort the results in descending order of modtime, and remember the modtimes
    map { $_->{modtime} = (stat($_->{dir}))[9]; $_ } #< decorate with modtimes
    grep $cond_match->($_),
    map $class->new(%$_),
    grep $_,
    map parse_release_dirname($_),
    glob(File::Spec->catfile($dirname,'*'));

}

=head2 releases_root_dir

  Usage: my $dir = CXGN::ITAG::Release->releases_root_dir
  Desc : get/set the directory where this object expects ITAG releases
         to be stored on the filesystem
  Args : (optional) directory to set
  Ret  : currently set directory
  Side Effects: dies if directory does not exist

=cut

# store the actual root dir here (see Class::Data::Inheritable)
__PACKAGE__->mk_classdata('_rootdir' =>
			  do { #< defaults to <conf:ftpsite_root>/<conf:itag_release_base>
                              my $cfg = CXGN::ITAG::Config->load;
                              my $ftpsite = $cfg->{'ftpsite_root'}
                                  or die "configuration variable 'ftpsite_root' not set";
                              my $itag_rel = $cfg->{'itag_release_base'}
                                  or die "configuration variable 'itag_release_base' not set";
                              File::Spec->catdir($ftpsite,$itag_rel);
			  }
			 );

sub releases_root_dir {
  my ($class,$newdir) = @_;
  $class->_rootdir($newdir) if $newdir;
  my $dir = $class->_rootdir;
  return $dir;
}


=head1 OBJECT METHODS

=head2 get_file_info

  Usage: my $f = $r->get_file_info('combi_gff3');
  Desc : get the full path to where a file *should* be located for
         this release.  It might not actually be there, you should
         check.
  Args : one of the file shortnames listed below
  Ret  : hashref of { file => file path string,
                      desc => text description of the file,
                      type => file type, one of
                              qw/gff3 tab fasta txt/,
                      seq_type => type of sequences contained or based on,
                                  one of qw/ genomic protein cdna cds /,
                    }
  Side Effects: none

=head3 File Shortnames

  genomic_fasta
  cdna_fasta
  cds_fasta
  protein_fasta
  contig_gff3
  models_gff3
  genefinders_gff3
  functional_prot_gff3
  cdna_algn_gff3
  combi_genomic_gff3
  contig_members_sol
  sol_gb_mapping
  readme
  ...

=cut

sub get_file_info {
  my ($self,$shortname) = @_;

  my $files = $self->get_all_files;

  $self->{known_files} ||= {map {$_=>1} keys %$files};
  croak "unknown release file shortname '$shortname'" unless $self->{known_files}->{$shortname};

  return $files->{$shortname};
}

=head2 get_all_files

  Usage: my $all = $rel->get_all_files
  Desc : get hashref of info for all files
  Args : none
  Ret  : hashref as
         {  shortname => file info hashref,
            ...
         }
  Side Effects: none

=cut

sub get_all_files {
  my ($self) = @_;

  $self->{files} ||=
    {
     genomic_fasta        => { file => $self->_assemble_release_filename("genomic_contigs","fasta"),
			       desc => 'fasta-format sequence file of genomic contig sequences.',
			       type => 'fasta',
			       seq_type => 'genomic',
			     },
     cdna_fasta           => { file => $self->_assemble_release_filename("cdna","fasta"),
			       desc => 'fasta-format sequence file of cDNA sequences.',
			       type => 'fasta',
			       seq_type => 'cdna',
			     },
     cds_fasta            => { file => $self->_assemble_release_filename("cds","fasta"),
			       desc => 'fasta-format sequence file of CDS sequences.',
			       type => 'fasta',
			       seq_type => 'cds',
			     },
     protein_fasta        => { file => $self->_assemble_release_filename("proteins","fasta"),
			       desc => 'fasta-format sequence file of protein sequences.',
			       type => 'fasta',
			       seq_type => 'protein',
			     },
     contig_gff3          => { file => $self->_assemble_release_filename("contig_assembly","gff3"),
			       desc => 'GFF version 3 file showing the assembly of clones in each contig.',
			       type => 'gff3',
			       seq_type => 'genomic',
			     },
     models_gff3          => { file => $self->_assemble_release_filename("gene_models","gff3"),
			       desc => 'GFF version 3 file containing gene models in this release.',
			       type => 'gff3',
			       seq_type => 'genomic',
			     },
     genefinders_gff3     => { file => $self->_assemble_release_filename("de_novo_gene_finders","gff3"),
			       desc => 'GFF version 3 file containing predictions from several de novo gene finders.  These were integrated into the final gene models by EuGene.',
			       type => 'gff3',
			       seq_type => 'genomic',
			     },
     sgn_data_gff3        => { file => $self->_assemble_release_filename("sgn_data","gff3"),
			       desc => 'GFF version 3 file containing alignments to sequences related to data on SGN.  Currently contains alignments to SGN unigenes, SGN marker sequences, and SGN locus sequences.',
			       type => 'gff3',
			       seq_type => 'genomic',
			     },
     functional_prot_gff3 => { file => $self->_assemble_release_filename("prot_functional","gff3" ),
			       desc => 'GFF version 3 file containing functional annotations to protein sequences.',
			       type => 'gff3',
			       seq_type => 'protein',
			     },
     cdna_algn_gff3      =>  { file => $self->_assemble_release_filename("cdna_alignments","gff3" ),
			       desc => 'GFF version 3 file containing alignments of existing EST and cDNA sequences to the genome.',
			       type => 'gff3',
			       seq_type => 'genomic',
			     },
     combi_genomic_gff3   => { file => $self->_assemble_release_filename("all_genomic","gff3"),
			       desc => 'GFF version 3 file containing all genomic annotations in this release.',
			       type => 'gff3',
			       seq_type => 'genomic',
			     },
     contig_members_sol   => { file => $self->_assemble_release_filename("contig_members_sol_ids","tab"),
			       desc => 'tab-delimited file listing which clone sequences comprise each contig (with clones listed by their SOL IDs).',
			       type => 'tab',
			       seq_type => '',
			     },
     sol_gb_mapping       => { file => $self->_assemble_release_filename("sol_gb_id_map","tab"),
			       desc => 'tab-delimited file giving the mapping between SOL-style sequence IDs (e.g. C01HBa0001A01) and GenBank accessions (e.g. AC123456).  Some BAC sequences may not yet have GenBank accessions registered with SGN, these are denoted by NO_ACCESSION_REGISTERED.',
			       type => 'tab',
			       seq_type => '',
			     },
     readme               => { file => $self->_assemble_release_filename("README","txt"),
			       desc => 'Overview of the release, with some statistics.',
			       type => 'txt',
			       seq_type => '',
			     },
     gbrowse_genomic_conf   => { file => $self->_assemble_release_filename("genomic","conf"),
				 desc => 'GBrowse configuration file used at SGN to display annotations on genomic sequences from this release',
				 type => 'conf',
				 seq_type => 'genomic',
			     },
     gbrowse_protein_conf   => { file => $self->_assemble_release_filename("protein","conf"),
				 desc => 'GBrowse configuration file used at SGN to display annotations on protein sequences from this release',
				 type => 'conf',
				 seq_type => 'protein',
			       },
#      readme_file_desc     => { file => assemble_release_filename("README_file_descriptions","txt"),
# 			       desc => 'this file',
# 			     },

    };
  lock_hash(%$_) foreach values %{$self->{files}};
  lock_hash(%{$self->{files}});
  return $self->{files};
}


sub _assemble_release_filename {
  my ($self,$name,$type) = @_;
  $name or croak "must pass name to assemble_release_filename()\n";
  $type or croak "must pass type to assemble_release_filename()\n";

  return File::Spec->catfile($self->dir,$self->release_tag.'_'.$name.'.'.$type);
}


=head2 release_number

  Usage: my $num = $rel->release_number
  Desc : get the release number of this release
  Args : none
  Ret  : integer release number, or 0 if no
         number (e.g. for a devel release)
  Side Effects: none

=cut

sub release_number {
  my ($self) = @_;
  return $self->{releasenum};
}


=head2 release_tag

  Usage: my $tag = $rel->release_tag
  Desc : get the 'release tag' for this release, which is a short
         string identifier for the release.  Example: ITAG1 would be
         the first official ITAG release
  Args : none
  Ret  : release tag string
  Side Effects: none

=cut

sub release_tag {
  my ($self) = @_;
  my $pre   = $self->is_pre_release ? '_pre'    : '';
  if( $self->is_devel_release ) {
    return "ITAG_devel$pre";
  } else {
    return "ITAG$self->{releasenum}${pre}";
  }
}

=head2 dir

  Usage: -d $rel->dir or die "no dir found for release!"
  Desc : get the full path for this release's directory,
         which may not exist if this object is new
  Args : (optional) true value if should return just the
         basename
  Ret  : string pathname, or basename if flag passed
  Side Effects: none

=cut

sub dir {
  my ($self,$arg) = @_;
  $arg && confess 'dir() takes no arguments';
  return $self->{dir} ||= File::Spec->catfile( $self->releases_root_dir, $self->dir_basename);
}

=head2 dir_basename

  Usage: my $bn = $rel->dir_basename
  Desc : get the basename of this analysis's dir, e.g. 'ITAG2_release'
  Args : none
  Ret  : string basename
  Side Effects: none

  Right now, this is actually just the release_name()

=cut

sub dir_basename {
  my ($self) = @_;
  return $self->{basename} || $self->release_name;
}

=head2 mkdir

  Usage: $rel->mkdir() or die "could not create release directory!";
  Desc : attempt to make the directory returned by dir()
  Args : none
  Ret  : true  if the creation was successful or the dir already existed,
         false if the creation was not successful, or the dir exists
               but is not writable
  Side Effects: might make a directory on the filesystem

=cut

sub mkdir {
  my ($self) = @_;

  my $d = $self->dir;

  if( -d $d ) {
    return -w $d;
  }

  return mkdir($d);
}


=head2 dir_modtime

  Usage: my $time = $rel->dir_modtime
  Desc : get the modification time of the release dir,
         as returned by (stat($dir))[9]
  Args : none
  Ret  : integer modification time
  Side Effects: does a stat()

=cut

sub dir_modtime {
  my ($self) = @_;
  return (stat($self->dir))[9];
}

=head2 text_description

  Usage: my $desc = $rel->text_description()
  Desc : get a short (max 80 chars) description of this release,
         useful for listings and such
  Args : none
  Ret  : text string
  Side Effects: none

=cut

sub text_description {
  my ($self) = @_;
  if( $self->is_devel_release ) {
    return 'ITAG development snapshot';
  } else {
    my $pre = $self->is_pre_release ? 'pre-' : '';
    my $num = $self->release_number;
    return "ITAG $num ${pre}release";
  }
}


=head2 tarfile_name

  Usage: my $t = $rel->tarfile_name
  Desc : get the full path to this release's tarfile
         may not exist if this object is new
  Args : none
  Ret  : string pathname
  Side Effects: none

=cut

sub tarfile_name {
  my ($self) = @_;
  my ($bn,$rootdir) = fileparse($self->dir);
  $rootdir ||= '..';
  return File::Spec->catfile($rootdir,$self->release_name.'.tar.gz');
}

=head2 make_tarfile

  Usage: $rel->make_tarfile
  Desc : run tar and gzip to tar up this release into the correct file
  Args : none
  Ret  : true on success
  Side Effects: dies on failure

=cut

sub make_tarfile {
  my ($self) = @_;

  my $tar = CXGN::Tools::Run->run( 'tar',
				   '-c',
				   -f => '-',
				   -C => $self->releases_root_dir,
				   $self->dir_basename,
				 );
  my $gz = CXGN::Tools::Run->run( 'gzip',
				  '-cn',
				  '--rsyncable',
				  $tar->out_file,
				);

  my $cp = CXGN::Tools::Run->run( 'cp',
				  '-v',
				  $gz->out_file,
				  $self->tarfile_name,
				);

  return 1 if -f $self->tarfile_name;
  return 0;
}


# helper method to get the name of this release, currently something
# like 'ITAG1_release' or 'ITAG_devel_release'

=head2 release_name

  Usage: my $name = $rel->release_name
  Desc : like release_tag(), except includes '_release'
         at the end of the returned string
  Args : none
  Ret  : string release name
  Side Effects: none

=cut

sub release_name {
  shift->release_tag.'_release';
}


=head2 has_gbrowse_dbs

  Usage: die 'not available in gbrowse' unless $rel->has_gbrowse_dbs
  Desc : check whether this release has been loaded into gbrowse
  Args : none
  Ret  : true if a corresponding gbrowse DB exists, false otherwise
  Side Effects: might open a new database connection

=cut

our @dbtypes = qw/genomic protein/;
sub has_gbrowse_dbs {
  my ($self) = @_;

  foreach my $type (our @dbtypes) {
    #warn "checking for existence of $db_name_root $type:\n";
    return 0 unless $self->get_annotation_db($type);
    #warn "present.\n";
  }
  return 1;
}

=head2 get_annotation_db

  Usage: my $bdbgff = $rel->get_annotation_db('protein');
  Desc : get a Bio::DB::GFF or Bio::DB::SeqFeature::Store
         object for accessing annotations from this release.
         Can be either 'genomic' annotations (annotations on
         genomic sequences) or 'protein' annotations (on
         protein sequences).
  Args : either 'genomic' or 'protein'
  Ret  : a Bio::DB::GFF or Bio::DB::SeqFeature::Store object,
         or nothing if the annotation db is not available

=cut

sub get_annotation_db {
  my ($self,$type) = @_;
  grep $type eq $_,our @dbtypes
    or croak "type arg to get_annotation_db() must be either 'genomic' or 'protein'";
  my $bdb = CXGN::DB::GFF::Versioned->new( -db => $self->dbname_root($type) );
  $bdb = $bdb->bdb
    or return;
  return $bdb;
}

=head2 dbname_root

  Usage: my $r = $rel->dbname_root('protein')
  Desc : get the root name of the postgres database used for storing
         annotations from this release of the given type, e.g.
         'ITAG_devel_release.protein'
  Args : 'genomic' or 'protein'
  Ret  : string dbname root

=cut

sub dbname_root {
  my ($self,$type) = @_;
  grep $type eq $_,our @dbtypes
    or croak "type arg to get_annotation_db() must be either 'genomic' or 'protein'";
  return $self->release_name.'.'.$type;
}


=head2 list_current_gbrowse_dbs

  Usage: my @dbs = $rel->list_current_gbrowse_dbs
  Desc : get a list of the current versions of this
         release's gbrowse databases
  Args : none
  Ret  : list of strings, or nothing if not loaded
  Side Effects: might open a new database connection
  Example:
     foreach my $d ($rel->list_current_gbrowse_dbs) {
       print "need to dump '$d'\n";
     }

=cut

sub list_current_gbrowse_dbs {
  my ($self) = @_;

  return map {
    #warn "checking for existence of $db_name_root $type:\n";
    my $bdb = CXGN::DB::GFF::Versioned->new( -db => $self->dbname_root($_) );
    return unless $bdb->bdb;
    $bdb->dbname
  } our @dbtypes;
}

=head2 is_devel_release

  Usage: print 'devel' if $rel->is_devel_release
  Desc : get whether this is a development release
  Args : none
  Ret  : boolean showing whether this is a development release
  Side Effects: none

=cut

sub is_devel_release {
  my ($self) = @_;
  return 1 if $self->{devel};
  return 0;
}

=head2 is_pre_release

  Usage: print 'prerelease' if $rel->is_pre_release
  Desc : get whether this is a pre-release
  Args : none
  Ret  : boolean showing whether this is a pre-release
  Side Effects: none

=cut

sub is_pre_release {
  my ($self) = @_;
  return 1 if $self->{pre};
  return 0;
}

=head2 is_current_release

  Usage: print 'this is the current' if $rel->is_current_release
  Desc : get whether this release is the current one (linked to be an ITAG_current symlink)
  Args : none
  Ret  : boolean showing whether this is the current release
  Side Effects: none

=cut

sub is_current_release {
  my ( $self, ) = @_;
  my $curr_link = File::Spec->catdir( $self->dir, '..', 'ITAG_current' );
#  warn "$curr_link\n";
  return 0 unless -l $curr_link;

  my $p1 = File::Spec->catdir( $self->dir, '..', readlink($curr_link) );
  my $p2 = $self->dir;
  my $i1 = (stat( $p1 ))[1];
  my $i2 = (stat( $p2 ))[1];
#  warn "$p1 - $i1, $p2 - $i2\n";
  return 1 if $i1 == $i2;

  return 0;
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
