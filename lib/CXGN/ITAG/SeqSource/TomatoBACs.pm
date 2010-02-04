package CXGN::ITAG::SeqSource::TomatoBACs;
use strict;
use warnings;
use English;
use Carp;

use Bio::PrimarySeq;

use CXGN::DB::Connection;
use CXGN::Tools::Class qw/parricide/;

=head1 NAME

CXGN::ITAG::SeqSource::TomatoBACs - ITAG pipeline sequence
source that provides tomato BAC sequences, using the BACs' genbank
accessions as sequence identifiers.

=head1 SYNOPSIS

  my $src = CXGN::ITAG::SeqSource::TomatoBACs->new(
                unfinished => 1,
                phase      => 3,
            );

=head1 DESCRIPTION

=head1 BASE CLASS(ES)

L<CXGN::ITAG::SeqSourceI>

=cut

use base qw/CXGN::ITAG::SeqSourceI Class::Accessor/;

=head1 METHODS

Provides the methods specified in the SeqSourceI interface, using
tomato BAC sequences.

By default, provides only finished tomato BAC sequences.  However,
if you pass unfinished => 1, will also provide unfinished BAC
sequences.

=cut

#Class::Accessor provides the new() method
__PACKAGE__->mk_accessors(qw/unfinished phase/);

sub all_seq_names {
  my ($self) = @_;

  #look up all available BAC genbank accessions in the chado DB and return them
  my $and_finished = $self->unfinished ? '' : <<EOS;
  (select count(*) from featureprop fp where fp.feature_id = f.feature_id and fp.value = '1' and fp.type_id in(select cvterm_id from cvterm where name='finished_seq')) >= 1
and
  f.residues not like '\%X%'
EOS
  my $phase = $self->phase;
  my $and_phase = $phase ? <<EOS : '';
  (select count(*) from featureprop fp where fp.feature_id = f.feature_id and fp.value = '$phase' and fp.type_id in(select cvterm_id from cvterm where name='htgs_phase')) >= 1
EOS

  my $where = join ' AND ', grep $_, $and_finished, $and_phase;
  $where = "WHERE $where" if $where;

  my $q = <<EOQ;
   select accession
   from
   ( select distinct on(cf.clone_id) cf.clone_id,feat.feature_id, dbx.accession, feat.residues
     from db
     join dbxref dbx using(db_id)
     join feature_dbxref fd on fd.dbxref_id=dbx.dbxref_id
     join feature feat using(feature_id)
     join clone_feature cf using(feature_id)
     join cvterm cvt on feat.type_id=cvt.cvterm_id
     where
        cvt.name = 'BAC_clone'
     and
       db.name = 'DB:GenBank_Accession'
     order by cf.clone_id,feat.timelastmodified DESC
   ) as f
   $where
EOQ

  our $LAST_QUERY = $q;
  my $sth = $self->db_cxgn->prepare_cached($q);

  $sth->execute();
  my $seqnames = $sth->fetchall_arrayref;
  #use Data::Dumper;
  #die Dumper $seqnames;

  #now return only the most recent for each of the seqnames
  my @names = sort {
    my ($aa,$av) = split /\./,$a;
    my ($ba,$bv) = split /\./,$b;
    $aa cmp $ba
      || $av <=> $bv
  }
  map $_->[0],
  @$seqnames;

  my %bases = map { my $i = $_; $i =~ s/\.\d+$//; $i => $_ } @names;

  return values %bases;
}

__PACKAGE__->set_sql(get_bac_seq => <<EOQ);
   select dbx.accession,f.name,f.residues
   from db
   join dbxref dbx using(db_id)
   join feature_dbxref fd on fd.dbxref_id=dbx.dbxref_id
   join feature f using(feature_id)
   where
       db.name = 'DB:GenBank_Accession'
     and dbx.accession = ?
EOQ

sub get_seq {
  my ($self, $seqname) = @_;

  #look up that accession in the chado db, stuff it in a
  #Bio::PrimarySeq, and return it
  my $sth = $self->sql_get_bac_seq;
  $sth->execute($seqname);
  my $seqs = $sth->fetchrow_arrayref;

  my $newseq = Bio::PrimarySeq->new( -id   => $seqs->[0],
				     -desc => $seqs->[1],
				     -seq  => $seqs->[2],
				   );

  return ( $newseq,
	   [{ ostart => 1,
	      oend => $newseq->length,
	      objname => 'tomato-clone-'.($newseq->display_id),
	      partnum => 1,
	      linenum => 1,
	      type => 'F',
	      typedesc => 'finished',
	      ident => $newseq->display_id,
	      cstart => 1,
	      cend => $newseq->length,
	      orient => '+',
	    }],
	 );
}

=head1 AUTHOR(S)

Robert Buels

=cut

sub DESTROY {
  my $this = shift;
  our @ISA;
  return parricide($this,@ISA);
}

###
1;#do not remove
###
