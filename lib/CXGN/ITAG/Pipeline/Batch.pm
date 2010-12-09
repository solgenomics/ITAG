=head1 NAME

CXGN::ITAG::Pipeline::Batch - a batch in an ITAG pipeline

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 METHODS

=cut

package CXGN::ITAG::Pipeline::Batch;
use English;
use Carp;
use File::Spec;
use Memoize;

use UNIVERSAL qw/isa/;

use CXGN::Tools::Identifiers qw/identifier_namespace/;
use CXGN::TomatoGenome::ChromosomeAssemblies qw/named_contigs/;
use CXGN::TomatoGenome::BACPublish qw/seq_name_to_genbank_acc/;

use constant DEFAULT_BATCH_SIZE => 95;

sub _new {
  my ($class,$args) = @_;

  my $self = bless $args,$class;

  #set defaults
  $self->{pipeline} or croak "must provide a pipeline object";
  $self->{size} ||= DEFAULT_BATCH_SIZE;

  return $self;
}
sub _validate {
  my ($self) = @_;

  $self->{batchnum} >= 1
    or confess "invalid batchnum $self->{batchnum}";

  $self->{pipeline} && $self->{pipeline}->isa('CXGN::ITAG::Pipeline')
    or confess "invalid pipeline ($self->{pipeline})";

  !defined($self->{seqlist})
    or ref($self->{seqlist}) eq 'ARRAY'
      or confess "if set, seqlist must be an arrayref";

  all(map {isa($_,'Bio::SeqI')} @{$self->{seqlist}})
      or confess "seqlist must be composed of all Bio::SeqI objects";
}

=head2 open

  Desc  :
  Args  :
  Params:
  Side Effects:

=cut

sub open {
  my ($class,%args) = @_;
  my $self = $class->_new(\%args);

  $self->{batchnum} ||= $self->{pipeline}->_largest_existing_batchnum;

  return unless -d $self->dir;

  $self->seqlist
    or confess "batch $self->{batchnum} contains no sequences!";

  return $self;
}

=head2 create

  Status  :
  Usage   :
  Returns :
  Args    :
  Side Eff:
  Example :

  Briefly describe the function here.

=cut

sub create {
  my ($class,%args) = @_;
  my $self = $class->_new(\%args);
  #figure out the number of the new batch
  $self->{batchnum} ||= ($self->{pipeline}->_largest_existing_batchnum || 0) + 1;

  -d $self->dir
    and croak "batch $self->{batchnum} already exists";

  $args{size} && $args{size} >= 0
    or croak 'must give a size for the new batch of at least 1 seq';

  # validate seq_source argument
  require Scalar::Util;
  $args{seq_source}
      or croak 'must provide a seq source object or seq source class name';
  unless( Scalar::Util::blessed($args{seq_source}) ) {
      eval "require $args{seq_source}";
      $@ and die "could not load seq source $args{seq_source}:\n$@";
      $args{seq_source} = $args{seq_source}->new
  }

  $args{seq_source}->isa('CXGN::ITAG::SeqSourceI')
      or croak 'seq source must inherit from CXGN::ITAG::SeqSourceI';

  unless($self->{seqlist}) {
    #get list of sequences in all the other batches and index them.
    #also lists the member clones of tomato contig sequences
    my %taken_seqs;
    unless( $args{include_all} ) {
      %taken_seqs = map {
	my $other_batchnum = $_;
	
	map { $_=>1 } #< turn it into a hash
	map { #< add the member seqs of any tomato contigs
	  my $seqname = $_;
	  if( $args{no_bacs_in_contigs} && identifier_namespace($seqname) eq 'tomato_bac_contig' ) {
	    our %all_contigs;
	    %all_contigs = map {named_contigs($_, include_old => 1)} 1..12 unless %all_contigs;
	    my $rec = $all_contigs{$seqname} || die "contig '$seqname' not found in any AGP files!";
	    my @members =
	      map { seq_name_to_genbank_acc($_) || () } #< convert them into genbank accessions or
                                                        #  remove them if they have none
	      map {$_->{ident}}                         #< turn into list of idents
	      @$rec;                                    #< list of AGP lines from contig
	    #warn "adding ".join(',',$seqname,@members)."\n";
	    $seqname,@members
	  } else {
	    $seqname
	  }
	}
	#get a list of sequences in thase batch
	$class->open(pipeline => $self->{pipeline}, batchnum => $other_batchnum)->seqlist;
      } grep {
	$_ != $self->{batchnum}
      } $self->{pipeline}->list_batches;
    }

    #find sequences that aren't in the other batches
    my @new_seqlist;
    foreach my $seqname ($args{seq_source}->all_seq_names) {
      push @new_seqlist,$seqname unless $taken_seqs{$seqname};
      last if @new_seqlist >= $args{size};
    }

    unless(@new_seqlist) {
      warn "Could not create new batch $self->{batchnum}, no more sequences available.\n";
      return;
    }

    mkdir($self->dir)
      or croak "could not create batch dir ".$self->dir.": $!";

    #warn if the batch is smaller than requested
    unless($args{size} == @new_seqlist) {
      my $numseqs = @new_seqlist;
      my $and_not_part_of = $args{include_all} ? '' : ' and not part of any other batch';
      warn "WARNING: Requested new batch size of $args{size}, but only $numseqs were available$and_not_part_of.  Creating new batch of size $numseqs.\n";
      $self->{size} = scalar @new_seqlist;
    }

    #set this new seqlist as our seqlist and write it out
    $self->seqlist(@new_seqlist);
  }

  #initialize the filesystem space for each analysis that has a def
  #file
  foreach my $atag ($self->pipeline->list_analyses) {
    $self->pipeline->analysis($atag)->create_dirs($self);
  }

  return $self;
}

=head2 recreate

  Status  :
  Usage   :
  Returns :
  Args    :
  Side Eff:
  Example :

  Briefly describe the function here.

=cut

sub recreate {
  my ($self) = @_;

  foreach my $atag ($self->pipeline->list_analyses) {
    my $a = $self->pipeline->analysis($atag);
    next if $a->is_disabled;
    $a->create_dirs($self);
  }
}

=head2 delete

  Status  :
  Usage   :
  Returns :
  Args    :
  Side Eff:
  Example :

  Briefly describe the function here.

=cut

sub delete {
  my ($self) = @_;

  system 'mv', $self->dir, $self->pipeline->_deleted_subdir;
  $CHILD_ERROR
    and croak "could not mv ",$self->dir," -> ",$self->pipeline->_deleted_subdir;
}

=head2 batchnum

  Status  :
  Usage   :
  Returns :
  Args    :
  Side Eff:
  Example :

  Briefly describe the function here.

=cut

#read-only accessors
sub batchnum {
  shift->{batchnum}
}

=head2 size

  Status  :
  Usage   :
  Returns :
  Args    :
  Side Eff:
  Example :

  Briefly describe the function here.

=cut

sub size {
  shift->{size}
}

=head2 pipeline

  Status  :
  Usage   :
  Returns :
  Args    :
  Side Eff:
  Example :

  Briefly describe the function here.

=cut

sub pipeline {
  shift->{pipeline}
}

=head2 seqlist

  Status  :
  Usage   :
  Returns :
  Args    :
  Side Eff:
  Example :

  Briefly describe the function here.

=cut

#get/set the list of sequences in the batch
sub seqlist {
  my ($self,@list) = @_;

  my $seqlist_filename = File::Spec->catfile($self->dir,'seqlist.txt');

  #do we need to write out the seq list?
  if(@list) {
    ref and confess "seqlist should be an array of strings, found a '$_'\n" foreach @list;
    $self->{seqlist} = \@list;
    CORE::open my $seq_fh, ">$seqlist_filename"
      or die "Could not write to seqlist file $seqlist_filename: $!";
    foreach my $seqname (@list) {
      print $seq_fh "$seqname\n";
    }
  }

  #do we need to read in the seq list?
  if(! $self->{seqlist}  && -f $seqlist_filename) {
    CORE::open my $seqs_fh, $seqlist_filename
      or die "could not open seqlist file $seqlist_filename for reading: $!";
    my @seqlist;
    while(my $seqname = <$seqs_fh>) {
      $seqname =~ s/\s//g;
      push @seqlist,$seqname;
    }
    $self->{seqlist} = \@seqlist;
  }

#   use Data::Dumper;
#   warn 'seqlist returning '.Dumper($self->{seqlist});

  return @{$self->{seqlist}};
}

=head2 dir

  Status  :
  Usage   :
  Returns :
  Args    :
  Side Eff:
  Example :

  Briefly describe the function here.

=cut

#return the directory for this batch
memoize 'dir';
sub dir {
  my ($self) = @_;
  return File::Spec->catdir($self->{pipeline}->_pipedir,sprintf('batch%03d',$self->{batchnum}));
}

=head2 list_batches

  Status  :
  Usage   :
  Returns :
  Args    :
  Side Eff:
  Example :

  Briefly describe the function here.

=cut

#return a list of all the batch numbers present in a given pipeline,
#sorted in ascending numerical order
sub list_batches {
  my ($class,$pipeline,%args) = @_;
  my @batchlist = glob $pipeline->_pipedir.'/batch[0-9][0-9][0-9]/';
  if($args{include_deleted}) {
    push @batchlist, glob $pipeline->_deleted_subdir.'/batch[0-9][0-9][0-9]/';
  }
  return sort {$a <=> $b}  map {
    my ($num) = /batch(\d+)\/$/;
    $num
  } @batchlist;
}


=head1 MAINTAINER

Robert Buels

=head1 AUTHOR

Robert Buels, E<lt>rmb32@cornell.eduE<gt>

=head1 COPYRIGHT & LICENSE

Copyright 2009 Boyce Thompson Institute for Plant Research

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut


###
1;#do not remove
###
