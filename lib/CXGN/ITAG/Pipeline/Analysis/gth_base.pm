package CXGN::ITAG::Pipeline::Analysis::gth_base;
use strict;
use warnings;
use English;
use Carp;

use Fatal qw/ open close chdir /;
use File::Basename qw/basename/;
use File::Path qw/ mkpath rmtree /;

use Bio::FeatureIO;
use Bio::SeqIO;

=head1 NAME

CXGN::ITAG::Pipeline::Analysis::gth_base - pipeline analysis base
class for genomethreader-based analyses

=head1 BASE CLASS(ES)

L<CXGN::ITAG::Pipeline::Analysis>

=cut

use base qw/CXGN::ITAG::Pipeline::Analysis/;

=head1 METHODS

=head2 locally_runnable

returns 1

=cut

#this analysis is locally runnable
sub locally_runnable { 1 }

#and this is the routine to run it
sub run {
    my ($self,$batch) = @_;


    # for each seq file, set up and run an analysis job on the
    # cluster, format its output, and queue it up to be copied into
    # place
    my @jobs;
    my @ops;
    foreach my $seqname ($batch->seqlist) {
        my ($outfile,$gff3_out_file) = map $self->cluster_temp($seqname,$_), 'xml','gff3';

        my $seqfile = $self->_seq_file($batch,$seqname);

        my $gs_est_job = $self->cluster_run_class_method( $batch,
							  $seqname,
							  'run_gth',
                                                          $seqname,
                                                          $seqfile,
							  $outfile,
                                                          $gff3_out_file,
							  #{ err_file => $errfile,
							  #},
                                                         );

        my ($out_dest,$gff_out_dest) = $self->files_for_seq($batch,$seqname);
	push @jobs, $gs_est_job;
	push @ops, [$gff3_out_file => $gff_out_dest ], [ $outfile => $out_dest ];

	# check for any failed jobs, so that we don't have to submit
	# all the jobs before we die.  submitting all the jobs could
	# take a very long time if there are many jobs, because the
	# job submission will block and wait for the torque server to
	# clear a bit if there are too many jobs in its queue.
	$_->alive for @jobs;
    }

    # wait for all the jobs to finish (also running their
    # on_completion hooks)
    sleep 10 while grep $_->alive, @jobs;

    #atomically move the results into position
    $self->atomic_move(@ops);
}

# class method, meant to be run on a cluster node through perl -e
sub run_gth {
    my ( $class, $seqname, $seqfile, $outfile, $gff3_out_file ) = @_;

    my $work_dir = $class->local_temp($seqname);
    mkdir $work_dir;

    local $SIG{TERM} = sub {
	rmtree( $work_dir );
	exit;
    };

    my $un_xed_seqs = $class->local_temp( $seqname, "un_xed.seq" );
    $class->make_un_xed_seqfile( $seqfile => $un_xed_seqs );

    my $cdna_file = Bio::SeqIO->new( -format => 'fasta',
                                     -file   => $class->_cdna_file,
                                   );
    # make xml outfile just a stub
    { open my $outfile_fh, '>', $outfile
          or die "$! writing outfile $outfile";
      $outfile_fh->print( <<EOF );
<?xml version="1.0" encoding="ISO-8859-1"?>
<!-- GenomeThreader XML output no longer provided, only GFF3 -->
EOF
    }

    open my $gff3_out_fh, '>', $gff3_out_file
        or die "$! writing gff3 file $gff3_out_file";
    { my @seqregions = $class->seqregion_statements( $un_xed_seqs );
      $gff3_out_fh->print("##gff-version 3\n@seqregions");
    }

    uniq_flush();
    while( my $seq = $cdna_file->next_seq ) {
        $class->build_and_run_batch( $work_dir, $un_xed_seqs, $seq, $cdna_file, $gff3_out_fh );
    }

    rmtree( $work_dir );
}


sub build_and_run_batch {
    my ($class, $work_dir, $un_xed_seqs, $first_seq, $cdna_file, $gff3_out_fh ) = @_;

    # make a temp dir for this batch
    my $tempdir = File::Temp->newdir( DIR => $work_dir );

    my $temp_cdna = "$tempdir/cdna";
    my $batch_seqio = Bio::SeqIO->new( -format => 'fasta', -file => ">$temp_cdna" );
    $batch_seqio->write_seq( $first_seq );

    # write out the proper number of sequences for this batch, it's ok
    # if it's a little over
    my $batch_size_target = 10_000_000; #< run the sequences 10MB at a time
    my $batch_size = $first_seq->length;
    while ( $batch_size < $batch_size_target && ( my $seq = $cdna_file->next_seq ) ) {
        $batch_seqio->write_seq( $seq );
        $batch_size += $seq->length;
    }

    # run the batch
    my $temp_gff3 = "$tempdir/gff3";
    #warn "running batch in $tempdir";

    my @cmd =
        ( 'gth',
          '-gff3out',
          '-intermediate',
          '-force',
          -autointroncutout  => 400,
          '-gcmaxgapwidth'   => 40000, #< max gap in any alignment, this is 2x the largest intron size we have seen so far
          -minalignmentscore => '0.90',
          -mincoverage       => '0.90',
          -seedlength        => 16,
          -o                 => "$temp_gff3",
          -species           => 'arabidopsis',
          -cdna              => "$temp_cdna",
          -genomic           => $un_xed_seqs,
         );
    ### command: @cmd
    system @cmd
        and die "$! running @cmd";

    #now convert the gthxml to gff3
    #$class->_gthxml_to_gff3( "$temp_out", $un_xed_seqs, $seqname, "$temp_gff3" );

    #append the gff3 to the overall gff3 file
    $class->_reformat_gff3( $temp_gff3, $gff3_out_fh );
}

sub _reformat_gff3 {
    my ( $class, $temp_gff3, $gff3_out_fh ) = @_;

    my %type_map = (
        gene                    => 'match',
        exon                    => 'match_part',
        five_prime_splice_site  => 'five_prime_splice_site',
        three_prime_splice_site => 'three_prime_splice_site',
       );

    my %id_map;

    open my $g, $temp_gff3 or die 'cannot open temp gff3??';
    while( my $line = <$g> ) {
        if( $line =~ /^\s*(#+)/ ) {
            next unless length $1 == 3;
            gff3_uniq_print($gff3_out_fh, $line);
        } else {
        #parse line
            chomp $line;
            my @f = split /\t/,$line,9;
            my %attrs = map { split /=/, $_, 2 } split /;/, pop @f;
            if ( $attrs{Target} && !$attrs{Name} ) {
                my ($tn) = $attrs{Target} =~ /(\S+)/;
                $attrs{Name} = $tn;
            }

            # map source
            $f[1] = $class->_source;

            # map type
            die "unknown type $f[2]"
                unless exists $type_map{ $f[2] };
            next unless $type_map{ $f[2] };
            $f[2] = $type_map{ $f[2] };

            #rebuild line
            $line = join "\t", @f, join ';',map { "$_=$attrs{$_}" } sort keys %attrs;

            gff3_uniq_print( $gff3_out_fh, $line."\n" );
        }
    }
}


sub seqregion_statements {
    my ($class, $seqfile) = @_;
    my $i = Bio::SeqIO->new( -file => $seqfile, -format => 'fasta');
    my @regions;
    while( my $s = $i->next_seq ) {
        push @regions, join(' ', '##sequence-region', $s->id, 1, $s->length )."\n";
    }
    return @regions;
}

sub gff3_uniq_print {
    my ($filehandle, $line) = @_;

    if( $line =~ /^\s*(#+)/ ) {
        uniq_flush() if $line =~ /^\s*###\s+$/;
    } elsif( $line =~ /\S/ ) {
        $line =~ s/(?<=Parent=)([^;\s]+)/uniq_parent_id($1)/e;
        $line =~ s/(?<=ID=)([^;\s]+)/uniq_id($1)/e;
    }
    $filehandle->print( $line );
}
{
    my %uniq_ctr;
    my %uniq_mapping;
    sub uniq_id {
        my ($id) = @_;
        my $idx = $uniq_ctr{$id}++;
        my $new = $id.'_'.$idx;
        $uniq_mapping{$id} = $new;
        return $new;
    }
    sub uniq_parent_id {
        my ($id) = @_;
        return $uniq_mapping{$id} || die "parent element with ID='$id' not found!\n";
    }
    sub uniq_flush {
        %uniq_mapping = ();
    }
}


sub make_un_xed_seqfile {
  my ($class,$seqfile,$un_xed_seqs) = @_;

  #don't make the file twice
  unless( -f $un_xed_seqs ) {
    open my $xfile, $seqfile
      or confess "could not open '$seqfile' for reading ($!)";
    open my $nfile, ">$un_xed_seqs"
      or confess "could not open un-xed-seqs file '$un_xed_seqs' for writing";
    while (my $line = <$xfile>) {
      unless($line =~ /^\s*[>#]/) {    #don't munge identifier or comment (comment?) lines
	$line =~ tr/X/N/;
      }
      print $nfile $line;
    }
    close $nfile;
    close $xfile;
  }

  return $un_xed_seqs;
}

#get the genomic file to run this against
sub _seq_file {
  my ($self,$batch,$seqname) = @_;
  my ($seq_file) = $batch->pipeline->analysis('seq')->files_for_seq($batch,$seqname);
  $seq_file && -f $seq_file
    or die "expected sequence file '$seq_file' not found";
  return $seq_file;
}

=head1 AUTHOR(S)

Robert Buels

=cut

###
1;#do not remove
###
