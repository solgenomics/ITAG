package CXGN::ITAG::Pipeline::Analysis::OutputValidator::md5sums;

use base 'CXGN::ITAG::Pipeline::Analysis::OutputValidator';

use File::Basename;

sub is_intensive {1}

# OLD DOCS FOR THIS METHOD
# =head2 check_md5sum

#   Usage: print 'yep' if $an->check_md5sum;
#   Desc : check the md5sums file in the dir,
#          if included
#   Args : none
#   Ret  : ('none') if no md5sums file,
#          list of error strings if mismatches found,
#          empty list if md5sums all check OK
#   Side Effects: looks in the filesystem, runs the
#                 'md5sum' program.

# =cut

sub run_online {
  my ($self,$analysis,$batch) = @_;
  my $mfile = $analysis->_filename('md5sum',$batch);
  my $dir = $analysis->work_dir($batch);
  if(-f $mfile) {
    my @errors;
    CORE::open my $m,$mfile
      or die "could not open md5sum file $mfile: $!\n";
    while(my $line = <$m>) {
      my ($msum,$filename) = split /\s+/,$line;
      my $full = "$dir/$filename";
      unless(-f $full) {
	push @errors,"$filename : mentioned in md5sums file, but is not present\n";
	next;
      }
      CORE::open my $fh, $full
	or do { push @errors, "$filename : cannot open file for reading\n";
		next;
	      };
      my $md5 = Digest::MD5->new;
      $md5->addfile($fh);
      my $fsum = $md5->hexdigest;
      unless(lc($fsum) eq lc($msum)) {
	push @errors, "$filename : MD5 sum mismatch\n";
	next;
      }
      close $fh;
    }
    close $m;
    return @errors;
  }
  return;
}

sub validates_file { 1 } #< can validate all files

sub needs_update {
  my ($self,$analysis,$batch) = @_;

  my $mfile = $analysis->_filename('md5sum',$batch);
  my $mfile_time = (stat $mfile)[9];
  return 0 unless $mfile_time;

  my @files = $self->files_to_validate($analysis,$batch);
  my $min_rs = time;
  foreach my $file (@files) {
      my $reportfile = $self->report_filename($analysis,$batch,$file);

      my $rs = (stat($reportfile))[9];
      my $fs = (stat($file))[9];
      return 1 if
	   !$rs                # file has not been checked
	|| $rs < $fs           # file is newer than its report
	|| $mfile_time >= $rs; # md5sum file is newer than the report that checks it

  } 

  return 0;
}

1;
