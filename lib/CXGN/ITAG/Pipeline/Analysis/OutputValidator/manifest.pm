package CXGN::ITAG::Pipeline::Analysis::OutputValidator::manifest;
use strict;
use warnings;

use base 'CXGN::ITAG::Pipeline::Analysis::OutputValidator';
use File::Basename;
sub is_intensive {0}
sub run_online {
  my ($self,$analysis,$batch) = @_;
  my $mfile = $analysis->_filename('manifest',$batch);
  my $dir = $analysis->work_dir($batch);
  if(-f $mfile) {
    my @errors;
    CORE::open my $m,$mfile
      or die "could not open manifest file $mfile: $!";
    while(my $line = <$m>) {
      my ($filename,$msize) = split /\s+/,$line;
      my $full = "$dir/$filename";
      unless(-f $full) {
	push @errors,"$filename : mentioned in manifest file, but is not present";
	next;
      }
      my $fsize = -s $full;
      unless($fsize == $msize) {
	push @errors,"$filename : size mismatch (manifest: $msize bytes, actual: $fsize bytes)";
	next;
      }
    }
    close $m;
    return @errors;
  }
  return;
}


1;
