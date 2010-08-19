package CXGN::ITAG::Pipeline::Analysis::OutputValidator::fasta_file_format;
use base 'CXGN::ITAG::Pipeline::Analysis::OutputValidator';
use File::Basename;
use Carp;

sub is_intensive { 1 }

sub validates_file {
  my ($self,$filename) = @_;
  return $filename =~ /\.fasta$/;
}

sub run_offline {
  my ($self,$analysis,$batch) = @_;
  my @fasta_files = $self->files_to_validate($analysis,$batch)
    or return; #< don't do any of this unless we have some fasta files

  my @file_errors;
  my %failures;
 FILE:
  foreach my $file (@fasta_files) {
    my $bn = basename($file);
#    warn "validating $bn...\n";
    my $outbn = $self->cache_dir($analysis,$batch)."/$bn";
    my $reportfile = "$outbn.report";

    @file_errors = ();
    open my $f,$file or confess "$! opening $file for reading";
    my %not_iupac_pats = ( dna     => qr/([^ACGTURYKMSWBDHVN]+)/i,
			   protein => qr/([^GAVLIPFYCMHKRWSTDENQBZ\.X\*]+)/i,
			   rna     => qr/([^ACGTURYKMSWBDHVN]+)/i,
			 );
    my @possible_alphabets = keys %not_iupac_pats;
  LINE:
    while( my $line = <$f> ) {
      if( $line =~ />/ ) {
	#there can be no spaces on either side of a >
	if ( $line =~ /.>/ || $line =~ />\s/ ) {
	  push @file_errors, "line $.: whitespace not allowed next to '>'";
	}
	$line =~ />\s*\S+(\n|\s+.+\n)/
	  or push @file_errors, "line $.: malformed definition line";
      } else {
	chomp $line;
	$line =~ /\S/ or push @file_errors, "line $.: blank lines not allowed in FASTA files";
	my @new_possible_alphabets = grep {  $line !~ $not_iupac_pats{$_} } @possible_alphabets;
	unless(@new_possible_alphabets) {
	  push @file_errors, "line $.: invalid characters in sequence";
	} else {
	  @possible_alphabets = @new_possible_alphabets;
	}
      }
    }
    close $f;
  } continue {
    #touch the report file, so we know the validation has run
    CORE::open my $report,'>',$self->report_filename($analysis,$batch,$file);
    if(@file_errors) {
      my $bn = basename($file);

      # record that this file had a failure
      $failures{$file} = time;

      # dump the errors into the report file
      print $report map "$_\n",@file_errors;
    }
  }

  $self->write_failures( $analysis, $batch, \%failures );

  return;
}


1;
