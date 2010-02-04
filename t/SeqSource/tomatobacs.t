#!/usr/bin/perl
use strict;
use warnings;
use English;

use Data::Dumper;

use CXGN::DB::Connection;

BEGIN {
  our @tests = (
		[ phase => 2 ],
		[unfinished => 1, phase => 2],
		[],
		[unfinished => 1],
	       );
}
our @tests;
use Test::More tests => 1+2*scalar @tests;
BEGIN {
  use_ok(  'CXGN::ITAG::SeqSource::TomatoBACs'  )
    or BAIL_OUT('could not include the module being tested');
}

my $prev_return = 0;
foreach my $pset ( @tests ) {
  my $s = CXGN::ITAG::SeqSource::TomatoBACs->new( @$pset );
  isa_ok( $s => 'CXGN::ITAG::SeqSource::TomatoBACs' );

  my $last_last = $CXGN::ITAG::SeqSource::TomatoBACs::LAST_QUERY;
  #diag $last_last;
  my @names =  $s->all_seq_names;
  local $Data::Dumper::Indent = 0;
  local $Data::Dumper::Terse = 1;
  ok @names >= $prev_return, Dumper($pset)." params yields ".@names.", which is at least as many seqs as the last run ($prev_return)"
    or diag "PREVIOUS Q: $last_last\nAND THIS Q: $CXGN::ITAG::SeqSource::TomatoBACs::LAST_QUERY";
  #diag map "$_\n",@names;
  #diag "(".scalar(@names)." seqs)\n";
  $prev_return = @names;
}
