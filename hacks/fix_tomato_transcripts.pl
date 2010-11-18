#!/usr/bin/env perl
use strict;
use warnings;

use Data::Dump 'dump';
use Storable 'dclone';

use List::Util qw/ min max /;

# hash the lines by Target
my %seen;
my %lines;
while( <> ) {
    unless( /^(\S+)\tITAG_transcripts_\S+.+Target=(\S+)/ ) {
        next;
    }
    #warn "$1:$2 $_";
    next if $seen{$_}++;

    my @f = split /\t/, $_, 9;
    push @{ $lines{"$1:$2"} }, \@f;
}

#warn dump( [values %lines] );

my @groups = map {
    my ($pts,$pte);
    my $pfe;
    my $curr_group = [];
    my @groups = ($curr_group);
    for my $l ( @$_ ) {
        my ($fs,$fe) = @{$l}[3,4];
        my $strand = $l->[6];
	my ($ts,$te) = $l->[8] =~ /Target=\S+ (\d+) (\d+)/ or die;

	if(    $pte
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

#warn dump( \@groups );

# make feature groups out of them
for my $g ( @groups ) {
    if( @$g > 1 ) {
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
        print(join "\t", @$_) for $superfeature, @$g;
    }
    else {
        print join "\t", @$_ for  @$g;
    }
}
