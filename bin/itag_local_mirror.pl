#!/usr/bin/perl

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use warnings;

use English;
use Carp;
use FindBin;

use Getopt::Std;
use Pod::Usage;

#use Data::Dumper;


######### DEFAULTS ###########

use constant ITAG_PATH_DEFAULT => '/data/shared/tomato_genome/itagpipeline';
use constant DEFAULT_HOST => 'eggplant.sgn.cornell.edu';

use constant REMOTE_ITAG_PATH => '/data/shared/tomato_genome/itagpipeline';
use constant ITAG_DIR_NAME => 'itag';

##############################

our %opt;
getopts('p:h:',\%opt) or pod2usage(1);

$opt{p} ||= ITAG_PATH_DEFAULT;
-d $opt{p} or die "directory $opt{p} does not exist\n";
-w $opt{p} or die "directory $opt{p} is not writable\n";

$opt{h} ||= DEFAULT_HOST;

for my $dep ( qw/ rsync tar ssh / ) {
  `which $dep` or die "$dep not found in path, please install it\n";
}

my $glob = "$opt{p}/itag/pipe*/batch*/seq/*.fasta";

if( glob $glob ) {
  exec 'rsync',
    '-avz',
    '--no-p',
    '--delete',
    "$opt{h}:".REMOTE_ITAG_PATH.'/'.ITAG_DIR_NAME,
    $opt{p};
} else {
  system "ssh $opt{h} 'cd ".REMOTE_ITAG_PATH."; tar -czf - ".ITAG_DIR_NAME."' | tar -C $opt{p} -xzvf - \n";
}

exit;

__END__

=head1 NAME

itag_local_mirror.pl - use rsync to make or update a local mirror of
itag pipeline files

=head1 SYNOPSIS

  itag_local_mirror.pl [options]

  Options:

   -p <path>
      set the base path for the local ITAG mirror.
      Defaults to /data/shared/tomato_genome/itagpipeline

   -h <host>
      host to rsync from.
      Defaults to upload.sgn.cornell.edu

=head1 MAINTAINER

Robert Buels

=head1 AUTHOR(S)

Robert Buels, E<lt>rmb32@cornell.eduE<gt>

=head1 COPYRIGHT & LICENSE

Copyright 2009 The Boyce Thompson Institute for Plant Research

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut
