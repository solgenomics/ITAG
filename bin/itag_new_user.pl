#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

#use Data::Dumper;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script <username> <groupname>

  Make a new ITAG pipeline user and group, and add the itagpipeline
  user to the new group.  Lets you pick a random password for the user.

  Options:

  none

EOU
}
sub HELP_MESSAGE {usage()}

sub psystem(@) {
  print join(' ',@_),"\n";
  system @_;
  die "$_[0] failed: $!" if $CHILD_ERROR;
}

die "you probably want to run this as root\n"
  unless $EFFECTIVE_USER_ID == 0;

our %opt;
getopts('',\%opt) or usage();

my ($user,$group) = @ARGV;
$group ||= $user;
$user && $group or usage;

$user =~ /^[a-z][a-z\d]+$/
  or die "invalid username $user\n";
$group =~ /^[a-z][a-z\d]+$/
  or die "invalid groupname $group\n";

#run pwgen
psystem 'pwgen';

#make the group
psystem( 'addgroup', $group );

psystem( 'adduser',
	'--shell' => '/usr/sbin/scponlyc',
	'--ingroup' => $group,
        '--home' => "/data/shared/tomato_genome/itagpipeline//home/$user",
	$user
      );

#put a symlink into the user's home dir
system ln => '-s', '../../itag' => "/data/shared/tomato_genome/itagpipeline//home/$user/itag";

#add itagpipeline to the group
psystem( 'adduser','itagpipeline',$group);

#add this user to the itag_all group
psystem( 'adduser',$user,'itagall');

#refresh the /etc/passwd in the chroot jail, it must list all the itag users
#in order for scp to work
psystem( "grep itag /etc/passwd > /data/shared/tomato_genome/itagpipeline/etc/passwd" );

#add the user to the sshd_config authorized users
my @lines = do {
  open my $sshd_conf,'/etc/ssh/sshd_config';
  <$sshd_conf>
};
open my $sshd_conf,">/etc/ssh/sshd_config"
  or die "could not open /etc/ssh/sshd_config for writing: $!";
foreach (@lines) {
  if(/^\s*AllowUsers\s/ && ! /\b$user\b/) {
    chomp;
    my @f = split;
    print $sshd_conf join(' ',@f,$user),"\n";
  } else {
    print $sshd_conf $_;
  }
}

#now restart sshd
system '/etc/init.d/ssh' => 'restart';

