package CXGN::ITAG::CmdLine::Command::make_release;
sub abstract {
    'construct a genome release from the given pipeline and batch numbers' }

use Moose;
use namespace::autoclean;
use autodie qw/ :all /;
extends qw(CXGN::ITAG::CmdLine::Command);

use Carp;
use POSIX;

use Moose::Util::TypeConstraints;
use MooseX::Types::Moose qw/ ArrayRef Int Str Bool /;
use Template;

##### attributes

has 'version' => (
    documentation => 'ITAG pipeline version to use',
    traits        => [qw(Getopt)],
    cmd_aliases   => 'v',
    is            => 'ro',
    isa           => Int,
    predicate     => 'has_version',
);

has 'pipeline' => (
    is  => 'ro',
    isa => 'CXGN::ITAG::Pipeline',
    lazy_build => 1,
   ); sub _build_pipeline {
       my $self = shift;
       return CXGN::ITAG::Pipeline->open(
           basedir => $self->itag_dir,
           ( $self->has_version ? (version => $self->version) : () ),
          );
   }

has 'target_dir' => (
    documentation => 'target directory for file output',
    traits        => [qw(Getopt)],
    cmd_aliases   => 't',

    is       => 'ro',
    isa      => 'Path::Class::Dir',
    coerce   => 1,
    required => 1,
   );

has 'release_num' => (
    documentation => 'release number for the release we are making',
    traits        => [qw(Getopt)],
    cmd_aliases   => 'r',

    is       => 'ro',
    isa      => Int,
    required => 1,
   );


{ my $rt = subtype __PACKAGE__.'::Type::ReleaseType'
      => as Str
      => where { $_ eq 'dev' || $_ eq 'pre' || $_ eq 'normal' };

  coerce $rt, from Str, via { lc $_ };

  has 'release_type' => (
      documentation => "type of release, either 'dev', 'pre', or 'normal', default normal",
      traits        => [qw(Getopt)],
      cmd_aliases   => 'P',

      is      => 'ro',
      isa     => $rt,
      default => 'normal',
     );

}

has 'dev_release' => (
    documentation => 'release number for the release we are making',
    traits        => [qw(Getopt)],
    cmd_aliases   => 'P',

    is      => 'ro',
    isa     => Bool,
    default => 0,
   );

has 'release' => (
    is  => 'ro',
    isa => 'CXGN::ITAG::Release',
    lazy_build => 1,
   ); sub _build_release {
       my $self = shift;
       # set the root dir to our target directory so the object will make it in the right place#
       CXGN::ITAG::Release->releases_root_dir( $self->target_dir->stringify );
       CXGN::ITAG::Release->new( releasenum => $self->release_num,
                                 pre   => $self->release_type eq 'pre',
                                 devel => $self->release_type eq 'dev',
                                );
   }

{ my $bns = subtype as ArrayRef[Int];
  coerce $bns, from Str, via { [ / (\d+) /xg ]}; # just extracts all runs of digits

  # this is an arrayref of integers
  has 'batch_nums' => (
      documentation => 'the batch numbers to include in this release, comma-separated',
      traits        => [qw(Getopt)],
      cmd_aliases   => 'b',

      is       => 'ro',
      isa      => $bns,
      required => 1,
      coerce   => 1,
     );
}


sub batches {
    my $self = shift;
    return map {
        $self->pipeline->batch($_)
            || die "batch number $_ does not exist in pipeline ".$self->pipeline->version."\n";
    } @{ $self->batch_nums };
}

########## main methods

with
    __PACKAGE__.'::Validation', #< validate_args and all the arg checking is in here
    __PACKAGE__.'::Readme',     #< readme generation and template
    __PACKAGE__.'::Dumping',    #< the actual data dumping
    ;

sub execute {
    my ( $self, $opt, $args ) = @_;

    # all this will do is assemble the release files in the target dir
    # with the given release tag

    $self->dump_data;

    #write the readme file
    $self->write_readme;
}
