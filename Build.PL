use strict;
use warnings;
use Module::Build;

my $build = Module::Build->new(
    dist_name     => 'ITAG',
    dist_author   => 'Sol Genomics Network',
    dist_abstract => 'data model, management, and analysis modules for the ITAG pipeline',
    license       => 'perl',
    create_makefile_pl => 'passthrough',

    # current version of our core libraries distribution
    dist_version => '0.1',

    recursive_test_files => 1,

    # and now the list of perl module dependencies
    requires => {
        perl                         => '5.10.0',
        'File::Temp'                 => 0,
        'List::MoreUtils'            => 0,
        'Memoize'                    => 0,

        'Bio::GFF3'                  => '0.4',
        # right now, we have to depend on whole bioperl, like this
        'Bio::Root::Version'         => '1.005002102',
        'Class::Accessor'            => 0,
        'Class::Accessor::Fast'      => 0,
        'Class::Data::Inheritable'   => 0,
        'Config::INI::Reader'        => 0,
        'Config::INI::Writer'        => 0,
        'Digest::MD5'                => 0,
        'File::NFSLock'              => 0,
        'File::Temp'                 => 0,
        'File::Slurp'                => 0,
        'Hash::Util'                 => 0,
        'List::MoreUtils'            => 0,
        'Mail::Sendmail'             => 0,
        'Memoize'                    => 0,
        'Module::Find'               => 0,
        'Moose'                      => 0,
        'MooseX::Types::Path::Class' => 0,
        'namespace::autoclean'       => 0,
        'Path::Class::File'          => 0,
        'Scalar::Util'               => 0,
        'Statistics::Descriptive'    => 0,
        'Storable'                   => 0,
        'Template'                   => 0,
        'Template::Plugin::Comma'    => 0,
        'Tie::Function'              => 0,
        'URI::Escape'                => 0,
    },

    build_requires => {
        'Capture::Tiny' => 0,
        'File::Temp'    => 0,
        'Mock::Quick'   => 0,
        'Path::Class'   => 0,
        'Test::More'    => 0,
        'Test::Most'    => 0,
    },
);
$build->create_build_script;
