use strict;
use Module::Build;

my $build = Module::Build->new(
    dist_name     => 'ITAG-analyses',
    dist_author   => 'Sol Genomics Network',
    dist_abstract => 'analysis modules for ITAG pipeline',
    license       => 'perl',
    create_makefile_pl => 'passthrough',

    # current version of our core libraries distribution
    dist_version => '0.1',

    recursive_test_files => 1,

    # and now the list of perl module dependencies
    requires => {
        perl              => '5.10.0',
        'File::Temp'      => 0,
        'List::MoreUtils' => 0,
        'Memoize'         => 0,
        # right now, we have to depend on whole bioperl, like this
       'Bio::Root::Version' => '==1.005002102',
#  	'Bio::FeatureIO' => 0,
# 	'Bio::SeqIO' => 0,
# 	'Bio::Tools::tRNAscanSE' => 0,

# 	'CXGN::BioTools::AGP' => 0,
# 	'CXGN::BlastDB' => 0,
# 	'CXGN::DB::Ima' => 0,
# 	'CXGN::ITAG::Pipeline::Analysis' => 0,
# 	'CXGN::Marker' => 0,
# 	'CXGN::Tools::Class' => 0,
# 	'CXGN::Tools::Identifiers' => 0,
# 	'CXGN::Tools::List' => 0,
# 	'CXGN::Tools::Run' => 0,
# 	'CXGN::Tools::Wget' => 0,
    },

    build_requires => {
        'File::Temp' => 0,
        'Test::Most' => 0,
    },
);
$build->create_build_script;
my $build = Module::Build->new(
    dist_name     => 'ITAG-analyses',
    dist_author   => 'Sol Genomics Network',
    dist_abstract => 'analysis modules for ITAG pipeline',
    license       => 'perl',
    create_makefile_pl => 'passthrough',

    # current version of our core libraries distribution
    dist_version => '0.1',

    recursive_test_files => 1,

    # and now the list of perl module dependencies
    requires => {

        # right now, we have to depend on whole bioperl, like this
       'Bio::Root::Version' => '==1.005002102',
#  	'Bio::FeatureIO' => 0,
# 	'Bio::SeqIO' => 0,
# 	'Bio::Tools::tRNAscanSE' => 0,

# 	'CXGN::BioTools::AGP' => 0,
# 	'CXGN::BlastDB' => 0,
# 	'CXGN::DB::Ima' => 0,
# 	'CXGN::ITAG::Pipeline::Analysis' => 0,
# 	'CXGN::Marker' => 0,
# 	'CXGN::Tools::Class' => 0,
# 	'CXGN::Tools::Identifiers' => 0,
# 	'CXGN::Tools::List' => 0,
# 	'CXGN::Tools::Run' => 0,
# 	'CXGN::Tools::Wget' => 0,

	'File::Temp' => 0,
	'List::MoreUtils' => 0,
	'Memoize' => 0,
    },

    build_requires => {
	'File::Temp' => 0,
	'Test::More' => 0,
    },
);
$build->create_build_script;

my $build = Module::Build->new(
    dist_name     => 'ITAG-analyses',
    dist_author   => 'Sol Genomics Network',
    dist_abstract => 'analysis modules for ITAG pipeline',
    license       => 'perl',
    create_makefile_pl => 'passthrough',

    # current version of our core libraries distribution
    dist_version => '0.1',

    recursive_test_files => 1,

    # and now the list of perl module dependencies
    requires => {
        autodie => 0,

        # right now, we have to depend on whole bioperl, like this
       'Bio::Root::Version' => '==1.005002102',
#  	'Bio::FeatureIO' => 0,
# 	'Bio::SeqIO' => 0,
# 	'Bio::Tools::tRNAscanSE' => 0,

# 	'CXGN::BioTools::AGP' => 0,
# 	'CXGN::BlastDB' => 0,
# 	'CXGN::DB::Ima' => 0,
# 	'CXGN::ITAG::Pipeline::Analysis' => 0,
# 	'CXGN::Marker' => 0,
# 	'CXGN::Tools::Class' => 0,
# 	'CXGN::Tools::Identifiers' => 0,
# 	'CXGN::Tools::List' => 0,
# 	'CXGN::Tools::Run' => 0,
# 	'CXGN::Tools::Wget' => 0,

	'File::Temp' => 0,
	'List::MoreUtils' => 0,
	'Memoize' => 0,
    },

    build_requires => {
        aliased => 0,
	'File::Temp' => 0,
	'Test::More' => 0,
    },
);
$build->create_build_script;
