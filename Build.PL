use strict;
use Module::Build;

my $build = Module::Build->new(
    dist_name     => 'ITAG',
    dist_author   => 'Sol Genomics Network',
    dist_abstract => 'data model and management scripts for ITAG pipeline',
    license       => 'perl',
    create_makefile_pl => 'passthrough',

    # current version of our core libraries distribution
    dist_version => '0.1',

    recursive_test_files => 1,

    # and now the list of perl module dependencies
    requires => {
        perl => '5.10',
        # right now, we have to depend on whole bioperl, like this
#        'Bio::Root::Version' => '1.006000',
# 	'Bio::FeatureIO' => 0,
# 	'Bio::Index::Fasta' => 0,
# 	'Bio::PrimarySeq' => 0,
# 	'Bio::SeqIO' => 0,
# 	'Bio::Tools::tRNAscanSE' => 0,

# 	'CXGN::BioTools::AGP' => 0,
# 	'CXGN::BlastDB' => 0,
# 	'CXGN::Config' => 0,
# 	'CXGN::DB::Connection' => 0,
# 	'CXGN::DB::DBICFactory' => 0,
# 	'CXGN::DB::GFF::Versioned' => 0,
# 	'CXGN::DB::Ima' => 0,
# 	'CXGN::IndexedLog' => 0,
# 	'CXGN::ITAG::SeqSource::TomatoContigs' => 0,
# 	'CXGN::Marker' => 0,
# 	'CXGN::Page::FormattingHelpers' => 0,
# 	'CXGN::TomatoGenome::BACPublish' => 0,
# 	'CXGN::TomatoGenome::ChromosomeAssemblies' => 0,
# 	'CXGN::Tools::Class' => 0,
# 	'CXGN::Tools::File' => 0,
# 	'CXGN::Tools::Identifiers' => 0,
# 	'CXGN::Tools::List' => 0,
# 	'CXGN::Tools::Run' => 0,
# 	'CXGN::Tools::Script' => 0,
# 	'CXGN::Tools::Wget' => 0,
# 	'CXGN::VHost' => 0,

	'Class::Accessor' => 0,
	'Class::Accessor::Fast' => 0,
	'Class::Data::Inheritable' => 0,
	'Digest::MD5' => 0,
	'File::NFSLock' => 0,
	'File::Temp' => 0,
	'Hash::Util' => 0,
	'List::MoreUtils' => 0,
	'Mail::Sendmail' => 0,
	'Memoize' => 0,
	'Module::Find' => 0,
        'Moose' => 0,
        'MooseX::Types::Path::Class' => 0,
        'namespace::autoclean' => 0,
	'Path::Class::File' => 0,
	'Scalar::Util' => 0,
	'Storable' => 0,
        'Template' => 0,
        'Template::Plugin::Comma' => 0,
	'Tie::Function' => 0,
	'URI::Escape' => 0,
    },

    build_requires => {
        'Capture::Tiny'  => 0,
	#'CXGN::DB::Connection' => 0,
	#'CXGN::Tools::Run' => 0,
	'File::Temp'  => 0,
        'Path::Class' => 0,
	'Test::More'  => 0,
	'Test::Warn'  => 0,
    },
);
$build->create_build_script;

