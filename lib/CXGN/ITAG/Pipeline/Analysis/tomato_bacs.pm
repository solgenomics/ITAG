package CXGN::ITAG::Pipeline::Analysis::tomato_bacs;
use strict;
use warnings;

use autodie ':all';

use base qw/CXGN::ITAG::Pipeline::Analysis::mummer_base/;

use List::MoreUtils 'uniq';

use CXGN::Tools::Wget qw/ wget_filter /;

sub mummer_params {

    return (
        [ qw( -mum -b -n -L ) ],

        {
            -l => 20_000,
        }
       );
}

sub gff3_source {
    'ITAG_tomato_bacs'
}

sub query_file_url {
    'ftp://ftp.solgenomics.net/genomes/Solanum_lycopersicum/bacs/curr/bacs.seq';
}

sub run_mummer {
    my $self = shift;
    my $seq_file = $_[2];
    $self->_load_deflines( $seq_file );
    $self->SUPER::run_mummer(@_);
}

# munge gff3 to add aliases to the attrs
sub munge_gff3 {
    my ( $class, $args, $gff3, $attrs ) = @_;
    my $name = $attrs->{Name};
    my $aliases = $class->_get_aliases( $name );
    if( @$aliases ) {
        $attrs->{Name} = shift @$aliases;
        $attrs->{Alias} = join ',', @$aliases if @$aliases;
        if( my $d = $class->_get_defline( $name ) ) {
            $attrs->{Note} ||= $d;
        }
    }
}

sub _accessions_file {
    my ( $class ) = @_;
    return wget_filter(
        'ftp://ftp.solgenomics.net/genomes/Solanum_lycopersicum/bacs/curr/bacs_accessions.txt'
            => $class->local_temp('accession_list.txt'),
       );
}

my %deflines;
sub _load_deflines {
    my ( $class, $seqfile ) = @_;
    open my $s, '<', $seqfile;
    while( my $l = <$s> ) {
        my ( $bac, $defline ) = $l =~ /^>(\S+)\s*(.+)/
            or next;
        chomp $defline;
        $deflines{$bac} = $defline;
    }
}
sub _get_defline {
    my ( $class, $bac ) = @_;
    $class->_load_deflines( $class->query_file ) unless %deflines;
    return $deflines{$bac};
}

my %aliases;
sub _load_aliases {
    my ( $class ) = @_;
    my $accessions = $class->_accessions_file;

    open my $accs, '<', $accessions;
    while( my $l = <$accs> ) {
        chomp $l;
        my ( $bac, $acc ) = split /\s+/, $l;
        ( my $unv_bac = $bac ) =~ s/(\.\d+)+$//;
        ( my $unv_acc = $acc ) =~ s/(\.\d+)+$//;
        my @aliases = uniq $bac, $unv_bac,  $acc, $unv_acc;
        $aliases{$bac} = \@aliases;
    }
}
sub _get_aliases {
    my ( $class, $bacname ) = @_;
    $class->_load_aliases unless %aliases;
    return $aliases{$bacname} if $aliases{$bacname};
    if( $bacname =~ s/(\.\d+)+$// ) {
        return [ $bacname ];
    } else {
        return [];
    }
}
sub _all_aliases {
    shift->_load_aliases unless %aliases;
    return %aliases;
}

1;
