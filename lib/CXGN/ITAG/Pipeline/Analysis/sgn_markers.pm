package CXGN::ITAG::Pipeline::Analysis::sgn_markers;
use strict;
use warnings;
use English;
use Carp;

use Bio::SeqIO;

use CXGN::Marker;
use CXGN::Tools::Wget qw/ wget_filter /;

=head1 NAME

CXGN::ITAG::Pipeline::Analysis::sgn_markers - pipeline analysis to
produce gapped alignments of SGN marker sequences

=head1 BASE CLASS(ES)

L<CXGN::ITAG::Pipeline::Analysis>

=cut

use base qw/CXGN::ITAG::Pipeline::Analysis::gth_base/;

=head1 METHODS

implements all abstract methods defined in
L<CXGN::ITAG::Pipeline::Analysis>

=cut

sub _source {
    'SGN_Markers'
}

sub _cdna_file {
    my $class = shift;
    wget_filter('cxgn-resource://sgn_marker_seqs'
                => $class->local_temp('sgn_markers.fasta'),
               );
}

# cached piece of class data, a db handle
{ my $dbh;
  sub _dbh {
      $dbh ||= CXGN::DB::Connection->new;
  }
}
sub _process_gff3_feature {
    my ($class, $feature) = @_;

    my $annotation = $feature->annotation;
    my $id = ($annotation->get_Annotations('Target'))[0]->target_id;

    my $marker =  $id =~ /SGN-M(\d+)/
        ? CXGN::Marker->new( $class->_dbh, $1 )
        : CXGN::Marker->new_with_name( $class->_dbh, $id );
    my @names = $marker->name_that_marker;
    push @names, 'SGN-M'.$marker->marker_id;
    my $primary_name = shift @names;

    $annotation->add_Annotation( Name  => Bio::Annotation::SimpleValue->new( -value => $primary_name ) );
    $annotation->add_Annotation( Alias => Bio::Annotation::SimpleValue->new( -value => $_ ) ) for @names;
    $annotation->add_Annotation( Note  => Bio::Annotation::SimpleValue->new( -value => "marker name(s): ".join(', ',$primary_name, @names) ) );

    return;
}


sub _parse_mode {
    'alignments_merged'
}

=head1 AUTHOR(S)

Robert Buels

=cut

###
1;#do not remove
###
