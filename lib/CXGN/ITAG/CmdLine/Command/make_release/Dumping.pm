package CXGN::ITAG::CmdLine::Command::make_release::Dumping;
use Moose::Role;

requires qw(
            pipeline
            release
           );


sub dump_data {
    my $self = shift;

    my @seqlist = $self->full_seqlist;
    my $output_files = $self->release->get_all_files;
    while( my ( $shortname, $file_info ) = each %$output_files ) {
        $self->find_dump_handler( $shortname, $file_info )
             ->dump( $self, \@seqlist );
    }
}

sub find_dump_handler {
    my ( $self, $shortname, $file_info ) = @_;
}


# object method to compare two gff3 lines
sub _gff3_cmp {
    # sort by reference sequence
    $_[1]->[0] cmp $_[2]->[0]
    # then starting coordinate
    || $_[1]->[3] <=> $_[2]->[3]
    # then analysis name
    || $_[1]->[1] cmp $_[2]->[1]
}


#############
# list of [ 'seqname', $batch_obj ] for all sequences we will be
# dumping, sorted by seq name descending
sub full_seqlist {
    my $self = shift;

    my @seqs =
        sort {
             $a->[0] cmp $b->[0]
          || die "duplicate sequence names '$a->[0]' found in batches $a->[1] and $b->[1]"
        }
        map  { my $batch = $_; map [$_,$batch], $batch->seqlist }
        $self->batches;
}



######## streaming stuff

# TODO: add merging and sorted-stream features to Data::Stream::Bulk,
# and refactor all this streaming stuff to use it

sub merge_sorted_streams {
    my ( $self, $cmp, @streams ) = @_;

    my @buffers = map little_buffer->new($_), @streams;

    return callback_stream->new(
        sub {
            @buffers =
                sort { $cmp->( $a->peek, $b->peek ) }
                grep { defined $_->peek }
                @buffers;

            return unless @buffers;
            return $buffers[0]->next;
        }
       );
}

package callback_stream;

sub new {
    my ( $class, $sub ) = @_;

    return bless $sub, $class;
}

sub get { shift->() }

package little_buffer;

sub new {
    my ($class,$thing) = @_;
    bless [ [] ,$thing ], $class;
}

sub next {
    my ($self) = @_;
    if( @{$self->[0]} ) {
        return shift @{ $self->[0] };
    } else {
        return $self->[1]->get;
    }
}

sub peek {
    my ($self) = @_;
    unless( @{$self->[0]} ) {
        @{$self->[0]} = $self->[1]->get;
    }
    return $self->[0][0];
}

