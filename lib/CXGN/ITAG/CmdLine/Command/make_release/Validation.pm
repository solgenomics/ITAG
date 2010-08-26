package CXGN::ITAG::CmdLine::Command::make_release::Validation;
use Moose::Role;

requires qw(
            batches
            target_path
            pipeline
            vsay
            vprint
           );

sub validate_args {
    my ( $self, $opt, $args ) = @_;

    $self->batches; #< die fast if we don't have all our batches

    { my $target_path = $self->target_path;
      -d $target_path or die "target path $target_path does not exist\n";
      -w $target_path or die "target path $target_path is not writable\n";
    }

    $self->_check_seq_inputs;
    $self->_check_required_analyses;

    $self->release->mkdir;
    -d $self->release->dir or die "$! creating build directory ".$release->dir;
    -w $self->release->dir or die "build directory ".$release->dir." is not writable\n";
}

########################

sub _check_seq_inputs {
    my ( $self ) = @_;

    for my $batch ( $self->batches ) {
        for my $aname (qw( renaming seq )) {
            for my $file ( $self->pipeline->analysis($aname)->->files_for_seq( $batch, $ctg->{name} ) ) {
                #check that all the files we need are readable
                -r or die "$_ file not readable\n";
            }
        }
    }
}


sub _check_required_analyses {
    # flag which analyses are required to be done in order for a batch to
    # be included in the release
    my @required_analyses = (
        #( map $_->{name}, grep {$_->{name} =~ /^blastp_/i} values %a), #< all the blastp analyses are required
        'renaming',
        'trnascanse',
        'sgn_markers',
        'infernal',
        'transcripts_sol',
        'transcripts_tomato',
       );
    foreach my $batch ( $self->batches ) {
        $self->vsay("checking analyses in batch $batch");
        foreach my $aname ( @required_analyses ) {
            $self->vprint("  checking $aname ... ");
            my $a = $self->pipeline->analysis( $aname )
                or die "unknown analysis $_";
            unless( $a->status( $batch ) eq 'done' ) {
                die "analysis $aname not done in batch $batch, cannot make release.\n";
            }
            $self->vsay('done.');
        }
    }
}


no Moose::Role;
1;
