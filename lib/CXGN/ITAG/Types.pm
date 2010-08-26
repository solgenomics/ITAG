package CXGN::ITAG::Types;

use MooseX::Types;

use MooseX::Types::Moose qw/Int HashRef ArrayRef/;

class_type 'CXGN::ITAG::Pipeline::Batch';

subtype 'CXGN::ITAG::Types::BatchList',
    as ArrayRef['CXGN::ITAG::Pipeline::Batch'];




