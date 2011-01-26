#!/bin/sh
cp TEMPLATE $1.pm;
perl -i -pe "s/ANALNAME/$1/" $1.pm;
