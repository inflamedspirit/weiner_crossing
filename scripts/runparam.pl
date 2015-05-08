#!/usr/bin/perl

use strict;

# Output files
our $bin_dir;
our $stdoutfile;
our $stderrfile;
our $run_number;
our $parametermotbase;

# Load Parameters
if(@ARGV != 1 ){
    die "Needs parameter.gen file as argument.";
}
my $paramfile = $ARGV[0];
require $paramfile;

# Launch Mot
system("$bin_dir/weiner_crossing $parametermotbase.$run_number 1> $stdoutfile 2> $stderrfile ");

