#!/usr/bin/perl

#
# ./scripts/setup_run.pl
# This script has two functions.
#
# (1) it generates $parameterbase.$j files with are used by runparam.pl script to launch the 
#     main simulation. These files contain the information on where stderr and stdout files
#     files are saved to (which are all of the form basename.jobnumber), as well as where
#     the $parametermotbase.$j file is stored.
#
# (2) it generates $parametermotbase.$j files which are used by the mot simulation (mot) and
#     store all the parameters used by the simulation in a simple format for loading.
#
# setup_run.pl should be called before launch_script.pl to run paramfiles.
#

use strict;
use warnings;

# File and job variables to be loaded:
our $num_jobs;

our $head_dir;
our $bin_dir;
our $data_dir;

our $run_title;

our $run_dir;
our $parameter_dir;
our $stderr_dir;
our $stdout_dir;
our $analysis_dir;

our $stdoutbase;
our $stderrbase;
our $parameterbase;
our $parametermotbase;

# Base parameters to be loaded:

# Simulation Parameters
our $seed;
our $num_steps;
our $num_walks;
our $discrete_steps;

# Regular Bin Parameters
our $logarithmic_binning;
our $bin_start;
our $bin_end;
our $bin_num;

# Crossings Parameters
our $crossings_logarithmic_binning;
our $crossings_bin_start;
our $crossings_bin_end;
our $crossings_bin_num;
our $crossings_max;

# Adjusted Parameters:
our $seed_increment;

our $do_step_increment;
our $step_multiplier;

our $do_walk_increment;
our $walk_multiplier;


# Load the PARAMETERS file:
if( @ARGV != 1 ) {
    die("Needs PARAMETERS.pl file as argument.");
}

my $paramfile = $ARGV[0];
require $paramfile;

# Create directories and run, unless directory exists.
if( -d "$run_dir" ) {
    die "Error: data directory exists. Delete $run_dir or rename \$run_title to something other than $run_title.\n";
} else { 
    mkdir("$run_dir") or die "Could not make directory";
    mkdir("$parameter_dir") or die "Could not make directory";
    mkdir("$stdout_dir") or die "Could not make directory";
    mkdir("$stderr_dir") or die "Could not make directory";
}

# Main loop:
my $j;
for ($j=1; $j<=$num_jobs; $j++) {

# Adjusted parameters
    $seed = $seed + $seed_increment;

    if( $do_step_increment == 1 ){
	$num_steps*=$step_multiplier;
    }
    if( $do_walk_increment == 1 ){
	$num_walks*=$walk_multiplier;
    }



# Print the $parameberbase.$j file
    open(my $paramfh, ">$parameterbase.$j");
    print $paramfh "# file automatically generated by makeparameters.pl\n";
    print $paramfh "\$bin_dir='$bin_dir';\n";
    print $paramfh "\$stdoutfile='$stdoutbase.$j';\n";
    print $paramfh "\$stderrfile='$stderrbase.$j';\n";
    print $paramfh "\$parametermotbase = '$parametermotbase';\n";
    print $paramfh "\$run_number='$j';\n";
    print $paramfh "\n";
    print $paramfh "1\n";
    print $paramfh "\n";
    close($paramfh);

# Print the $paramebermotbase.$j file
    open($paramfh, ">$parametermotbase.$j");
    print $paramfh "seed                              $seed\n";
    print $paramfh "num_steps                         $num_steps\n";
    print $paramfh "num_walks                         $num_walks\n";
    print $paramfh "discrete_steps                    $discrete_steps\n";
    print $paramfh "logarithmic_binning               $logarithmic_binning\n";
    print $paramfh "bin_start                         $bin_start\n";
    print $paramfh "bin_end                           $bin_end\n";
    print $paramfh "bin_num                           $bin_num\n";
    print $paramfh "crossings_logarithmic_binning     $crossings_logarithmic_binning\n";
    print $paramfh "crossings_bin_start               $crossings_bin_start\n";
    print $paramfh "crossings_bin_end                 $crossings_bin_end\n";
    print $paramfh "crossings_bin_num                 $crossings_bin_num\n";
    print $paramfh "crossings_max                     $crossings_max\n";
    
    close($paramfh);

}

#return 0; #this is for success?
