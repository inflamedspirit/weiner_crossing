# File and job variables:

our $run_title = "decreasing_dt";
our $num_jobs = 10;

our $head_dir = "/home/wwe/weiner_crossing";

our $root_paramfile = "$head_dir/PARAMETERS.pl";

our $bin_dir = "$head_dir/bin";
our $script_dir = "$head_dir/scripts";
our $data_dir = "$head_dir/data";
our $run_dir = "$data_dir/$run_title";
our $parameter_dir = "$run_dir/parameters";
our $stderr_dir = "$run_dir/stderr";
our $stdout_dir = "$run_dir/stdout";

our $stdoutbase = "$stdout_dir/stdout";
our $stderrbase = "$stderr_dir/stderr";

our $parameterbase = "$parameter_dir/parameters.gen";
our $parametermotbase = "$parameter_dir/parameters.mot";

# Functionality Parameters
our $dry_run=0;

# Simulation Parameters
our $seed=2;
our $num_steps=10**1.0;
our $num_walks=10**1.0;

# Regular Bin Parameters
our $logarithmic_binning=1;
our $bin_start=0.001;
our $bin_end=10**9.0;
our $bin_num=60;

# Crossings Parameters
our $crossings_logarithmic_binning=1;
our $crossings_bin_start=10**(-7.0);
our $crossings_bin_end=10.0;
our $crossings_bin_num=30;
our $crossings_max=10;

# Adjusted Parameters:
our $seed_increment=1;

our $do_step_increment=1;
our $step_multiplier=10;

our $do_walk_increment=0;
our $walk_multiplier=10;


return 1;
