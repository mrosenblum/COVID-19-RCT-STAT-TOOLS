__README for ordinal simulation code__

# R package 

The `drord` package can be installed from GitHub.

> devtools::install_github("benkeser/drord")

# Overview 

The code for this simulation were run on a 
Linux system with a Slurm Workload Manager. The scripts can be 
executed with different commands passed in from the system 
environment. We have also included a shell script (sce.sh) that shows the 
workflow of executing the R scripts. However, note that some of 
the options for batching the jobs are native to the host system 
they were executed on and thus will error if executed on other 
systems. To read more about this work flow see 
[this page](https://github.com/FredHutch/slurm-examples/tree/master/centipede). 

The basic idea is that the bash script submits a job that first
computes how many jobs will be executed (submitting `ordinal_simulation.R` script 
with argument `"listsize"`). Then the `ordinal_simulation.R` script is submitted 
with argument `"prepare"`, which in this case, does nothing. 
Then `ordinal_simulation.R` is submitted as a slurm array
with argument `"run"`, which executes a single iteration of the 
simulation. Finally `ordinal_simulation.R` is submitted with argument `"merge"`, 
which merges the results of the individual jobs, formats the output and
produces the tables included in the manuscript and supplement. 

# Questions

There is a substantial amount of code associated with this project and
there was a significant amount of computational burden in executing the
production size jobs. The code submitted along with the manuscript needs 
to be modified in some places in order to ease its use on different systems. 
I have not debugged the code across all different systems and I 
cannot guarantee that the code will run error-free on new systems. If you come 
across any issues, [please reach out to me by email](benkeser@emory.edu) 
and I am happy to help sort them out. 