__README for time-to-event simulation code__

# R package 

The `survtmlerct` package can be installed from GitHub.

> devtools::install_github("idiazst/survtmlerct")

# Overview 

The code for this simulation was run on a 
Linux system with a Slurm Workload Manager. The simulation can be replicated by running the shell script simula.sh available in each folder for the risk difference (prob) or for the restricted mean survival time (rmst). 
Some of the options for batching the jobs are native to the host system 
they were executed on and thus will error if executed on other 
systems. The basic idea is that the bash script requests an array task and submits a job
to each process in the array. The script code.r creates a set of tasks to be submitted to each array process, and contains a function funslave 
that splits the tasks across processes in the array according to the process identifier. The script simula.r is run on each process, it reads the process 
identifier using the system variable SLURM_ARRAY_TASK_ID, and then runs the tasks corresponding to that index. Finally, the script summarize.r 
can be used to compile the results and produce the LaTeX tables in the paper. 

# Questions

Please reach out to Iván Díaz by [email](ild2005@med.cornell.edu) if you have any questions about the `survtmlerct` package or the 
simulation. 
