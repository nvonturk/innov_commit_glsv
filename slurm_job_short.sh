#!/bin/bash -login

##### change these accordingly like you would your job #########

#SBATCH -J innov_commit # Name of job
#SBATCH -n 30 #Request tasks (cores)
#SBATCH -t 0-12:00 #Request runtime of days-hours:minutes
#SBATCH -C centos7 #Request only Centos7 nodes
#SBATCH -p sched_mit_hill #Run on sched_engaging_default partition
#SBATCH --mem-per-cpu=5GB #Request 4G of memory per CPU
#SBATCH -o Sbatch_output/output_%j.txt #redirect output to output_JOBID.txt
#SBATCH -e Sbatch_error/error_%j.txt #redirect errors to error_JOBID.txt
#SBATCH --mail-type=BEGIN,END,FAIL #Mail when job starts and ends
#SBATCH --mail-user=nvonturk@mit.edu #email recipient

./slurm_job_matlab.sh
