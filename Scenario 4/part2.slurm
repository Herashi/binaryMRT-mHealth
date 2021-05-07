#!/bin/bash
#SBATCH --job-name=combine
#SBATCH --time=1:00:00
#SBATCH --mail-type=begin,end,fail
#SBATCH --account=wdem2
#SBATCH --partition=standard
#SBATCH --mem=5g
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/herashi/MRT/slurm_output/slurm-%j.out

module load R
R CMD BATCH --no-save --no-restore /home/herashi/MRT/indirect_effect/100_25/part2.R


#!/bin/bash
#Put your slurm configurations here
#Slurm configurations begin with #SBATCH.
#To comment out a slurm configuration use two # symbols
##SBATCH This line is commented out

#Get email when the job reaches these states
#SBATCH --mail-type=BEGIN,END,FAIL

#Set your memory request
#SBATCH --mem=1G

#Set your job runtime
#Acceptable time formats include
# "minutes", "minutes:seconds", "hours:minutes:seconds"
# "days-hours", "days-hours:minutes", "days-hours:minutes:seconds"
#SBATCH --time=10

#Set your number of processing cores.
#Unless you have done work to parallelize your code, this should remain 1.
#SBATCH --cpus-per-task=1

#SBATCH --output=/home/herashi/binary_mrt/slurm_output/slurm-%j.out

#Put your code after this line
#If you are running and R job, write something like this:
R CMD BATCH --no-save --no-restore /home/herashi/binary_mrt/Scenario4/part2.R