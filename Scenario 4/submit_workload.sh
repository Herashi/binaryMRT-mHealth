#!/bin/bash
part1tmax30=$(sbatch --job-name=tmax30 --export=ALL,tmax=30 part1.slurm | awk '{print $NF}')
# part1tmax50=$(sbatch --job-name=tmax50 --export=ALL,tmax=50 part1.slurm | awk '{print $NF}')
sbatch --dependency=afterok:$part1tmax30 part2.slurm

# :$part1tmax50
