#!/bin/bash
#SBATCH -J icc
#SBATCH --time=0-01:00:00
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -p ncf
#SBATCH --account=mclaughlin_lab
# Outputs ----------------------------------
#SBATCH -o log/%x-%A-%a.out

source /users/jflournoy/code/R_3.5.1_modules.bash

Rscript --vanilla oldicc.R