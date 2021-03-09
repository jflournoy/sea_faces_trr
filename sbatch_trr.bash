#!/bin/bash
#SBATCH -J trr
#SBATCH --time=0-01:00:00
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=9G
#SBATCH -p ncf
#SBATCH --account=mclaughlin_lab
# Outputs ----------------------------------
#SBATCH -o log/%x-%A-%a.out

source /users/jflournoy/code/R_3.5.1_modules.bash

Rscript --vanilla bicc.R
