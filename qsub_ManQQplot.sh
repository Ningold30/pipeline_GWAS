#!/bin/bash
#PBS -N Manhattan_plot
#PBS -l mem=4gb,walltime=04:00:00,ncpus=1   
#PBS -o Manhattan_plot.o
#PBS -e Manhattan_plot.e

module load R/4.3.1

Rscript /mnt/backedup/home/nathanI/Scripts/plotting/Man_QQplot.R ${file1} ${dirR} ${title}
