#!/bin/bash

#$ -N IV-non-normal  # Job name
#$ -t 1:40     # Number of jobs
#$ -q long.q    # Queue. Use long.q for run time >8h and all.q otherwise
#$ -l h_vmem=5G # Memory limit, e.g. reserve 1 GB memory 
#$ -tc 128      # Max concurrent jobs
#$ -cwd         # Run in current directory
#$ -o output/IV-non-normal/   # Direct output to subdirectory
#$ -e output/IV-non-normal/   # Direct output to subdirectory

R CMD BATCH BATCH_IV-non-normal.R output/IV-non-normal/$JOB_NAME-I-$SGE_TASK_ID.Rout --no-restore --no-save

## go to directory    ## cd Cluster/LVMproject/article-smallSampleInference/
## clean outputs      ## rm -r ./output/IV-non-normal/*
## clean results      ## rm -r ./Results/IV-non-normal/*
## submission command ## qsub SUBM_IV-non-normal.sh

## submission output  ## Your job-array 11574.1-40:1 ("IV-non-normal") has been submitted
## submission time    ## 03/07/19 12:32 
#     user   system  elapsed 
# 5273.805   46.836 4877.423 

## documentation      ## https://ifsv.sund.ku.dk/biostat/biostatwiki/index.php/IT:Cluster : biostat wiki about the cluster
                      ## http://gridscheduler.sourceforge.net/htmlman/manuals.html : grid engine manual 
                      ## http://bayes/ganglia                                      : current load and history of the cluster

## commands           ## qstat         : view jobs of the user
                      ## qstat -u *   : view jobs of all users (the first column shows the job id)
                      ## qstat -j 1034 : show details of a job (or job array) with job id 1034 type     
                      ## qdel 1034     : delete the job with job id 1034 from the queue type
                      ## finger login  : get the name corresponding to the login

## status             ## qw : quewing
                      ##  r : running
                      ## dr : dual state (r)unning and being (d)eleted
