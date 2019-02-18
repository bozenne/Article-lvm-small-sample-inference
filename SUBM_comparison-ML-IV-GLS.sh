#!/bin/bash

#$ -N comparison-ML-IV-GLS  # Job name
#$ -t 1:40     # Number of jobs
#$ -q all.q    # Queue. Use long.q for run time >8h and all.q otherwise
#$ -l h_vmem=5G # Memory limit, e.g. reserve 1 GB memory 
#$ -tc 128      # Max concurrent jobs
#$ -cwd         # Run in current directory
#$ -o output/comparison-ML-IV-GLS/   # Direct output to subdirectory
#$ -e output/comparison-ML-IV-GLS/   # Direct output to subdirectory

R CMD BATCH BATCH_comparison-ML-IV-GLS.R output/comparison-ML-IV-GLS/$JOB_NAME-I-$SGE_TASK_ID.Rout --no-restore --no-save

## go to directory    ## cd Cluster/LVMproject/article-smallSampleInference/
## clean outputs      ## rm -r ./output/comparison-ML-IV-GLS/*
## clean results      ## rm -r ./Results/comparison-ML-IV-GLS/*
## submission command ## qsub SUBM_comparison-ML-IV-GLS.sh

## submission output  ## Your job-array 11388.1-40:1 ("comparison-ML-IV-GLS") has been submitted
## submission time    ## 02/08/19 10:38 
#     user   system  elapsed 
# 9554.332   96.733 8283.346 

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