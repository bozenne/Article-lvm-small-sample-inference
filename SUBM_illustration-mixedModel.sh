#!/bin/bash

#$ -N illustration-mixedModel  # Job name
#$ -t 1:40     # Number of jobs
#$ -q all.q    # Queue. Use long.q for run time >8h and all.q otherwise
#$ -l h_vmem=1G # Memory limit, e.g. reserve 1 GB memory 
#$ -tc 128      # Max concurrent jobs
#$ -cwd         # Run in current directory
#$ -o output/illustration-mixedModel/   # Direct output to subdirectory
#$ -e output/illustration-mixedModel/   # Direct output to subdirectory

R CMD BATCH BATCH_illustration-mixedModel.R output/illustration-mixedModel/$JOB_NAME-I-$SGE_TASK_ID.Rout --no-restore --no-save

## go to directory    ## cd Cluster/LVMproject/article-smallSampleInference/
## clean outputs      ## rm -r ./output/illustration-mixedModel/*
## clean results      ## rm -r ./Results/illustration-mixedModel/*
## submission command ## qsub SUBM_illustration-mixedModel.sh

## submission output  ## Your job-array 11186.1-40:1 ("illustration-mixedModel") has been submitted
## submission time    ## 01/31/19 3:52 
#    user  system elapsed 
# 748.370   3.610 754.717

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
