#!/bin/bash
###############################################################################
##                                                                           ##
##                           Campus Cluster                                  ##
##                    Sample OpenMP Job Batch Script                         ##
##                                                                           ##
## SLURM Options (To view, run the following command below)                  ##
##                                                                           ##
##     man sbatch                                                            ##
##                                                                           ##
###############################################################################
#                                                                            
#SBATCH --time=32:00:00                    # Job run time (hh:mm:ss)
#SBATCH --nodes=1                          # Number of nodes
#SBATCH --ntasks-per-node=20               # Number of task (cores/ppn) per node
#SBATCH --job-name=cohort200fin                # Name of batch job
#SBATCH --partition=ichass                 # Partition (queue)           
#SBATCH --output=openmp.o%j                # Name of batch job output file
##SBATCH --error=openmp.e%j                # Name of batch job error file
##SBATCH --mail-user=tunder@illinois.edu   # Send email notifications
##SBATCH --mail-type=BEGIN,END             # Type of email notifications to send
#                                                                            
###############################################################################

# Change to the directory from which the batch job was submitted
# Note: SLURM defaults to running jobs in the directory where
# they are submitted, no need for cd'ing to $SLURM_SUBMIT_DIR

#cd ${SLURM_SUBMIT_DIR}

# Add scripts to run here
bin/mallet train-topics --input ../code/cohort/data/cohort4.mallet --num-threads 20 --num-topics 200 --num-iterations 800 --num-top-words 200 --random-seed 99 --optimize-interval 20 --optimize-burn-in 100 --output-state ../code/cohort/finalmodel/k200state.gz --output-doc-topics ../code/cohort/finalmodel/k200doctopics.txt --output-topic-keys ../code/cohort/finalmodel/k200keys.txt --show-topics-interval 200 --diagnostics-file ../code/cohort/finalmodel/k200diagnostics.xml
