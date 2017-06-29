#!/bin/bash

#$ -S /bin/bash          # run with this shell
#$ -N temp_name    # this name shows in qstat
#$ -q rxn               # run in this Q
#$ -j y
#$ -cwd            # run the job out of the current directory;
                   # (the one from which you ran the script)
module load gmp/5.1.3
module load mpc/1.0.1 mpfr/3.1.2
module load Cluster_Defaults
module load binutils/2.23.2
module load gdb/7.8
module load openmpi-1.8.3/gcc-4.8.2
module load gcc/4.8.2
#One of the above modules is necessary to keep BNG running after 100 iterations or so. 
#I don't know which module is necessary

TAU="5000"
#tau for both Voronoi and transition steps of the simulation
TAU_NAME="5000_062917"
#suffix of the output file
SPECIES="8"
#number of columns in the BNG file to read

BINS="250"
#number of voronoi regions
LOOPS="1000"
#number of WE loops for the transition matrix
VORLOOPS="1"
#number of WE loops for the voronoi regions
REPS="200"
FILE_ROOT="MISA"  #PREFIX for the output file
BATCH="60" #number of BNG iterations in each job of the job-array (DON'T CHANGE)
RSTART="1" #current iteration to restart from, set to 1 if new simulation
OUTPUT="/dfs2/elread/rxn-share/WE_HPC_BATCH_BNG_CODES" #LOCATION WHERE YOU WANT TO SAVE THE OUTPUTS

ME="${USER}"
TMP="/fast-scratch/${ME}/stage_name"  #TEMP LOCATION, MUST BE IN FAST-SCRATCH. MAKE SURE THE NAME DOESN'T OVERLAP WITH ANYTHING ELSE RUNNING
mkdir $TMP

INFILEMAT="trans_NANOG_tau20_3SPEC_f50_062017.mat" #INPUT file for the voronoi, replicas and replica weights. Not necessary for initialization
INFILEVOR="2" #current replica data to restart from
FULLRSTART="0" #FLAG FOR restarting transition matrix (1 for full restart, 0 for continuation)


module load MATLAB/r2015a #loading matlab
cp ${OUTPUT}/run_network ${TMP}/
cp ${OUTPUT}/sub_BNG_batch_template.sh ${TMP}/
cp ${OUTPUT}/sub_BNG_MISAEx_init.sh ${TMP}/

mcc -m General_WE_trans_062617_MISA_init.m
./General_WE_trans_062617_MISA_init $VORLOOPS $LOOPS $TMP $BINS $REPS $FILE_ROOT $TAU $TAU_NAME $OUTPUT $SPECIES $BATCH $RSTART

#mcc -m General_WE_trans_062617_restart.m #for restarting
#./General_WE_trans_062617_restart $VORLOOPS $LOOPS $TMP $BINS $REPS $FILE_ROOT $TAU $TAU_NAME $OUTPUT $INFILEMAT $INFILEVOR $SPECIES $BATCH $RSTART $FULLRSTART $EIGV


