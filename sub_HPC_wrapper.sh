#!/bin/bash

#$ -S /bin/bash          # run with this shell
#$ -N f50_tau20_fp    # this name shows in qstat
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
TAU="20"
TAU_NAME="20_3SPEC_f50_062017"
SPECIES="3"
#BNGfile="NANOG_WOLYNES_f_100.net"

BINS="250"
LOOPS="1000"
VORLOOPS="1"
REPS="200"
FILE_ROOT="NANOG"  #NAME OF THE FILE
BATCH="60"
RSTART="102"
OUTPUT="/dfs2/elread/rxn-share/adaptive_WE_HPC" #LOCATION WHERE YOU WANT TO SAVE THE OUTPUTS
EIGV="MISA_VVs.mat"

ME="${USER}"
TMP="/fast-scratch/${ME}/stage_f50_tau20_062017"  #TEMP LOCATION, MUST BE IN FAST-SCRATCH. MAKE SURE THE NAME DOESN'T OVERLAP WITH ANYTHING ELSE RUNNING
mkdir $TMP

INFILEMAT="trans_NANOG_tau20_3SPEC_f50_062017.mat" 
INFILEVOR="2"
FULLRSTART="0"

#tar -C $TMP -xf input/${INFILEVOR}.tar.gz
echo $INFILEMAT
echo $INFILEVOR
module load MATLAB/r2015a
cp /dfs2/elread/rxn-share/BioNetGen-2.2.6-stable/bin/run_network ${TMP}/
cp /dfs2/elread/rxn-share/adaptive_WE_HPC/sub_BNG_scripts_batch_051217.sh ${TMP}/


#mcc -m General_trans_062617_NANOG_HPC_fullspec_f50_restart_rw.m
./General_trans_062617_NANOG_HPC_fullspec_f50_restart_rw $VORLOOPS $LOOPS $TMP $BINS $REPS $FILE_ROOT $TAU $TAU_NAME $OUTPUT $INFILEMAT $INFILEVOR $SPECIES $BATCH $RSTART $FULLRSTART $EIGV


