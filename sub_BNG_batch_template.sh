#!/bin/bash      
#$ -ckpt restart
#$ -S /bin/bash          # run with this shell  
#$ -N temp    # this name shows in qstat   
#$ -q rxn,pub64,free*               # run in this Q
#$ -cwd
#$ -j y   
# RESTART_EMAIL_SUMMARY
#$ -l h=!compute-1-10&!compute-11-10&!compute-8-15      #current buggy compute nodes         
module load gmp/5.1.3   
module load mpc/1.0.1 mpfr/3.1.2   
module load Cluster_Defaults   
module load binutils/2.23.2   
module load gdb/7.8   
module load openmpi-1.8.3/gcc-4.8.2   
module load gcc/4.8.2   #ONE OF THESE MODULES IS NECESSARY FOR BNG. NOT SURE WHICH
OUT_NUM="0" #CURRENT 
TMP="/scratch/tsem1/TEMP_F10_NEW"
TSTEP="5.000000" 
NUMSIMSTEPS="2"   
PREV_ID="0" 
NEW_ID="1" 
START="1"
STOP="30"







#$ -t 1-1000 
> ${TMP}/tempout_${SGE_TASK_ID}.txt
for i in $(seq ${START} ${STOP});
  do
	CURRLINE="$(($((${STOP}*$((${SGE_TASK_ID}-1))))+${i}))"
	FILENAME="$(sed "${CURRLINE}q;d" curr_rep_ind.txt)"   
	PREV_FILE="${TMP}/t${PREV_ID}/newsim_${FILENAME}_end.net"   
	output5="${TMP}/t${NEW_ID}/newsim_${CURRLINE}"    
	${TMP}/run_network -o ${output5} -e -p ssa -h $RANDOM --cdat 0 -g $PREV_FILE $PREV_FILE $TSTEP $NUMSIMSTEPS &
	wait     
	TEMPout=""$(tail ${output5}.gdat -n 1)" "$(sed "${CURRLINE}q;d" ${TMP}/curr_rep_data.txt)""     
	echo $TEMPout >> ${TMP}/tempout_${SGE_TASK_ID}.txt
   done  
#cat ${TMP}/tempout_${SGE_TASK_ID}.txt >> ${TMP}/curr_out_rep_${OUT_NUM}.txt
echo ${SGE_TASK_ID} >> ${TMP}/finished.txt 
   
   
   
   
  
 
