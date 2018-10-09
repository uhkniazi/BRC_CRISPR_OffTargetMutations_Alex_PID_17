# Autogenerated script from bwa_array_job.R
# date Tue Oct  9 12:42:42 2018
# make sure directory paths exist before running script
#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 1
#$ -cwd
#$ -N bwa_mem-array
#$ -j y
#$ -l h_vmem=19G
#$ -t 1-2



module load bioinformatics/bwa/0.7.15



# Parse parameter file to get variables.
number=$SGE_TASK_ID
paramfile=bwa_param.txt
 
inr1=`sed -n ${number}p $paramfile | awk '{print $1}'`
inr2=`sed -n ${number}p $paramfile | awk '{print $2}'`
outsam=`sed -n ${number}p $paramfile | awk '{print $3}'`

# 9. Run the program.
bwa mem /users/k1625253/brc_scratch/Data/MetaData/GenomeIndex/hg38_bwa/hg38.fa $inr1 $inr2 > $outsam



