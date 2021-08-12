#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=12:00:00,h_data=128M
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 12
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea
#$ -t 1-15:1

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh
## Edit the line below as needed:
module load gcc/9.3.0

## substitute the command to run your code
## in the two lines below:
##echo '/usr/bin/time -v hostname'
##/usr/bin/time -v hostname
filename=~/Chemistry/Code/Basic/Scripts/params17.dat
basedir="PhaseDiagramDesign12"
if [ -e ${filename}   ]; then
   # use the unix command sed -n ${line_number}p to read by line
   den=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $1}'`
   i1=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $2}'` 
   i2=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $3}'` 
   i3=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $4}'`
   i4=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $5}'`
   m1=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $6}'`
   m2=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $7}'`
   rate=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $8}'`
   ae=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $9}'`
   ie=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $10}'`
   echo "read file correctly" 
else
   den=0.01
   i1=0
   i2=0
   i3=0
   i4=0
   m1=1000
   m2=1000
   rate=0.01
   ae=-20.
   e=2.
   echo "did not read file correctly"
fi
dirwemake="den=${den}_i1=${i1}_i2=${i2}_i3=${i3}_i4=${i4}_m1=${m1}_m2=${m2}_rate=${rate}_ae=${ae}_ie=${ie}"
mkdir /u/scratch/d/dinoo/${basedir}/${dirwemake}
cp ~/Chemistry/Code/mainlocal.cpp /u/scratch/d/dinoo/${basedir}/${dirwemake}
g++ -fopenmp ~/Chemistry/Code/mainlocal.cpp -o /u/scratch/d/dinoo/${basedir}/${dirwemake}/angron
cd /u/scratch/d/dinoo/${basedir}/${dirwemake}
export OMP_NUM_THREADS=12
./angron 10000000 $den $i1 $i2 $i3 $i4 $m1 $m2 $rate $ae $ie >log
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####