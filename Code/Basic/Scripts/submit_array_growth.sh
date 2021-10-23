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
#$ -pe shared 8
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea
#$ -t 1-60:1

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
filename=~/Chemistry/Code/Basic/Scripts/paramsgrowth3.dat
basedir="GrowthRun3"
if [ -e ${filename}   ]; then
   # use the unix command sed -n ${line_number}p to read by line
   den=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $1}'`
   i1=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $2}'` 
   i2=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $3}'` 
   i3=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $4}'`
   i4=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $5}'`
   echo "read file correctly" 
else
   den=0.01
   i1=1.0
   i2=12.0
   i3=0.927
   i4=4
   echo "did not read file correctly"
fi
dirwemake="den=${den}_d=${i1}_e=${i2}_a=${i3}_arms=${i4}"
mkdir /u/scratch/d/dinoo/${basedir}/${dirwemake}
cp ~/Chemistry/Code/main.cpp /u/scratch/d/dinoo/${basedir}/${dirwemake}
g++ -fopenmp ~/Chemistry/Code/main.cpp -o /u/scratch/d/dinoo/${basedir}/${dirwemake}/angron
cd /u/scratch/d/dinoo/${basedir}/${dirwemake}
export OMP_NUM_THREADS=8
./angron 10000000 $den $i1 $i2 $i3 $i4 >log
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####