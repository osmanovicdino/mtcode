#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=8:00:00,h_data=1024M
#$ -t 1-200:1

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
filename=~/Chemistry/Code/Basic/Scripts/paramsrandoms.dat
basedir="PhaseDiagramDesignAll"
if [ -e ${filename}   ]; then
   # use the unix command sed -n ${line_number}p to read by line
   den=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $1}'`
   i1=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $2}'` 
   i2=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $3}'` 
   i3=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $4}'`
   i4=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $5}'`
   i5=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $6}'`
   m1=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $7}'`
   m2=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $8}'`
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
   ae=-20.
   e=2.
   wi=0
   echo "did not read file correctly"
fi
dirwemake="den=${den}_i1=${i1}_i2=${i2}_i3=${i3}_i4=${i4}_i5=${i5}_m1=${m1}_m2=${m2}_ae=${ae}_ie=${ie}_num=${SGE_TASK_ID}"
mkdir /u/scratch/d/dinoo/${basedir}/${dirwemake}
cp ~/Chemistry/Code/mainlocalRandom.cpp /u/scratch/d/dinoo/${basedir}/${dirwemake}
g++ ~/Chemistry/Code/mainlocalRandom.cpp -o /u/scratch/d/dinoo/${basedir}/${dirwemake}/angron
cd /u/scratch/d/dinoo/${basedir}/${dirwemake}
./angron 1000000 $den $i1 $i2 $i3 $i4 $i5 $m1 $m2 $ae $ie >log
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####