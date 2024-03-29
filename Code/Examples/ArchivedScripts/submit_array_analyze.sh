#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=00:15:00,h_data=256M
## Modify the parallel environment
## and the number of cores as needed:
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea
#$ -t 1-180:1

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
filename=~/Chemistry/Code/Basic/Scripts/paramsgrowth9analyze.dat
basedir="GrowthRun5"
if [ -e ${filename}   ]; then
   # use the unix command sed -n ${line_number}p to read by line
   den=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $1}'`
   d=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $2}'` 
   i=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $3}'` 
   an=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $4}'`
   ar=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $5}'`
   strdir=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $6}'`
   echo "read file correctly" 
else
   echo "did not read file correctly"
   exit 1
fi
dirwemake=$strdir
cp ~/Chemistry/Code/mainanalyze.cpp /u/scratch/d/dinoo/${basedir}/${dirwemake}
# g++ ~/Chemistry/Code/mainanalyze.cpp -o /u/scratch/d/dinoo/${basedir}/${dirwemake}/angron_analyze
cp ~/Chemistry/Code/angron_analyze /u/scratch/d/dinoo/${basedir}/${dirwemake}/angron_analyze
cd /u/scratch/d/dinoo/${basedir}/${dirwemake}
./angron_analyze $den $d $i $an $ar $strdir >log
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####