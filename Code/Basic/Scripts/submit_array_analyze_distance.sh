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
#$ -pe shared 8
#$ -t 1-36:1

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
filename=~/Chemistry/Code/Basic/Scripts/paramsgrowth3analyze.dat
basedir="PhaseDiagramDesign8"
if [ -e ${filename}   ]; then
   # use the unix command sed -n ${line_number}p to read by line
   strdir=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $1}'`
   l=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $2}'`
   d=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $3}'`
   echo "read file correctly" 
else
   echo "did not read file correctly"
   exit 1
fi
dirwemake=$strdir
# cp ~/Chemistry/Code/mainanalyzedist.cpp /u/scratch/d/dinoo/${basedir}/${dirwemake}
# g++ ~/Chemistry/Code/mainanalyze.cpp -o /u/scratch/d/dinoo/${basedir}/${dirwemake}/angron_analyze
cp ~/Chemistry/Code/studygrowth /u/scratch/d/dinoo/${basedir}/studygrowth${SGE_TASK_ID
cd /u/scratch/d/dinoo/${basedir}/
./studygrowth $strdir $l $d >log${SGE_TASK_ID}
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####