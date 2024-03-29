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
#$ -t 1-20:1

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh
## Edit the line below as needed:
module load gcc/10.2.0

## substitute the command to run your code
## in the two lines below:
##echo '/usr/bin/time -v hostname'
##/usr/bin/time -v hostname
filename=~/Chemistry/Code/Basic/Scripts/params22.dat
basedir="GasLiquid3"
if [ -e ${filename}   ]; then
   # use the unix command sed -n ${line_number}p to read by line
   m1=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $1}'`
   m2=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $2}'`
   i1=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $3}'`
   echo "read file correctly" 
else
   m1=500;
   m2=500;
   i1=12.;
   echo "did not read file correctly"
fi
dirwemake="m1=${m1}_m2={$m2}_i1={$i1}"
mkdir /u/scratch/d/dinoo/${basedir}/${dirwemake}
cp ~/Chemistry/Code/mainlocal2.cpp /u/scratch/d/dinoo/${basedir}/${dirwemake}
g++ -fopenmp -std=c++17 ~/Chemistry/Code/mainlocal2.cpp -o /u/scratch/d/dinoo/${basedir}/${dirwemake}/angron
cd /u/scratch/d/dinoo/${basedir}/${dirwemake}
export OMP_NUM_THREADS=12
./angron 20000000 0.005 12. $i1 $m1 $m2 >log
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####