#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=12:00:00,h_data=1G
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 8
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh
## Edit the line below as needed:
module load gcc/7.5.0

## substitute the command to run your code
## in the two lines below:
##echo '/usr/bin/time -v hostname'
##/usr/bin/time -v hostname
den=0.02
i1=8.0
i2=0.0
i3=0.0
i4=0.0
m1=0
m2=0
dirwemake="den=${den}_i1=${i1}_i2=${i2}_i3=${i3}_i4=${i4}_m1=${m1}_m2=${m2}"
mkdir /u/scratch/d/dinoo/Tests/${dirwemake}
cp ~/Chemistry/Code/main.cpp /u/scratch/d/dinoo/Tests/${dirwemake}
g++ -std=c++11 -fopenmp ~/Chemistry/Code/mainlocal.cpp -o /u/scratch/d/dinoo/Tests/${dirwemake}/angron
cd /u/scratch/d/dinoo/Tests/${dirwemake}
export OMP_NUM_THREADS=8
./angron 10000000 $den $i1 $i2 $i3 $i4 $m1 $m2 >log
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####