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
#$ -pe shared 12
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
module load gcc/9.3.0

## substitute the command to run your code
## in the two lines below:
##echo '/usr/bin/time -v hostname'
##/usr/bin/time -v hostname
den=0.01
i1=15.
i2=20.
i3=18.
dirwemake="den=${den}_i1=${i1}_i2=${i2}_i3=${i3}"
mkdir /u/scratch/d/dinoo/PhaseDiagramBivalent/${dirwemake}
cp ~/Chemistry/Code/main.cpp /u/scratch/d/dinoo/PhaseDiagramBivalent/${dirwemake}
g++ ~/Chemistry/Code/main.cpp -o /u/scratch/d/dinoo/PhaseDiagramBivalent/${dirwemake}/angron
cd /u/scratch/d/dinoo/PhaseDiagramBivalent/${dirwemake}
export OMP_NUM_THREADS=36
./angron 10000000 $den $i1 $i2 $i3 >log
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####