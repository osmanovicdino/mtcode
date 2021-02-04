#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=01:00:00,h_data=1G
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 36
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
den=0.01
i1=10.
i2=10.
i3=10.
dirwemake="den=${den}_i1=${i1}_i2=${i2}_i3=${i3}"
mkdir /u/scratch/d/dinoo/PhaseDiagramBivalent/${dirwemake}
cp ~/Chemistry/Code/main.cpp /u/scratch/d/dinoo/PhaseDiagramBivalent/${dirwemake}
g++ -fopenmp -std=c++11 -pg -no-pie ~/Chemistry/Code/main.cpp -o /u/scratch/d/dinoo/PhaseDiagramBivalent/${dirwemake}/angron
cd /u/scratch/d/dinoo/PhaseDiagramBivalent/${dirwemake}
./angron $den $i1 $i2 $i3 >log
gprof angron gmon.out > analysis2.txt
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####