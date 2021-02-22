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
#$ -t 1-45:1

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
if [ -e  ~/Chemistry/Code/Basic/Scripts/params6.dat ]; then
   # use the unix command sed -n ${line_number}p to read by line
   den=`sed -n ${SGE_TASK_ID}p ~/Chemistry/Code/Basic/Scripts/params6.dat | awk '{print $1}'`
   i1=`sed -n ${SGE_TASK_ID}p ~/Chemistry/Code/Basic/Scripts/params6.dat | awk '{print $2}'` 
   i2=`sed -n ${SGE_TASK_ID}p ~/Chemistry/Code/Basic/Scripts/params6.dat | awk '{print $3}'` 
   i3=`sed -n ${SGE_TASK_ID}p ~/Chemistry/Code/Basic/Scripts/params6.dat | awk '{print $4}'`
   echo "read file correctly" 
else
   den=0.01
   i1=0
   i2=0
   i3=0
   echo "did not read file correctly"
fi
dirwemake="den=${den}_i1=${i1}_i2=${i2}_i3=${i3}"
mkdir /u/scratch/d/dinoo/PhaseDiagramBivalent/${dirwemake}
cp ~/Chemistry/Code/main.cpp /u/scratch/d/dinoo/PhaseDiagramBivalent/${dirwemake}
g++ -std=c++11 -fopenmp ~/Chemistry/Code/main.cpp -o /u/scratch/d/dinoo/PhaseDiagramBivalent/${dirwemake}/angron
cd /u/scratch/d/dinoo/PhaseDiagramBivalent/${dirwemake}
export OMP_NUM_THREADS=8
./angron 10000000 $den $i1 $i2 $i3 >log
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####