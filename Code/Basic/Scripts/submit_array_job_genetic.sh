#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=23:59:59,h_data=512M
## Modify the parallel environment
#$ -t 1-100:1


# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh
## Edit the line below as needed:
module load gcc/11.3.0

## substitute the command to run your code
## in the two lines below:
##echo '/usr/bin/time -v hostname'
##/usr/bin/time -v hostname

dirwemake="den${SGE_TASK_ID}"
subdir="GeneticTry4"
mkdir /u/scratch/d/dinoo/${subdir}/${dirwemake}
cp ~/Chemistry/Code/mainNanotubeBox.cpp /u/scratch/d/dinoo/${subdir}/${dirwemake}
g++ ~/Chemistry/Code/mainNanotubeBox.cpp -std=c++17 -o /u/scratch/d/dinoo/${subdir}/${dirwemake}/angron
cd /u/scratch/d/dinoo/${subdir}/${dirwemake}

./angron >log
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####