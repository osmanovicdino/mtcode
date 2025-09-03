#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=12:00:00,h_data=512M
## Modify the parallel environment
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea
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
mkdir /u/scratch/d/dinoo/GeneticTry4/${dirwemake}
cp ~/Chemistry/Code/main.cpp /u/scratch/d/dinoo/GeneticTry4/${dirwemake}
g++ ~/Chemistry/Code/mainNanotubeBox.cpp -o /u/scratch/d/dinoo/GeneticTry4/${dirwemake}/angron
cd /u/scratch/d/dinoo/GeneticTry4/${dirwemake}

./angron >log
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####