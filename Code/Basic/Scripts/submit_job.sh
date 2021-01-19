#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=1:00:00,h_data=1G
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
module load gcc/4.9.3

## substitute the command to run your code
## in the two lines below:
##echo '/usr/bin/time -v hostname'
##/usr/bin/time -v hostname
mkdir ~/Chemistry/Results1
cp ~/Chemistry/Code/main.cpp ~/Chemistry/Results1/
g++ -fopenmp -std=c++11 ~/Chemistry/Code/main.cpp -o ~/Chemistry/Results1/angron
cd ~/Chemistry/Results1/
./angron
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####