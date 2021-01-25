#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=2:00:00,h_data=8G
## Modify the parallel environment
## and the number of cores as needed:
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea
#$ -t 1-48:1

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh
## Edit the line below as needed:
module load gcc/4.9.3
module load mathematica/12.1


if [ -e  ~/Chemistry/Code/Basic/Scripts/params.dat ]; then
   # use the unix command sed -n ${line_number}p to read by line
   dir=`sed -n ${SGE_TASK_ID}p ~/Chemistry/Code/Basic/Scripts/directories.txt`
   echo "read file correctly" 
else
   echo "did not read file correctly"
fi
## substitute the command to run your code
## in the two lines below:
##echo '/usr/bin/time -v hostname'
##/usr/bin/time -v hostname
##mkdir ~/Chemistry/Results1
##cp ~/Chemistry/Code/main.cpp ~/Chemistry/Results1/
##g++ -fopenmp -std=c++11 ~/Chemistry/Code/main.cpp -o ~/Chemistry/Results1/angron
##cd ~/Chemistry/Results1/
##./angron
all_lines=`ls -d ${dir}/pos*.csv`
for item in $all_lines; 
do
	echo "starting file\n";
	echo "$item";
   echo "\n";
	./Code/Plotting/PlotFrame.wls "$item" ${dir}/col.csv
   echo "\n"
done

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####