#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=12:00:00,h_data=16G
## Modify the parallel environment
## and the number of cores as needed:
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea
#$ -t 1-50:1

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh
## Edit the line below as needed:
module load gcc/7.5.0
module load mathematica/12.1
module load ffmpeg

if [ -e  ~/Chemistry/Code/Basic/Scripts/directories.txt ]; then
   # use the unix command sed -n ${line_number}p to read by line
   dir=`sed -n ${SGE_TASK_ID}p ~/Chemistry/Code/Basic/Scripts/directories.txt`
   echo "read file correctly" 
else
   dir="/u/scratch/d/dino/"
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
cd $dir
all_lines=`ls -d ./pos*.csv`
for item in $all_lines; 
do
   echo "starting file\n";
   echo "$item";
   echo "\n";
        wolframscript -file /u/home/d/dinoo/Chemistry/Code/Plotting/PlotFrame2.wl "$item";
   echo "\n";
done

ffmpeg -i pos_beta\=1_i\=%05d.jpg -vcodec libx264 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -pix_fmt yuv420p test.mp4

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####