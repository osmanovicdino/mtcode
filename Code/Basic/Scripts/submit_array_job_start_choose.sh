#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=00:59:59,h_data=128M
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 6
#$ -t 1-25:1

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
#filename=~/Chemistry/Code/Basic/Scripts/paramsSA6.dat
#basedir="SelfAssembly6"

filename=~/Chemistry/Code/Basic/Scripts/parameter_settings20.dat
if [ -e ${filename}   ]; then
   # use the unix command sed -n ${line_number}p to read by line
   wt=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $1}'`
   wt2=`sed -n ${SGE_TASK_ID}p ${filename} | awk '{print $2}'`
   echo "read file correctly" 
else
   wt=~
   echo "did not read file correctly"
fi
basedir="SelfAssembly32"
dirwemake="try=${SGE_TASK_ID}"
mkdir /u/scratch/d/dinoo/${basedir}/${dirwemake}
cp ~/Chemistry/Code/mainNanotubeElasticShell4.cpp /u/scratch/d/dinoo/${basedir}/${dirwemake}
cp ~/Chemistry/Code/IsocohedronI.csv /u/scratch/d/dinoo/${basedir}/${dirwemake}
cp ~/Chemistry/Code/${wt2} /u/scratch/d/dinoo/${basedir}/${dirwemake}/IsocohedronP.csv
cp ~/Chemistry/Code/Basic/InitialConditions/${wt} /u/scratch/d/dinoo/${basedir}/${dirwemake}/param.csv
g++ -fopenmp -std=c++17 ~/Chemistry/Code/mainNanotubeElasticShell4.cpp -o /u/scratch/d/dinoo/${basedir}/${dirwemake}/angron
cd /u/scratch/d/dinoo/${basedir}/${dirwemake}
export OMP_NUM_THREADS=6
./angron 'param.csv' >log
# cp ~/Chemistry/Code/mainNanotubeElasticShell.cpp /u/scratch/d/dinoo/${basedir}/${dirwemake}
# cp ~/Chemistry/Code/params.csv /u/scratch/d/dinoo/${basedir}/${dirwemake}
# cp ~/Chemistry/Code/IsocohedronP.csv /u/scratch/d/dinoo/${basedir}/${dirwemake}
# g++ -fopenmp -std=c++17 ~/Chemistry/Code/mainNanotubeElasticShell.cpp -o /u/scratch/d/dinoo/${basedir}/${dirwemake}/angron
# cd /u/scratch/d/dinoo/${basedir}/${dirwemake}


# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####