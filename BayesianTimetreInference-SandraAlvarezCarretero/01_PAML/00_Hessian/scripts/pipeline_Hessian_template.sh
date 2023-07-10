#!/bin/bash

#==========================================================================================#
# Contact Sandra Alvarez-Carretero for any doubts about this script: sandra.ac93@gmail.com #
#==========================================================================================#

# 1. Find global dirs for paths
pipeline_dir=$( pwd )
main_dir=$( echo $pipeline_dir | sed 's/\/main\/..*/\/main\//' )

# 2. Set a `for` loop to run the tasks with as
# many alignments as defined by the user
for TASK_ID in `seq 1 NUM`
do

# ------------------------------------- #
# Creating file structure to run BASEML #
# ------------------------------------- # 

# 3. Move to main wd
cd $main_dir/Hessian/$TASK_ID
home_dir=$( pwd ) 

# 4. Create specific log file
exec 3>&1> >(while read line; do echo "$line" >> $pipeline_dir/log.hessian.dir$TASK_ID.txt; done;) 2>&1
start=`date`
echo Job starts":" $start

# 5. Start analysis
echo The analyses will take place in directory $home_dir
printf "\n"
# Move to analysis dir
cd $home_dir
# Soft link the tmp* files here 
ln -s $home_dir/prepare_baseml/tmp0001.ctl .
ln -s $home_dir/prepare_baseml/tmp0001.trees .
ln -s $home_dir/prepare_baseml/tmp0001.txt .

# 6. Run BASEML
printf "\nRunning BASEML to calculate the Hessian and the gradient...\n"
baseml tmp0001.ctl

# 7. Close
printf "\n"
echo BASEML FINISHED"!"
printf "\n"
end=`date`
echo Job ends":" $end
done
