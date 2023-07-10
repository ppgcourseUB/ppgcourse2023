#!/bin/bash

#==========================================================================================#
# Contact Sandra Alvarez-Carretero for any doubts about this script: sandra.ac93@gmail.com #
#==========================================================================================#

# 1. First arg is dataset
dir=$( echo DIR )       # 1, 2, 3... According to the dataset being analysed

# 2. Find global dirs for paths
pipeline_dir=$( pwd )
main_dir=$( echo $pipeline_dir | sed 's/\/main\/..*/\/main\//' )

# 3. Set a `for` loop to run the tasks with as
#    many alignments as defined by the user
for TASK_ID in `seq 1 NUM`
do             

# 3.1. Creating file structure to run MCMCtree 
cd $main_dir/Hessian/$dir
baseml_dir=$( pwd )
name_inBV=$( echo in.BV )
path_inBV=$( echo $baseml_dir/$name_inBV )
cd $main_dir/MCMCtree/$TASK_ID/PHY/CAL/CLK
home_dir=$( pwd )
cd $main_dir/control_files
ctl_dir=$( pwd )
cd $main_dir/alignments/$dir 
aln_dir=$( pwd )
name_aln=`ls *phy`
path_aln=$( echo $aln_dir/$name_aln | sed 's/\_/\\\_/g' | sed 's/\//\\\//g' | sed 's/\./\\\./g' )
cd $main_dir/trees/calibrated
tree_dir=$( pwd )
name_tree=`ls *_CAL.tree`
path_tree=$( echo $tree_dir/$name_tree | sed 's/\_/\\\_/g' | sed 's/\//\\\//g' | sed 's/\./\\\./g' )

# 3.2. Create specific log file
exec 3>&1> >(while read line; do echo "$line" >> $pipeline_dir/log.MCMCtree.dir$dir"_r"$TASK_ID".txt"; done;) 2>&1
start=`date`
echo Job starts":" $start

# 3.3. Start analysis
echo The analyses will take place in directory $home_dir
printf "\n"
# Move to analysis dir
cd $home_dir
# Copy control file
cp $ctl_dir/tmp_CAL.ctl $home_dir/mcmctree_CAL_PHY_$dir"_r"$TASK_ID".ctl"
sed -i 's/ALN/'${path_aln}'/' $home_dir/mcmctree_CAL_PHY_$dir"_r"$TASK_ID".ctl"
sed -i 's/TREE/'${path_tree}'/' $home_dir/mcmctree_CAL_PHY_$dir"_r"$TASK_ID".ctl"
sed -i 's/usedata..*/usedata\ \=\ 2\ \.\/in\.BV\ 1/' $home_dir/mcmctree_CAL_PHY_$dir"_r"$TASK_ID".ctl"
sed -i 's/\ ndata..*/\ ndata\ \=\ NUMPARTS/' $home_dir/mcmctree_CAL_PHY_$dir"_r"$TASK_ID".ctl"
if [[ CLK =~ "GBM" ]]
then 
sed -i 's/clock..*/clock\ \=\ 3/' $home_dir/mcmctree_CAL_PHY_$dir"_r"$TASK_ID".ctl"
elif [[ CLK =~ "ILN" ]]
then 
sed -i 's/clock..*/clock\ \=\ 2/' $home_dir/mcmctree_CAL_PHY_$dir"_r"$TASK_ID".ctl"
fi 

# Soft link the in.BV file here
ln -s $path_inBV $home_dir/in.BV

# 3.4. Run MCMCtree
printf "\nRunning MCMCtree for divergence times estimation ... ...\n"
cd $home_dir
CMDRUN *.ctl

# 3.5. Close
printf "\n"
echo MCMCtree FINISHED"!"
printf "\n"
end=`date`
echo Job ends":" $end

# 3.6. Return to pipedir
cd $pipeline_dir

done
