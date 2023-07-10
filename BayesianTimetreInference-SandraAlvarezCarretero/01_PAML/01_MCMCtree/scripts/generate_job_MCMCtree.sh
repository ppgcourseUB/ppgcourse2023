#!/bin/bash

# Get args
dir=$1     # Alignment #1, #2, #3... 1 == conc | 2 == part
clock=$2   # `GBM` or `ILN`
ndat=$3    # 1, 2, 3... As many blocks as partitions in the alignment
pipeloc=$4 # Path to MCMCtree pipeline dir
runmcmc=$5 # Command to execute MCMCtree
nchains=$6 # Number of MCMCs to run
caltype=$7 # Type of calibration: `seccal` or `fosscal`

# Replace type of alignment
if [[ $dir -eq 1 ]]
then
partaln=$( echo "conc" )
else
partaln=$( echo "part" )
fi

# Replace vars in template bash script for job array
cp pipeline_MCMCtree_template.sh $pipeloc/$partaln/$caltype/$clock/pipeline_$clock"_"$caltype"_"$partaln".sh"
sed -i 's/PHY/'${partaln}'/g' $pipeloc/$partaln/$caltype/$clock/pipeline_$clock"_"$caltype"_"$partaln".sh"
sed -i 's/CAL/'${caltype}'/g' $pipeloc/$partaln/$caltype/$clock/pipeline_$clock"_"$caltype"_"$partaln".sh"
sed -i 's/DIR/'${dir}'/g' $pipeloc/$partaln/$caltype/$clock/pipeline_$clock"_"$caltype"_"$partaln".sh"
sed -i 's/CLK/'${clock}'/g' $pipeloc/$partaln/$caltype/$clock/pipeline_$clock"_"$caltype"_"$partaln".sh"
sed -i 's/NUMPARTS/'${ndat}'/' $pipeloc/$partaln/$caltype/$clock/pipeline_$clock"_"$caltype"_"$partaln".sh"
sed -i 's/CMDRUN/'${runmcmc}'/' $pipeloc/$partaln/$caltype/$clock/pipeline_$clock"_"$caltype"_"$partaln".sh"

if [[ $nchains -eq 1 ]]
then 
sed -i 's/for\ TASK\_ID\ in..*//' $pipeloc/$partaln/$caltype/$clock/pipeline_$clock"_"$caltype"_"$partaln".sh"
sed -i 's/^done//' $pipeloc/$partaln/$caltype/$clock/pipeline_$clock"_"$caltype"_"$partaln".sh"
sed -i 's/^do//' $pipeloc/$partaln/$caltype/$clock/pipeline_$clock"_"$caltype"_"$partaln".sh"
else 
sed -i 's/NUM/'${nchains}'/' $pipeloc/$partaln/$caltype/$clock/pipeline_$clock"_"$caltype"_"$partaln".sh"
fi


