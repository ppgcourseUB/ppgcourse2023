#### RUN DIFFERENT CONDITIONS. OBTAIN SFS IN SYN AND NONSYN. CALCULATE MKTasymptotic. ####
#### GENERAK CONDITIONS:
#### Two populations: ONE target plus ONE outgroup.
#### 1.The initial population run for 5Ne generations to achieve some equilibrium mutation-selection-drift
#### 2.Split target and outgroup for 10Ne generations
#### 3.Possible change in number of individuals in target population after 5Ne+10Ne
#### 4.Possible selective sweeps (n simulataneous sweeps) 0.5Ne before end simulation

#header for run in slurm
echo \#!/bin/bash > ./run_slim_conditions.sh
echo \# >> ./run_slim_conditions.sh
echo \#SBATCH --job-name=run_slim_sim_01 >> ./run_slim_conditions.sh
echo \#SBATCH -o %j.out >> ./run_slim_conditions.sh
echo \#SBATCH -e %j.err >> ./run_slim_conditions.sh
echo \#SBATCH --ntasks=9 >> ./run_slim_conditions.sh
echo \#SBATCH --mem=6GB >> ./run_slim_conditions.sh
echo \#SBATCH --partition=normal >> ./run_slim_conditions.sh
echo \# >> ./run_slim_conditions.sh
echo module load SLiM >> ./run_slim_conditions.sh
echo >> ./run_slim_conditions.sh

#fixed paraneters
Ne=750; L=500000; ngenes=100;
mut_rate=1e-6;
ind_sample_size=25; out_sample_size=1;

# CONDITION 0:
#Neutral. No change Ne.
FILEOUT="'./00_slim_SFS_SNM.txt'"
Neb=750; nsweeps=0;
rec_rate=1e-4;
rate_ben=0; s_backg_ben=0;
rate_del=0; s_backg_del=0;
rate_neut=1;

echo srun --ntasks 1 --exclusive --mem-per-cpu=1GB slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim\& >> ./run_slim_conditions.sh

# CONDITION 1:
#Strong BACKGROUND SELECTION. No beneficial selection. No change Ne. No sweep
FILEOUT="'./01_slim_SFS_BGS.txt'"
Neb=750; nsweeps=0;
rec_rate=1e-4
rate_ben=0; s_backg_ben=0;
rate_del=9; s_backg_del=-0.005;
rate_neut=1;

echo srun --ntasks 1 --exclusive --mem-per-cpu=1GB slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim\& >> ./run_slim_conditions.sh

# CONDITION 2:
#No background selection. BENEFICIAL SELECTION. No change Ne. No sweep
FILEOUT="'./02_slim_SFS_PS.txt'"
Neb=750; nsweeps=0;
rec_rate=1e-4
rate_ben=0.2; s_backg_ben=0.002;
rate_del=0; s_backg_del=0;
rate_neut=1;

echo srun --ntasks 1 --exclusive --mem-per-cpu=1GB slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim\& >> ./run_slim_conditions.sh

# CONDITION 3:
#Strong BACKGROUND SELECTION. MIDDLE PROPORTION BENEFICIAL SELECTION. No change Ne. No sweep.
FILEOUT="'./03_slim_SFS_BGS_PS.txt'"
Neb=750; nsweeps=0;
rec_rate=1e-4
rate_ben=0.2; s_backg_ben=0.005;
rate_del=8; s_backg_del=-0.002;
rate_neut=1;

echo srun --ntasks 1 --exclusive --mem-per-cpu=1GB slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim\& >> ./run_slim_conditions.sh

echo wait >> ./run_slim_conditions.sh
