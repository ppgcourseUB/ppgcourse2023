#!/bin/bash
#
#SBATCH --job-name=run_slim_sim_01
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --ntasks=9
#SBATCH --mem=12GB
#SBATCH --partition=normal
#
module load SLiM

srun --ntasks 1 --exclusive --mem-per-cpu=1GB slim -t -m -d "Ne=1000" -d "L=500000" -d "Neb=1000" -d "mut_rate=1e-6" -d "rec_rate=1e-4" -d "ngenes=100" -d "rate_ben=0" -d "rate_del=0" -d "s_backg_ben=0" -d "s_backg_del=0" -d "nsweeps=0" -d "freq_sel_init=0.05" -d "freq_sel_end=0.95" -d "s_beneficial=0.0" -d "ind_sample_size=25" -d "out_sample_size=1" -d "file_output1='./00_slim_SFS_SNM.txt'" ./slim_template.slim&
srun --ntasks 1 --exclusive --mem-per-cpu=1GB slim -t -m -d "Ne=1000" -d "L=500000" -d "Neb=1000" -d "mut_rate=1e-6" -d "rec_rate=1e-4" -d "ngenes=100" -d "rate_ben=0" -d "rate_del=9" -d "s_backg_ben=0" -d "s_backg_del=-0.005" -d "nsweeps=0" -d "freq_sel_init=0.05" -d "freq_sel_end=0.95" -d "s_beneficial=0.0" -d "ind_sample_size=25" -d "out_sample_size=1" -d "file_output1='./01_slim_SFS_BGS.txt'" ./slim_template.slim&
srun --ntasks 1 --exclusive --mem-per-cpu=1GB slim -t -m -d "Ne=1000" -d "L=500000" -d "Neb=1000" -d "mut_rate=1e-6" -d "rec_rate=1e-4" -d "ngenes=100" -d "rate_ben=0.20" -d "rate_del=0" -d "s_backg_ben=0.002" -d "s_backg_del=0" -d "nsweeps=0" -d "freq_sel_init=0.05" -d "freq_sel_end=0.95" -d "s_beneficial=0.1" -d "ind_sample_size=25" -d "out_sample_size=1" -d "file_output1='./02_slim_SFS_PS.txt'" ./slim_template.slim&
srun --ntasks 1 --exclusive --mem-per-cpu=1GB slim -t -m -d "Ne=1000" -d "L=500000" -d "Neb=1000" -d "mut_rate=1e-6" -d "rec_rate=1e-4" -d "ngenes=100" -d "rate_ben=0.20" -d "rate_del=8" -d "s_backg_ben=0.002" -d "s_backg_del=-0.005" -d "nsweeps=0" -d "freq_sel_init=0.05" -d "freq_sel_end=0.95" -d "s_beneficial=0.1" -d "ind_sample_size=25" -d "out_sample_size=1" -d "file_output1='./03_slim_SFS_BGS_PS.txt'" ./slim_template.slim&
wait
