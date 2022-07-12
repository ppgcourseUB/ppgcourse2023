#!/bin/bash
#
#SBATCH --job-name=9slims
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --ntasks=9
#SBATCH --mem=12GB
#SBATCH --partition=normal
#
module load SLiM

srun --ntasks 1 --exclusive --mem-per-cpu=1GB slim -t -m -d "Ne=500" -d "L=500000" -d "Neb=500" -d "mut_rate=1e-6" -d "rec_rate=1e-4" -d "ngenes=100" -d "rate_ben=0" -d "rate_del=0" -d "s_backg_ben=0" -d "s_backg_del=0" -d "nsweeps=0" -d "freq_sel_init=0.05" -d "freq_sel_end=0.95" -d "s_beneficial=0.1" -d "ind_sample_size=25" -d "out_sample_size=1" -d "file_output1='./00_slim_SFS_SNM.txt'" ./slim_template.slim
srun --ntasks 1 --exclusive --mem-per-cpu=1GB slim -t -m -d "Ne=500" -d "L=500000" -d "Neb=500" -d "mut_rate=1e-6" -d "rec_rate=1e-4" -d "ngenes=100" -d "rate_ben=0" -d "rate_del=8" -d "s_backg_ben=0" -d "s_backg_del=-0.1" -d "nsweeps=0" -d "freq_sel_init=0.05" -d "freq_sel_end=0.95" -d "s_beneficial=0.1" -d "ind_sample_size=25" -d "out_sample_size=1" -d "file_output1='./01_slim_SFS_BGS.txt'" ./slim_template.slim
srun --ntasks 1 --exclusive --mem-per-cpu=1GB slim -t -m -d "Ne=500" -d "L=500000" -d "Neb=500" -d "mut_rate=1e-6" -d "rec_rate=1e-4" -d "ngenes=100" -d "rate_ben=0.005" -d "rate_del=0" -d "s_backg_ben=0.005" -d "s_backg_del=0" -d "nsweeps=0" -d "freq_sel_init=0.05" -d "freq_sel_end=0.95" -d "s_beneficial=0.1" -d "ind_sample_size=25" -d "out_sample_size=1" -d "file_output1='./02_slim_SFS_WBGS.txt'" ./slim_template.slim
srun --ntasks 1 --exclusive --mem-per-cpu=1GB slim -t -m -d "Ne=500" -d "L=500000" -d "Neb=200" -d "mut_rate=1e-6" -d "rec_rate=1e-4" -d "ngenes=100" -d "rate_ben=0.05" -d "rate_del=8" -d "s_backg_ben=0" -d "s_backg_del=-1.0" -d "nsweeps=0" -d "freq_sel_init=0.05" -d "freq_sel_end=0.95" -d "s_beneficial=0.1" -d "ind_sample_size=25" -d "out_sample_size=1" -d "file_output1='./03_slim_SFS_BGS_RED.txt'" ./slim_template.slim
srun --ntasks 1 --exclusive --mem-per-cpu=1GB slim -t -m -d "Ne=500" -d "L=500000" -d "Neb=2500" -d "mut_rate=1e-6" -d "rec_rate=1e-4" -d "ngenes=100" -d "rate_ben=0.05" -d "rate_del=8" -d "s_backg_ben=0" -d "s_backg_del=-1.0" -d "nsweeps=0" -d "freq_sel_init=0.05" -d "freq_sel_end=0.95" -d "s_beneficial=0.1" -d "ind_sample_size=25" -d "out_sample_size=1" -d "file_output1='./04_slim_SFS_BGS_EXP.txt'" ./slim_template.slim
srun --ntasks 1 --exclusive --mem-per-cpu=1GB slim -t -m -d "Ne=500" -d "L=500000" -d "Neb=500" -d "mut_rate=1e-6" -d "rec_rate=1e-4" -d "ngenes=100" -d "rate_ben=0.05" -d "rate_del=8" -d "s_backg_ben=0.005" -d "s_backg_del=-1.0" -d "nsweeps=0" -d "freq_sel_init=0.05" -d "freq_sel_end=0.95" -d "s_beneficial=0.1" -d "ind_sample_size=25" -d "out_sample_size=1" -d "file_output1='./05_slim_SFS_BGS_PSL.txt'" ./slim_template.slim
srun --ntasks 1 --exclusive --mem-per-cpu=1GB slim -t -m -d "Ne=500" -d "L=500000" -d "Neb=500" -d "mut_rate=1e-6" -d "rec_rate=1e-4" -d "ngenes=100" -d "rate_ben=0.5" -d "rate_del=8" -d "s_backg_ben=0.005" -d "s_backg_del=-1.0" -d "nsweeps=0" -d "freq_sel_init=0.05" -d "freq_sel_end=0.95" -d "s_beneficial=0.1" -d "ind_sample_size=25" -d "out_sample_size=1" -d "file_output1='./06_slim_SFS_BGS_PSM.txt'" ./slim_template.slim
srun --ntasks 1 --exclusive --mem-per-cpu=1GB slim -t -m -d "Ne=500" -d "L=500000" -d "Neb=500" -d "mut_rate=1e-6" -d "rec_rate=1e-4" -d "ngenes=100" -d "rate_ben=2" -d "rate_del=8" -d "s_backg_ben=0.005" -d "s_backg_del=-1.0" -d "nsweeps=0" -d "freq_sel_init=0.05" -d "freq_sel_end=0.95" -d "s_beneficial=0.1" -d "ind_sample_size=25" -d "out_sample_size=1" -d "file_output1='./07_slim_SFS_BGS_PSH.txt'" ./slim_template.slim
srun --ntasks 1 --exclusive --mem-per-cpu=1GB slim -t -m -d "Ne=500" -d "L=500000" -d "Neb=500" -d "mut_rate=1e-6" -d "rec_rate=1e-4" -d "ngenes=100" -d "rate_ben=0.5" -d "rate_del=8" -d "s_backg_ben=0.005" -d "s_backg_del=-1.0" -d "nsweeps=10" -d "freq_sel_init=0.05" -d "freq_sel_end=0.95" -d "s_beneficial=0.1" -d "ind_sample_size=25" -d "out_sample_size=1" -d "file_output1='./08_slim_SFS_BGS_PSM_SW.txt'" ./slim_template.slim
wait
