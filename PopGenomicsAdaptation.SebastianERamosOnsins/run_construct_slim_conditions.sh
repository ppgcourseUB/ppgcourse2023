#### RUN 15 DIFFERENT CONDITIONS. OBTAIN SFS IN SYN AND NONSYN. CALCULATE MKTasymptotic. ####
#### GENERAK CONDITIONS:
#### Two populations: ONE target plus ONE outgroup.
#### 1.The initial population run for 5Ne generations to achieve some equilibrium mutation-selection-drift
#### 2.Split target and outgroup for 10Ne generations
#### 3.Possible change in number of individuals in target population after 5Ne+10Ne
#### 4.Possible selective sweeps (n simulataneous sweeps) 0.5Ne before end simulation

#fixed paraneters
Ne=1000; L=50000; ngenes=50;
mut_rate=1e-6;
ind_sample_size=25; out_sample_size=1;

# CONDITION 0:
#Neutral. No change Ne. No sweep
FILEOUT="'./slim_SFS_SNM_00.txt'"
Neb=1000; nsweeps=0;
rec_rate=1e-4;
rate_ben=0; s_backg_ben=0;
rate_del=0; s_backg_del=0;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim > ./run_slim_conditions.sh

# CONDITION 1:
#Strong BACKGROUND SELECTION. No beneficial selection. No change Ne. No sweep
FILEOUT="'./slim_SFS_BGS_01.txt'"
Neb=1000; nsweeps=0;
rec_rate=1e-4
rate_ben=0; s_backg_ben=0;
rate_del=8; s_backg_del=-1.0;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim >> ./run_slim_conditions.sh

# CONDITION 2:
#Weaker BACKGROUND SELECTION. No beneficial selection. No change Ne. No sweep
FILEOUT="'./slim_SFS_WBGS_02.txt'"
Neb=1000; nsweeps=0;
rec_rate=1e-4
rate_ben=0; s_backg_ben=0;
rate_del=8; s_backg_del=-0.25;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim >> ./run_slim_conditions.sh

# CONDITION 3:
#Strong background selection. No beneficial selection. No change Ne. No sweep. NO RECOMBINATION WITHIN
FILEOUT="'./slim_SFS_BGS_R0_03.txt'"
Neb=1000; nsweeps=0;
rec_rate=1e-4
rate_ben=0; s_backg_ben=0;
rate_del=8; s_backg_del=-1.0;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim >> ./run_slim_conditions.sh

# CONDITION 4:
#Strong background selection. No beneficial selection. POPULATION REDUCTION. No sweep.
FILEOUT="'./slim_SFS_BGS_RED_04.txt'"
Neb=200; nsweeps=0;
rec_rate=1e-4
rate_ben=0; s_backg_ben=0;
rate_del=8; s_backg_del=-1.0;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim >> ./run_slim_conditions.sh

# CONDITION 5:
#Strong background selection. No beneficial selection. POPULATION EXPANSION. No sweep.
FILEOUT="'./slim_SFS_BGS_EXP_05.txt'"
Neb=5000; nsweeps=0;
rec_rate=1e-4
rate_ben=0; s_backg_ben=0;
rate_del=8; s_backg_del=-1.0;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim >> ./run_slim_conditions.sh

# CONDITION 6:
#Strong background selection. SMALL PROPORTION BENEFICIAL SELECTION. No change Ne. No sweep.
FILEOUT="'./slim_SFS_BGS_PSL_06.txt'"
Neb=1000; nsweeps=0;
rec_rate=1e-4
rate_ben=0.1; s_backg_ben=0.005;
rate_del=8; s_backg_del=-1.0;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim >> ./run_slim_conditions.sh

# CONDITION 7:
#Strong background selection. MIDDLE PROPORTION BENEFICIAL SELECTION. No change Ne. No sweep.
FILEOUT="'./slim_SFS_BGS_PSM_07.txt'"
Neb=1000; nsweeps=0;
rec_rate=1e-4
rate_ben=0.5; s_backg_ben=0.005;
rate_del=8; s_backg_del=-1.0;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim >> ./run_slim_conditions.sh

# CONDITION 8:
#Strong background selection. HIGH PROPORTION BENEFICIAL SELECTION. No change Ne. No sweep.
FILEOUT="'./slim_SFS_BGS_PSH_08.txt'"
Neb=1000; nsweeps=0;
rec_rate=1e-4
rate_ben=2; s_backg_ben=0.005;
rate_del=8; s_backg_del=-1.0;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim >> ./run_slim_conditions.sh

# CONDITION 9:
#Strong background selection. Beneficial selection. No change Ne. SWEEPS.
FILEOUT="'./slim_SFS_BGS_SW_09.txt'"
Neb=1000; nsweeps=3;
rec_rate=1e-4
rate_ben=0.5; s_backg_ben=0;
rate_del=8; s_backg_del=-1.0;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim >> ./run_slim_conditions.sh

# CONDITION 10:
#Weaker background selection. Beneficial selection. Pop Reduction. No Sweeps
FILEOUT="'./slim_SFS_WBGS_PSH_RED_SW_10.txt'"
Neb=100; nsweeps=3;
rec_rate=1e-4
rate_ben=0.5; s_backg_ben=0.005;
rate_del=8; s_backg_del=-0.25;

echo slim -t -m -d \"Ne=$Ne\" -d \"L=$L\" -d \"Neb=$Neb\" -d \"mut_rate=$mut_rate\" -d \"rec_rate=$rec_rate\" -d \"ngenes=$ngenes\" -d \"rate_ben=$rate_ben\" -d \"rate_del=$rate_del\" -d \"s_backg_ben=$s_backg_ben\" -d \"s_backg_del=$s_backg_del\" -d \"nsweeps=$nsweeps\" -d \"freq_sel_init=0.05\" -d \"freq_sel_end=0.95\" -d \"s_beneficial=0.1\" -d \"ind_sample_size=$ind_sample_size\" -d \"out_sample_size=$out_sample_size\" -d \"file_output1=$FILEOUT\" ./slim_template.slim >> ./run_slim_conditions.sh

