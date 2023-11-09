# load in packages
module load StdEnv/2020
module load gcc/9.3.0
module load r/4.3.1
module load proj/9.0.1
module load gdal/3.5.1
module load udunits/2.2.28

# fit models
out1=$(sbatch fit_model_PHMM.sh $1 --output=/dev/null)
sar1=($out1)
jid1=(${sar1[3]})

# fit models that failed to run before
out2=$(sbatch --dependency=afterany:$jid1 fit_model_PHMM_cleanup.sh $1 --output=/dev/null)
sar2=($out2)
jid2=(${sar2[3]})

# then plot the results
sbatch --dependency=afterany:$jid2 plot_PHMM.sh $1 --output=/dev/null

# make confusion matrices
out3=$(sbatch --dependency=afterany:$jid2 conf_mat_PHMM.sh $1 --output=/dev/null)
sar3=($out3)
jid3=(${sar3[3]})

# Summarize the resulting models
sbatch --dependency=afterany:$jid3 summarize_model_PHMM.sh $1 --output=/dev/null

# change permissions
lfs find ~/projects/*/ -group evsi8432
chown -h -R evsi8432:def-nheckman -- ~/projects/def-nheckman/evsi8432/
lfs find ~/projects/def-nheckman/evsi8432 -type d -print0 | xargs -0 chmod g+s
