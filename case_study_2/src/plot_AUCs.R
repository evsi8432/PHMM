# make plots of AUCs
means <- c()
sds <- c()
AUCs <- c()

param_directory <- "params"#"params_unshared_label_1"
plt_directory <- "plt"#"plt_unshared_label_1"
K <- "4"

# load in base AUCS
filename <- paste0(param_directory,"/AUC_base_",K,".csv")
AUCs <- c(AUCs,data.frame(fread(filename))[,2])
means <- c(means,mean(data.frame(fread(filename))[,2]))
sds <- c(sds,sd(data.frame(fread(filename))[,2]))

# load in other AUCs
for(log_alpha in c("-4","-3","-2","-1","0")){
  filename <- paste0(param_directory,"/AUC_delt_d_htv_jp_normed_",log_alpha,"_",K,".csv")
  AUCs <- c(AUCs,data.frame(fread(filename))[,2])
  means <- c(means,mean(data.frame(fread(filename))[,2]))
  sds <- c(sds,sd(data.frame(fread(filename))[,2]))
}

# make dataframe
df <- data.frame(alpha = rep(c("Baseline","0.0001","0.001","0.01","0.1","1"),each = K),
                 baseline = rep(c(T,F,F,F,F,F),each = K),
                 fold = rep(1:K,6),
                 AUC = AUCs)

df_means <- data.frame(alpha = c("Baseline","0.0001","0.001","0.01","0.1","1"),
                       baseline = c(T,F,F,F,F,F),
                       mean = means,
                       sd = sds)
# ggplot
plot0 <- ggplot(df, aes(x = alpha, y = AUC)) + 
  #geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 4, alpha = 0.5,
             width = 0.05, height = 0.00) +
  #geom_line(aes(group = fold),alpha = 0.5) +
  geom_point(data = df_means, 
             aes(x = alpha, y = mean),
             shape = "*", size = 20) + 
  geom_line(data = df_means, 
             aes(x = alpha, y = mean, group = baseline)) + 
  #geom_errorbar(data = df_means,
  #              aes(x = Alpha, y = mean, ymin=mean-sd, ymax=mean+sd), width=.2,
  #              position=position_dodge(0.05)) +
  geom_vline(aes(xintercept = 5.5)) +
  labs(x=TeX("\\alpha"),y="AUC") +
  theme_classic() +
  theme(text = element_text(size=16))

#plot0

#plot0 <- ggplot(df_means, aes(x=Alpha, y=mean, group = 1)) + 
#  geom_line() +
#  geom_point() +
#  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
#                position=position_dodge(0.05))

#plot0

ggsave(paste0(directory,"/plt/AUCs_",
              paste0(setdiff(names(dist),"knownState"),collapse = "_"),"_",
              round(log10(lambda),3),"_",
              K,".png"),
       plot0,
       width = 6, height = 4)
