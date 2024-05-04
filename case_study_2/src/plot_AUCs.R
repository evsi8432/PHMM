# make plots of AUCs
means <- c()
sds <- c()
AUCs <- c()

# load in base AUCS
filename <- "params/AUC_base_5.csv"
AUCs <- c(AUCs,data.frame(fread(filename))[,2])
means <- c(means,mean(data.frame(fread(filename))[,2]))
sds <- c(sds,sd(data.frame(fread(filename))[,2]))

# load in other AUCs
for(log_alpha in c("-Inf","-3","-2","-1","-0.301","0")){
  filename <- paste0("params/AUC_delt_d_rajp_htv_jp_normed_",log_alpha,"_5.csv")
  AUCs <- c(AUCs,data.frame(fread(filename))[,2])
  means <- c(means,mean(data.frame(fread(filename))[,2]))
  sds <- c(sds,sd(data.frame(fread(filename))[,2]))
}

# make dataframe
df <- data.frame(Alpha = rep(c("Baseline","0","0.001","0.01","0.1","0.5","1"),each = 5),
                 AUC = AUCs)

df_means <- data.frame(Alpha = c("Baseline","0","0.001","0.01","0.1","0.5","1"),
                       mean = means,
                       sd = sds)
# ggplot
plot0 <- ggplot(df, aes(x = Alpha, y = AUC)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1, alpha = 0.5,
              width=0.1, height=0) +
  geom_point(data = df_means, 
             aes(x = Alpha, y = mean),
             shape = "*", size = 5) + 
  geom_vline(aes(xintercept = 6.5)) +
  labs(x=TeX("\\alpha"),y="AUC") +
  theme_classic() +
  theme(text = element_text(size=16))
  #facet_wrap(~feature, ncol = 1, labeller = as_labeller(labs), scales = "free_y")

ggsave("plt/AUCs_delt_d_rajp_htv_jp_normed_4.png",
       plot0,
       width = 6, height = 4)
