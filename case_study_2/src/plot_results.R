

df <- data.frame(model = rep(c("Tennessen","PHMM"),each=2),
                 metric = rep(c("Precision","Recall"),2),
                 value = c(prc_base,rcl_base,prc_PHMM,rcl_PHMM))

plot0 <- ggplot(df,
                aes(x=metric,
                    y=value,
                    fill=model)) +
  geom_col(position="dodge") #+
#facet_wrap(~metric,ncol=1)

model_name <- paste0("model_comparison_",
                     setdiff(names(dist),"knownState"),
                     ".rds")

ggsave(make_title(paste0(directory,"/plt/"),
                  "model_comparison.png"),
       plot = plot0,
       width = 8,
       height = 8,
       device='png', 
       dpi=500)
