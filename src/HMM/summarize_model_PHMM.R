library(momentuHMM)
library(dplyr)
library(mclust)
library(data.table)
library(ggplot2)

#setwd("/Users/evsi8432/Documents/Research/PHMM/src/bash")

# get command-line arguments
args <- commandArgs(trailingOnly=TRUE)

opt_file <- args[1]
args <- as.integer(args[2])

# get options
source(paste0('../opt/',opt_file))

# Set Seed
set.seed(1)

# Select Model
models <- c("no","fixed_1","fixed_2","half_random","random")

dir.create(directory, showWarnings = FALSE)
dir.create(paste0(directory,"/params"), showWarnings = FALSE)
dir.create(paste0(directory,"/plt"), showWarnings = FALSE)

make_title <- function(start,end){
  title <- paste0(start,statesPerBehaviour[1])
  for(nstates in statesPerBehaviour[2:3]){
    title <- paste0(title,"-",nstates)
  }
  for(feature in names(dist)){
    title <- paste0(title,"-",feature)
  }
  if(length(sex) > 1){
    title <- paste0(title,"_all")
  } else {
    title <- paste0(title,"_",sex)
  }
  title <- paste0(title,"_",end)
  return(title)
}

whales <- c()
if("Male" %in% sex){
  whales <- c(whales,"D21","I107","L87","L88")
}
if("Female" %in% sex) {
  whales <- c(whales,"A100","A113","D26","I129","I145","R48","R58")
}

# get metrics from all models
df <- data.frame(model = c(),
                 behaviour = c(),
                 metric = c(),
                 value = c())


for(model in models){
  
  # get conf matrix for each whale
  conf_matrix <- matrix(rep(0,3*3),nrow=3,ncol=3)
    
  for(holdout_whale in whales){
  
    conf_matrix_whale <- read.csv(make_title(paste0(directory,"/params/"),
                                             paste0(model,"-",
                                                    holdout_whale,"-",
                                                    "confusion_matrix.csv")))
    
    
    print(conf_matrix_whale)
    conf_matrix <- conf_matrix + conf_matrix_whale[1:3,2:4]
  }
  
  rownames(conf_matrix) <- c("True Resting", "True Travelling", "True Foraging")
  colnames(conf_matrix) <- c("Predicted Resting", "Predicted Travelling", "Predicted Foraging")
  
  write.csv(conf_matrix,make_title(paste0(directory,"/params/"),
                                   paste0(model,"-",
                                          "confusion_matrix_all.csv")))
  
  se_r <- sum(conf_matrix[ 1, 1]) / sum(conf_matrix[ 1,])
  sp_r <- sum(conf_matrix[-1,-1]) / sum(conf_matrix[-1,])
  se_t <- sum(conf_matrix[ 2, 2]) / sum(conf_matrix[ 2,])
  sp_t <- sum(conf_matrix[-2,-2]) / sum(conf_matrix[-2,])
  se_f <- sum(conf_matrix[ 3, 3]) / sum(conf_matrix[ 3,])
  sp_f <- sum(conf_matrix[-3,-3]) / sum(conf_matrix[-3,])
  
  df_model <- data.frame(model = rep(model,6),
                         behaviour = rep(c("Resting","Travelling","Foraging"),each=2),
                         metric = rep(c("Sensitivity","Specificity"),3),
                         value = c(se_r,sp_r,se_t,sp_t,se_f,sp_f))
  
  df <- rbind(df,df_model)

}

plot0 <- ggplot(df,
                aes(x=behaviour,
                    y=value,
                    fill=model)) +
  geom_col(position="dodge") +
  facet_wrap(~metric,ncol=1)

ggsave(make_title(paste0(directory,"/plt/"),
                  "model_comparison.png"),
       plot = plot0,
       width = 8,
       height = 8,
       device='png', 
       dpi=500)
