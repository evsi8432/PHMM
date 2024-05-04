library(momentuHMM)
library(dplyr)
library(mclust)
library(data.table)
library(ggplot2)
library(latex2exp)

#setwd("/Users/evsi8432/Documents/Research/PHMM/src/bash")

# get command-line arguments
args <- commandArgs(trailingOnly=TRUE)
args <- c("logMDDD_1-1-1_dd-30_2023-10-23.R",NA)

opt_file <- args[1]
args <- as.integer(args[2])

# get options
source(paste0('../opt/',opt_file))

# Set Seed
set.seed(1)

# define models
source("../preprocessing/load_data.R")
ratio <- sum(!(Data$knownState %in% 4)) / nrow(Data)
ratio <- round(ratio,3)
models <- list()
#models[[1]] <- c("no",     1.0)
#models[[2]] <- c("random", 1.0)
models[[1]] <- c("fixed",  0.0)       # no weight
models[[2]] <- c("fixed",  0.5*ratio) 
models[[3]] <- c("fixed",  ratio)     # equal weight
models[[4]] <- c("fixed",  0.5 + 0.5*ratio)
models[[5]] <- c("fixed",  1.0)       # natural weight
n_models <- length(models) 

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

whales <- c("none","A100a","A100b","A113a","A113b",
            "D21a","D21b","D26a","D26b",
            "I107a","I107b","I129","I145a","I145b",
            "L87","L88","R48a","R48b","R58a","R58b")
whales <- c()

if("Male" %in% sex){
  whales <- c(whales,"D21a","D21b","I107a","I107b","L87","L88")
}
if("Female" %in% sex) {
  whales <- c(whales,"A100a","A100b","A113a","A113b","D26a","D26b",
              "I129","I145a","I145b","R48a","R48b","R58a","R58b")
}

# get metrics from all models
df <- data.frame(model = c(),
                 behaviour = c(),
                 metric = c(),
                 value = c())

for(model in models){
  
  # extract model and lambda
  lamb <- model[2]
  model <- model[1]
  
  # get conf matrix for each whale
  conf_matrix <- matrix(rep(0,3*3),nrow=3,ncol=3)
    
  for(holdout_whale in whales){
  
    conf_matrix_whale <- read.csv(make_title(paste0(directory,"/params/"),
                                             paste0(model,"-",
                                                    lamb,"-",
                                                    holdout_whale,"-",
                                                    "confusion_matrix.csv")))
    
    
    print(conf_matrix_whale)
    conf_matrix <- conf_matrix + conf_matrix_whale[1:3,2:4]
  }
  
  rownames(conf_matrix) <- c("True Resting", "True Travelling", "True Foraging")
  colnames(conf_matrix) <- c("Predicted Resting", "Predicted Travelling", "Predicted Foraging")
  
  write.csv(conf_matrix,make_title(paste0(directory,"/params/"),
                                   paste0(model,"-",
                                          lamb,"-",
                                          "confusion_matrix_all.csv")))
  
  se_r <- sum(conf_matrix[ 1, 1]) / sum(conf_matrix[ 1,])
  sp_r <- sum(conf_matrix[-1,-1]) / sum(conf_matrix[-1,])
  se_t <- sum(conf_matrix[ 2, 2]) / sum(conf_matrix[ 2,])
  sp_t <- sum(conf_matrix[-2,-2]) / sum(conf_matrix[-2,])
  se_f <- sum(conf_matrix[ 3, 3]) / sum(conf_matrix[ 3,])
  sp_f <- sum(conf_matrix[-3,-3]) / sum(conf_matrix[-3,])
  
  df_model <- data.frame(model = rep(lamb,6),
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
  labs(x = "Behaviour",
       y = "",
       fill = TeX("$\\alpha$")) +
  facet_wrap(~metric,ncol=1)

ggsave(make_title(paste0(directory,"/plt/"),
                  "model_comparison.png"),
       plot = plot0,
       width = 4,
       height = 4,
       device='png', 
       dpi=500)
