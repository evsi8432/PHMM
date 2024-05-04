conf_matrix_PHMM = matrix(0, nrow = 2, ncol = 2)
rownames(conf_matrix_PHMM) <- c("True Positive", "True Negative")
colnames(conf_matrix_PHMM) <- c("Pred Positive", "Pred Negative")

conf_matrix_base = matrix(0, nrow = 2, ncol = 2)
rownames(conf_matrix_base) <- c("True Positive", "True Negative")
colnames(conf_matrix_base) <- c("Pred Positive", "Pred Negative")

for(k in 1:K){
  conf_matrix_base <- conf_matrix_base + conf_matrices_base[[k]]
  conf_matrix_PHMM <- conf_matrix_PHMM + conf_matrices_PHMM[[k]]
}

prc_base <- conf_matrix_base["True Positive","Pred Positive"] / sum(conf_matrix_base[,"Pred Positive"])
rcl_base <- conf_matrix_base["True Positive","Pred Positive"] / sum(conf_matrix_base["True Positive",])

prc_PHMM <- conf_matrix_PHMM["True Positive","Pred Positive"] / sum(conf_matrix_PHMM[,"Pred Positive"])
rcl_PHMM <- conf_matrix_PHMM["True Positive","Pred Positive"] / sum(conf_matrix_PHMM["True Positive",])

write.csv(conf_matrix_base,paste0(directory,"/params/conf_mat_base.csv"))
write.csv(conf_matrix_PHMM,paste0(directory,"/params/conf_mat_",
                                  paste0(setdiff(names(dist),"knownState"),collapse = "_"),"_",
                                  round(log10(lambda),3),"_",
                                  K,".csv"))

write.csv(AUCs_base,paste0(directory,"/params/AUC_base_",K,".csv"))
write.csv(AUCs_PHMM,paste0(directory,"/params/AUC_",
                                      paste0(setdiff(names(dist),"knownState"),collapse = "_"),"_",
                                      round(log10(lambda),3),"_",
                                      K,".csv"))

# make plots of AUCs


