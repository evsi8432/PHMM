Data_for_test <- Data_fine[Data_fine$divenum %in% test_dives,]

# get jerk peak
dive_maxs <- Data_for_test %>% 
  group_by(divenum) %>%
  summarize(jp_normed = max(jp_normed),
            max_depth = max(ad))

# get roll at jerk peak
dive_maxs <- left_join(dive_maxs,
                       Data_for_test[c("divenum","jp_normed","rajp")])

names(dive_maxs) <- c("divenum","max_jp","max_depth","global_rajp")

# get average heading total variation at bottom
Data_for_test <- left_join(Data_for_test,
                            dive_maxs[c("divenum","max_depth")])

Data_for_test$bottom <- Data_for_test$ad > 0.7*Data_for_test$max_depth
Data_for_test_bot <- Data_for_test[Data_for_test$bottom,]

dive_htvs <- Data_for_test_bot %>% 
  group_by(divenum) %>%
  summarize(avg_htv_bot = mean(htv))


probs <- c()
labs <- c()

conf_matrices_base[[k]] = matrix(0, nrow = 2, ncol = 2)
rownames(conf_matrices_base[[k]]) <- c("True Positive", "True Negative")
colnames(conf_matrices_base[[k]]) <- c("Pred Positive", "Pred Negative")

for(divenum in test_dives){
  
  if(divenum %in% pos_dives){
    rownum <- 1
    labs <- c(labs,TRUE)
  } else if (divenum %in% neg_dives){
    rownum <- 2
    labs <- c(labs,FALSE)
  }
  
  htv_check <- dive_htvs$avg_htv_bot[dive_htvs$divenum == divenum] >= base_model$htv
  jp_check <- dive_maxs$max_jp[dive_maxs$divenum == divenum] >= base_model$jp_normed
  rajp_check <- abs(dive_maxs$global_rajp[dive_maxs$divenum == divenum]) >= base_model$rajp
  
  print(divenum)
  print(divenum %in% pos_dives)
  print(paste("htv:",htv_check))
  print(paste("jp:",jp_check))
  print(paste("rajp:",rajp_check))
  
if(htv_check & jp_check & rajp_check){
    colnum <- 1
    probs <- c(probs,1)
  } else {
    colnum <- 2
    probs <- c(probs,0)
  }
  
  conf_matrices_base[[k]][rownum,colnum] = conf_matrices_base[[k]][rownum,colnum] + 1
}

probs_base[[k]] <- probs
AUCs_base[k] <- roc(response = labs, predictor=probs, direction = "<")$auc
plot(roc(response = labs, predictor=probs))
