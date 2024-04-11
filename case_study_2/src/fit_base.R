Data_for_model <- Data_fine[Data_fine$divenum %in% train_dives &
                            Data_fine$divenum %in% pos_dives,]

# get max depths
dive_maxs <- Data_for_model %>% 
  group_by(divenum) %>%
  summarize(max_depth = max(ad))
Data_for_model <- left_join(Data_for_model,
                            dive_maxs[c("divenum","max_depth")])

# get bottom phase
Data_for_model$bottom <- Data_for_model$ad > 0.7*Data_for_model$max_depth
Data_for_model <- Data_for_model[Data_for_model$bottom,]

# get jerk peak and heading total variation
dive_sums <- Data_for_model %>% 
  group_by(divenum) %>%
  summarize(avg_htv_bot = mean(htv),
            jp_normed = max(jp_normed),
            max_depth = max(ad))

# get roll at jerk peak
dive_sums <- left_join(dive_sums,
                       Data_for_model[c("divenum","jp_normed","rajp")])

# summarize results in list (this defines the model)
base_model <- list("htv" = min(dive_sums$avg_htv_bot),
                   "jp_normed" = min(dive_sums$jp_normed),
                   "rajp" = min(abs(dive_sums$rajp)))
