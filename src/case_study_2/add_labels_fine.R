# set beginning of each dive
knownStates <- rep(NA,nrow(Data_fine))
start_inds <- which(Data_fine$divenum != lag(Data_fine$divenum))
knownStates[start_inds] <- 1

# get average time of each dive
avg_times <- Data_fine %>% 
  group_by(divenum) %>%
  summarise(avg_time = mean(stime))

# get eating in the second half
Data_fine <- left_join(x = Data_fine, y = avg_times)
Data_fine$eating_second_half <- Data_fine$eating * (Data_fine$stime > Data_fine$avg_time)

# foraging is either foraging or eating in the second half of the dive
Data_fine$foraging_signs <- (Data_fine$foraging %in% c(1,2)) | (Data_fine$eating_second_half %in% c(1,2))

dive_sums <- Data_fine %>% 
              group_by(divenum) %>%
              summarise(num_forg = sum(foraging_signs),
                        camon_perc = mean(camon))

pos_dives <- dive_sums$divenum[dive_sums$num_forg > 0]
neg_dives <- dive_sums$divenum[dive_sums$num_forg == 0 & dive_sums$camon_perc == 1.0]

# indiced of dive ends
pos_end_inds <- which(Data_fine$divenum != lead(Data_fine$divenum) & Data_fine$divenum %in% pos_dives)
neg_end_inds <- which(Data_fine$divenum != lead(Data_fine$divenum) & Data_fine$divenum %in% neg_dives)

knownStates[pos_end_inds] <- 6
knownStates[neg_end_inds] <- 5