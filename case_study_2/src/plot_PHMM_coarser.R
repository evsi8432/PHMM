### plot results ###

labs <- c(Elevation = "Depth (m)",
          maxDepth = "Maximum Depth (m)",
          diveDuration = "Dive Duration (s)",
          w_low = "Wiggliness (Low Frequency)",
          w_high = "Wiggliness (High Frequency)",
          logWLow = "Wiggliness (Low Frequency) (log10)",
          logWHigh = "Wiggliness (High Frequency) (log10)",
          logWTotal = "Wiggliness (Total) (log10)",
          postDiveInt = "Post Dive Interval (s)",
          htv = "Heading Total Variation (rad / s)",
          logHtv = "Heading Total Variation (rad) (log10)",
          logJpNorm = "Jerk Peak (normalized, log10)",
          jp_normed = "Jerk Peak (normalized)",
          delt_d = "Change in Depth (m)",
          roll = "Roll (rad)",
          rajp = "Roll @ Peak Jerk (rad)",
          p_catch = "Probability of Catch",
          VeDBA = "VeDBA")

plot_dives <- function(dives){
  
  for(dive in dives){
    
    if(dive %in% pos_dives){
      title0 <- "pos"
    } else if (dive %in% neg_dives){
      title0 <- "neg"
    } else {
      title0 <- "NA"
    }
    
    title1 <- tail(Data_less_fine_unlabelled$p_fish[Data_less_fine_unlabelled$ID == dive],1)
    title1 <- title1 + tail(Data_less_fine_unlabelled$p_catch[Data_less_fine_unlabelled$ID == dive],1)
    title1 <- round(title1, 3)
    
    df <- Data_less_fine_unlabelled[Data_less_fine_unlabelled$ID == dive,]
    df$Elevation <- -df$ad
    stime <- min(df$stime,na.rm=T)
    df$stime <- (df$stime - min(df$stime,na.rm=T))/60
    cols_to_plot <- c("Elevation","p_catch",setdiff(names(dist),c("knownState")))
    df_long <- df %>%
      pivot_longer(cols = cols_to_plot,
                   names_to = "feature")
    
    colors <- hcl(h = seq(15, 375, length = 6 + 1),
                  l = 65, c = 100)[1:6]
    colors <- c("1" = colors[1],
                "2" = colors[2],
                "3" = colors[3],
                "4" = colors[4],
                "5" = colors[5],
                "6" = colors[6])
    
    plot0 <- ggplot(df_long,aes(x=stime, y=value)) +
      geom_line() +
      geom_point(aes(color=factor(viterbi,levels=1:6),
                     group=divenum)) +
      geom_hline(yintercept = 0) +
      geom_vline(aes(xintercept = stime),
                 data = df %>% dplyr::filter(true_label == 4)) +
      labs(title = TeX(paste("$\\alpha =", lambda, "$, ",
                             "$P(X_{s,T_s} \\in \\{4,6\\} \\ | \\ Y_s) =", title1, "$")),
           color="", y="",
           x="Time (min)") +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
      scale_color_manual(labels = c("1" = "descent",
                                    "2" = "bottom",
                                    "3" = "chase",
                                    "4" = "capture",
                                    "5" = "ascent w/o fish",
                                    "6" = "ascent w/ fish"),
                         values = colors) +
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            text = element_text(size=16),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5)) +
      facet_wrap(~feature, ncol = 1, labeller = as_labeller(labs), scales = "free_y")
    
    ggsave(paste0(directory,"/plt/profile_",
                  paste0(setdiff(names(dist),"knownState"),collapse = "_"),"_",
                  round(log10(lambda),3),"_",
                  K,"_",
                  title0,"_",
                  dive,".png"), 
           plot0, 
           width = 6, height = 8)
    print(plot0)
    print(dive)
  }
}

#plot_dives(setdiff(dives,c(pos_dives,neg_dives)),"NA")
#plot_dives(neg_dives,"neg")
#plot_dives(pos_dives,"pos")
plot_dives(test_dives)
