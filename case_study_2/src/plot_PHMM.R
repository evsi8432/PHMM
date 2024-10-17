### plot results ###
options(scipen=999)

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

plot_dives <- function(dives,df){
  
  for(dive in dives){
    
    if(dive %in% pos_dives){
      title0 <- "pos"
    } else if (dive %in% neg_dives){
      title0 <- "neg"
    } else {
      title0 <- "NA"
    }
    
    title1 <- tail(df$p_fish[df$ID == dive],1)
    title1 <- title1 + tail(df$p_catch[df$ID == dive],1)
    title1 <- round(title1, 3)
    
    dive_df <- df[df$ID == dive,]
    dive_df$Elevation <- -dive_df$ad
    stime <- min(dive_df$stime,na.rm=T)
    dive_df$stime <- (dive_df$stime - min(dive_df$stime,na.rm=T))/60
    cols_to_plot <- c("Elevation","p_catch",setdiff(names(dist),c("knownState")))
    dive_df_long <- dive_df %>%
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
    
    if(lambda == 1.0){
      legend.pos <- "right"
      wid <- 8
    } else {
      legend.pos <- "none"
      wid <- 6
    }
    
    plot0 <- ggplot(dive_df_long,aes(x=stime, y=value)) +
      geom_line() +
      geom_point(aes(color=factor(viterbi,levels=1:6),
                     group=divenum,
                     shape=viterbi==4),
                 size = 2) +
      geom_hline(yintercept = 0) +
      geom_vline(aes(xintercept = stime),
                 data = dive_df %>% dplyr::filter(true_label == 4)) +
      labs(title = TeX(paste0("$\\alpha = ", lambda, "$")),
           color="", y="",
           x="Time (min)") +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
      scale_color_manual(limits = c("1","2","3","4","5","6"),
                         labels = c("1" = "descent",
                                    "2" = "bottom",
                                    "3" = "chase",
                                    "4" = "capture",
                                    "5" = "ascent w/o fish",
                                    "6" = "ascent w/ fish"),
                         values = colors) +
      scale_shape_manual(values = c(19,8),
                         guide = "none") +
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            text = element_text(size=16),
            legend.position = legend.pos,
            plot.title = element_text(hjust = 0.5)) +
      facet_wrap(~feature, ncol = 1, labeller = as_labeller(labs), scales = "free_y")
    
    ggsave(paste0(directory,"/plt/profile_",
                  paste0(setdiff(names(dist),"knownState"),collapse = "_"),"_",
                  round(log10(lambda),3),"_",
                  K,"_",
                  title0,"_",
                  dive,".png"), 
           plot0, 
           width = wid, height = 8)
    print(plot0)
    print(dive)
  }
}

if(K == 1){
  plot_dives(unique(df_unlabelled$ID),df_unlabelled)
} else {
  plot_dives(test_dives,df_unlabelled)
}

labs <- c(htv = "Heading Total Variation (rad/s, log)",
          jp_normed = "Jerk Peak (normalized, log)",
          delt_d = "Change in Depth (m)",
          rajp = "Roll @ Peak Jerk (rad)",
          p_catch = "Probability of Catch")

xs <- c()
states <- c()
features <- c()
values <- c()
npoints <- 1000

xlims <- list("delt_d" = c(-10,10),
              "htv"    = c(-5,2),
              "jp_normed" = c(-3,6))

colors <- hcl(h = seq(15, 375, length = 6 + 1),
              l = 65, c = 100)[1:6]

colors <- c("1" = colors[1],
            "2" = colors[2],
            "3" = colors[3],
            "4" = colors[4],
            "5" = colors[5])

if(K == 1){
  
  # make dataframe of densities
  for(feature in c("delt_d","htv","jp_normed")){
    features <- c(features,rep(feature,5*npoints))
    for (state in 1:5){
      states <- c(states,rep(state,npoints))
      
      xs0 <- seq(xlims[[feature]][1],
                 xlims[[feature]][2],
                 length.out=npoints)
      
      mu <- hmm$mle[[feature]][1,state]
      sd <- hmm$mle[[feature]][2,state]
      
      if(feature %in% c("htv","jp_normed")){
        shape <- (mu/sd)^2
        rate <- mu/(sd^2)
        values0 <- dgamma(exp(xs0),shape,rate) * exp(xs0)
      } else {
        values0 <- dnorm(xs0,mean = mu, sd = sd)
      }
      xs <- c(xs,xs0)
      values <- c(values,values0)
    }
  }

  if(lambda == 1.0){
    legend.pos <- "right"
    wid <- 8
  } else {
    legend.pos <- "none"
    wid <- 6
  }
  
  df_to_plot <- data.frame(x = xs,
                           feature = features,
                           state = factor(states),
                           value = values)
  
  lim_df <- data.frame(x = xs,
                       feature = features,
                       state = factor(states),
                       value = rep(1.0,length(xs)))
  
  plot0 <- ggplot(df_to_plot, aes(x = x, y = value, 
                         group = state,
                         color = state)) +
    geom_line() + 
    geom_blank(data=lim_df) +
    labs(title = TeX(paste("$\\alpha =", lambda, "$")),
         color="", y="Density",
         x="") +
    scale_color_manual(limits = c("1","2","3","4","5"),
                       labels = c("1" = "descent",
                                  "2" = "bottom",
                                  "3" = "chase",
                                  "4" = "capture",
                                  "5" = "ascent"),
                       values = colors) +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          text = element_text(size=16),
          legend.position = legend.pos,
          plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~feature,ncol = 1,labeller = as_labeller(labs), scales = "free")
  
  ggsave(paste0(directory,"/plt/emissions_",
                paste0(setdiff(names(dist),"knownState"),collapse = "_"),"_",
                round(log10(lambda),3),"_",
                K,".png"), 
         plot0, 
         width = wid, height = 6)
}
