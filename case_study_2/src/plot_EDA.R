plot <- F
plot_raw <- F

labs <- c(Elevation = "Depth (meters)",
          maxDepth = "Maximum Depth (m)",
          diveDuration = "Dive Duration (s)",
          w_low = "Wiggliness (Low Frequency)",
          w_high = "Wiggliness (High Frequency)",
          logWLow = "Wiggliness (Low Frequency) (log10)",
          logWHigh = "Wiggliness (High Frequency) (log10)",
          logWTotal = "Wiggliness (Total) (log10)",
          postDiveInt = "Post Dive Interval (s)",
          htv = "Heading Total Variation (m)",
          delt_d = "Change in Depth (m)",
          roll = "Roll (radians)",
          VeDBA = "VeDBA")

if(plot){
  
  # plot scatterplots
  plot1 <- ggplot(Data_fine,
                  aes(x=w_low,
                      y=w_high,
                      color=as.factor(foraging_signs))) +
    geom_point() +
    stat_density_2d(color="white") +
    scale_x_log10() + scale_y_log10() +
    labs(color="Labelled Behaviour") +
    facet_wrap(~foraging_signs , ncol = 3)
  
  print(plot1)
  
  plot1 <- ggplot(Data_fine,
                  aes(x=htv,
                      y=delt_d,
                      color=as.factor(foraging_signs))) +
    geom_point() +
    stat_density_2d(color="white") +
    scale_x_log10() +
    labs(color="Labelled Behaviour") +
    facet_wrap(~foraging_signs , ncol = 3)
  
  print(plot1)
  
  plot1 <- ggplot(Data_fine,
                  aes(x=abs_roll,
                      y=jp,
                      color=as.factor(foraging_signs))) +
    geom_point() +
    stat_density_2d(color="white") +
    scale_y_log10() +
    labs(color="Labelled Behaviour") + 
    facet_wrap(~foraging_signs , ncol = 3)
  
  print(plot1)
  
  plot1 <- ggplot(Data_fine,
                  aes(x=roll,
                      y=pitch,
                      color=as.factor(foraging_signs))) +
    geom_point() +
    stat_density_2d(color="white") +
    labs(color="Labelled Behaviour") + 
    facet_wrap(~foraging_signs , ncol = 3)
  
  print(plot1)
  
  plot1 <- ggplot(Data_fine,
                  aes(x=rtv,
                      y=ptv,
                      color=as.factor(foraging_signs))) +
    geom_point() +
    stat_density_2d(color="white") +
    scale_x_log10() + scale_y_log10() +
    labs(color="Labelled Behaviour") +
    facet_wrap(~foraging_signs , ncol = 3)
  
  print(plot1)
  
  for(dive in dives){
    df <- Data_fine[Data_fine$divenum == dive,]
    
    df$foraging_signs[is.na(df$foraging_signs)] <- "none"
    df$Elevation <- -df$ad
    stime <- min(df$stime,na.rm=T)
    df$stime <- (df$stime - min(df$stime,na.rm=T))/60
    
    cols_to_plot <- c("Elevation","delt_d","htv","VeDBA")
    df_long <- df %>%
      pivot_longer(cols = cols_to_plot,
                   names_to = "feature")
    
    plot0 <- ggplot(df_long,aes(x=stime, y=value)) +
      geom_line(aes(color=as.factor(foraging_signs),
                    group=divenum)) +
      geom_hline(yintercept = 0) +
      labs(title = paste(df$ID[1],
                         "camon:",df$camon[3],
                         "foraging:",sum(df$forage %in% c(1,2)),
                         "eating:",sum(df$eating %in% c(1,2))),
           color="", y="",
           x="Elapsed time (minutes)") +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            text = element_text(size=16)) +
      facet_wrap(~feature, ncol = 1, labeller = as_labeller(labs), scales = "free_y")
    
    print(plot0)
    print(dive)
  }
}

if(plot_raw){
  
  # import raw Data
  rawData <- data.frame(fread('../../dat/Final_rawData_Beth.csv'))
  rawData <- rawData[rawData$divenum %in% dives,]
  rawData$elev <- -rawData$p
  
  # define color scheme
  colors <- hcl(h = seq(15, 375, length = 6),
                l = 65, c = 100)[1:5]
  colors <- c(colors,"#000000")
  
  colors <- c("echo.rapid" = colors[1],
              "chase" = colors[2],
              "crunch" = colors[3],
              "scales" = colors[4],
              "foraging" = colors[5],
              "none" = colors[6])
  
  #cols_to_plot <- c("head","pitch","roll","elev")
  #col_title <- "hpr"
  #cols_to_plot <- c("Aw_1","Aw_2","Aw_3","elev")
  #col_title <- "Aw"
  cols_to_plot <- c("dyn.aw1","dyn.aw2","dyn.aw3","elev")
  col_title <- "dyn.w"
  #cols_to_plot <- c("VeDBA","elev","d_elev")
  #col_title <- "elev"
  
  for(dive in c(pos_dives,neg_dives)){
    
    print(dive)
    
    if(dive %in% pos_dives){
      titl <- "positive"
    } else if (dive %in% neg_dives){
      titl <- "negative"
    } else {
      titl <- "unknown"
    }
    
    df <- rawData[rawData$divenum == dive,]
    
    w <- min(50,floor(nrow(df)/2))
    df$d_elev <- c(NA,diff(stats::filter(df$elev, rep(1 / w, w), sides = 2)))
    
    # add foraging signs to raw Data
    df$foraging_signs <- NA
    Data_fine0 <- Data_fine[Data_fine$divenum == dive & 
                              !is.na(Data_fine$foraging_signs),
                            c("stime","etime","foraging_signs")]
    
    for(seg_ind in 1:nrow(Data_fine0)){
      stime <- Data_fine0$stime[seg_ind]
      etime <- Data_fine0$etime[seg_ind]
      foraging_sign <- Data_fine0$foraging_signs[seg_ind]
      df$foraging_signs[df$Time >= stime & df$Time <= etime] <- foraging_sign
      df$foraging_signs[df$Time >= stime & df$Time <= etime] <- foraging_sign
    }
    
    df <- df %>% pivot_longer(cols = cols_to_plot,
                              names_to = "feature")
    
    stime <- min(df$Time)
    #df$Time <- df$Time - stime
    
    plot1 <- ggplot(df,aes(x=Time, y=value, color=foraging_signs)) +
      geom_path(aes(group = 1)) +
      geom_hline(yintercept = 0) +
      scale_color_manual(values = colors) +
      scale_x_datetime(labels = scales::label_date("%H:%M:%S")) +
      labs(color="", y="",
           x="Local Time",
           title=paste(df$ID[1],dive,stime,titl)) +
      theme_classic() +
      theme(strip.background = element_blank(),
            strip.placement = "outside",
            text = element_text(size=16)) +
      facet_wrap(~feature, ncol = 1, scales = "free_y")
    
    ggsave(paste0("plt/",
                  df$ID[1],"_",
                  dive,"_",
                  col_title,
                  ".png"),
           plot1,
           width = 8, height = 6)
    
    print(plot1)
    
    for(axs in c()){
      
      acc <- rawData[rawData$divenum == dive,axs]
      
      # create spectrogram
      spec = specgram(x = acc,
                      n = 256,
                      Fs = 50
      )
      
      # discard phase information
      P = abs(spec$S)
      
      # normalize
      max_P <- max(P)
      #P = P/max_P
      
      # convert to dB
      P = 10*log10(P)
      
      # config time axis
      t = spec$t
      
      max_freq = 25
      eind <- ceiling(nrow(P) * max_freq / 25)
      
      # plot spectrogram
      imagep(x = t,
             y = spec$f[1:eind],
             z = t(P[1:eind,]),
             col = oce.colorsViridis,
             ylab = 'Frequency [Hz]',
             xlab = 'Time [s]',
             zlab = paste(axs,dive,"max intensity",max_P),
             zlim = c(-12,12),
             drawPalette = T,
             decimate = F)
    }
  }
}