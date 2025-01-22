l1norm_preprocess <- function(x)
{
  #calculate the distance
  x$L1Norm <- apply(x[,6:14], 1, function(y) as.numeric(norm(as.matrix(y),type = '1'))) #l1-norm of a vector 
  #calculte the mean
  x.ave <- aggregate(L1Norm ~ genotype+start, x, mean)
  #calcualte the standard deviation
  x.STD <- aggregate(L1Norm ~ genotype+start, x, sd)
  #add the standard error to the mean df
  SEM <- (x.STD$L1Norm)/(sqrt(nrow(x)/length(unique(x$start))))
  x.ave$SEM <- as.numeric(SEM)
  #add another column for timepoints
  x.ave$timepoint <- x.ave$start -  min(x.ave$start)
  return(x.ave)
}

l2norm_preprocess <- function(x)
{
  #calculate the distance
  x$L2Norm <- apply(x[,6:14], 1, function(y) norm(y,type = '2')) #l2-norm of a vector 
  #calculte the mean
  x.ave <- aggregate(L2Norm ~ genotype+start, x, mean)
  #calcualte the standard deviation
  x.STD <- aggregate(L2Norm ~ genotype+start, x, sd)
  #add the standard error to the mean df
  SEM <- (x.STD$L2Norm)/(sqrt(nrow(x)/length(unique(x$start))))
  x.ave$SEM <- as.numeric(SEM)
  #add another column for timepoints
  x.ave$timepoint <- x.ave$start -  min(x.ave$start)
  return(x.ave)
}

duration_preprocess <- function(x)
{
  #calculate the duration
  x$duration <- apply(x[,c('inadur', 'smldur', 'lardur')], 1, function(y) as.numeric(norm(as.matrix(y),type = '1'))) #total duration 
  #calculte the mean
  x.ave <- aggregate(duration ~ genotype+start, x, mean)
  #calcualte the standard deviation
  x.STD <- aggregate(duration ~ genotype+start, x, sd)
  #add the standard error to the mean df
  SEM <- (x.STD$duration)/(sqrt(nrow(x)/length(unique(x$start))))
  x.ave$SEM <- as.numeric(SEM)
  #add another column for timepoints
  x.ave$timepoint <- x.ave$start -  min(x.ave$start)
  return(x.ave)
}

distance_preprocess <- function(x)
{
  #calculate the distance
  x$distance <- apply(x[,c('inadist', 'smldist', 'lardist')], 1, function(y) as.numeric(norm(as.matrix(y),type = '1')))
  #calculte the mean
  x.ave <- aggregate(distance ~ genotype+start, x, mean)
  #calcualte the standard deviation
  x.STD <- aggregate(distance ~ genotype+start, x, sd)
  #add the standard error to the mean df
  SEM <- (x.STD$distance)/(sqrt(nrow(x)/length(unique(x$start))))
  x.ave$SEM <- as.numeric(SEM)
  #add another column for timepoints
  x.ave$timepoint <- x.ave$start -  min(x.ave$start)
  return(x.ave)
}

L1Norm_plot <- function(dat.ave, light)
{
  #specify several plot parameters
  break_points <- c(0,(1:10)*30)
  label_points <- c(0,(1:10)*30)-30
  dat.type <- dat.ave$genotype
  plot_color <- c("Red", "Black")
  #ymax change based on ligth conditions
  ymax <- if(light=="on") c(180, 225) else c(0, 4)
  #ggplot object
  response <- ggplot(data=dat.ave, aes_string(x="timepoint", y="L1Norm")) +
    labs(y = "Average L1 Norm (A.U.)", x = "time (s)") +
    geom_line(aes(color = dat.type), size = 0.6, alpha = 1) +
    geom_ribbon(aes(ymin = L1Norm - SEM, ymax = L1Norm + SEM, fill = dat.type), alpha = 0.2) +
    scale_color_manual(limits = unique(dat.type), values = plot_color,
                       labels = unique(dat.type)) +
    scale_fill_manual(limits = unique(dat.type), values = plot_color,
                      labels = unique(dat.type))+
    scale_x_continuous(breaks=break_points, labels=label_points, expand = c(0,0)) +
    coord_cartesian(ylim = ymax) +
    theme_bw(base_size = 32) +
    theme(panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.major = element_line(colour = "#bababa",size = 0.1),
          axis.text.x = element_text(color="#000000"), axis.text.y = element_text(color="#000000"),
          legend.position = "top",legend.title=element_blank(), legend.text = element_text(size = 15), legend.spacing = unit(1,"cm"),
          axis.title.x = element_text(size =20), axis.title.y = element_text(size =20, vjust = 0.3),
          plot.margin = unit(c(0, 1, 0.5, 0.5), "cm"))
}

L2Norm_plot <- function(dat.ave, light)
{
  #specify several plot parameters
  break_points <- c(0,(1:10)*30)
  label_points <- c(0,(1:10)*30)-30
  dat.type <- dat.ave$genotype
  plot_color <- c("Red", "Black")
  #ymax change based on ligth conditions
  ymax <- if(light=="on") c(180, 225) else c(0, 4)
  #ggplot object
  response <- ggplot(data=dat.ave, aes_string(x="timepoint", y="L2Norm")) +
    labs(y = "Average L2 Norm (A.U.)", x = "time (s)") +
    geom_line(aes(color = dat.type), size = 0.6, alpha = 1) +
    geom_ribbon(aes(ymin = L2Norm - SEM, ymax = L2Norm + SEM, fill = dat.type), alpha = 0.2) +
    scale_color_manual(limits = unique(dat.type), values = plot_color,
                       labels = unique(dat.type)) +
    scale_fill_manual(limits = unique(dat.type), values = plot_color,
                      labels = unique(dat.type))+
    scale_x_continuous(breaks=break_points, labels=label_points, expand = c(0,0)) +
    coord_cartesian(ylim = ymax) +
    theme_bw(base_size = 32) +
    theme(panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.major = element_line(colour = "#bababa",size = 0.1),
          axis.text.x = element_text(color="#000000"), axis.text.y = element_text(color="#000000"),
          legend.position = "top",legend.title=element_blank(), legend.text = element_text(size = 15), legend.spacing = unit(1,"cm"),
          axis.title.x = element_text(size =20), axis.title.y = element_text(size =20, vjust = 0.3),
          plot.margin = unit(c(0, 1, 0.5, 0.5), "cm"))
}

duration_plot <- function(dat.ave)
{
  #specify several plot parameters
  break_points <- c(0,(1:10)*30)
  label_points <- c(0,(1:10)*30)-30
  dat.type <- dat.ave$genotype
  plot_color <- c("Red", "Black")
  
  #ggplot object
  response <- ggplot(data=dat.ave, aes_string(x="timepoint", y="duration")) +
    labs(y = "Average duration (s)", x = "time (s)") +
    geom_line(aes(color = dat.type), size = 0.6, alpha = 1) +
    geom_ribbon(aes(ymin = duration - SEM, ymax = duration + SEM, fill = dat.type), alpha = 0.2) +
    scale_color_manual(limits = unique(dat.type), values = plot_color,
                       labels = unique(dat.type)) +
    scale_fill_manual(limits = unique(dat.type), values = plot_color,
                      labels = unique(dat.type))+
    scale_x_continuous(breaks=break_points, labels=label_points, expand = c(0,0)) +
    coord_cartesian(ylim = c(0,1)) +
    theme_bw(base_size = 32) +
    theme(panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.major = element_line(colour = "#bababa",size = 0.1),
          axis.text.x = element_text(color="#000000"), axis.text.y = element_text(color="#000000"),
          legend.position = "top",legend.title=element_blank(), legend.text = element_text(size = 15), legend.spacing = unit(1,"cm"),
          axis.title.x = element_text(size =20), axis.title.y = element_text(size =20, vjust = 0.3),
          plot.margin = unit(c(0, 1, 0.5, 0.5), "cm"))
}

distance_plot <- function(dat.ave)
{
  #specify several plot parameters
  break_points <- c(0,(1:10)*30)
  label_points <- c(0,(1:10)*30)-30
  dat.type <- dat.ave$genotype
  plot_color <- c("Red", "Black")
  
  #ggplot object
  response <- ggplot(data=dat.ave, aes_string(x="timepoint", y="distance")) +
    labs(y = "Average Distance (cm)", x = "time (s)") +
    geom_line(aes(color = dat.type), size = 0.6, alpha = 1) +
    geom_ribbon(aes(ymin = distance - SEM, ymax = distance + SEM, fill = dat.type), alpha = 0.2) +
    scale_color_manual(limits = unique(dat.type), values = plot_color,
                       labels = unique(dat.type)) +
    scale_fill_manual(limits = unique(dat.type), values = plot_color,
                      labels = unique(dat.type))+
    scale_x_continuous(breaks=break_points, labels=label_points, expand = c(0,0)) +
    coord_cartesian(ylim = c(0,0.4)) +
    theme_bw(base_size = 32) +
    theme(panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.major = element_line(colour = "#bababa",size = 0.1),
          axis.text.x = element_text(color="#000000"), axis.text.y = element_text(color="#000000"),
          legend.position = "top",legend.title=element_blank(), legend.text = element_text(size = 15), legend.spacing = unit(1,"cm"),
          axis.title.x = element_text(size =20), axis.title.y = element_text(size =20, vjust = 0.3),
          plot.margin = unit(c(0, 1, 0.5, 0.5), "cm"))
}

plot_wrapper <- function(dat.list, light, x)
{
  
  if(x == 'duration')
  {
    x_ave <- do.call("rbind", lapply(dat.list, duration_preprocess))
    x.plot <- duration_plot(x_ave)
    #save the plot
    out_file <- paste("output/", light, "_ave_", x,".pdf", sep = "")
  }
  if(x == 'distance')
  {
    x_ave <- do.call("rbind", lapply(dat.list, distance_preprocess))
    x.plot <- distance_plot(x_ave)
    #save the plot
    out_file <- paste("output/", light, "_ave_", x,".pdf", sep = "") 
  }
  pdf(out_file)
  print(x.plot)
  dev.off()
  
  return(x.plot)
}
#### Main ####
library(ggplot2)
setwd("/depot/yleung/data/Feature_selection_framework/ave_plot/")
#read all normalized data into a list
on <- read.csv("../input_data/Opt_on.csv")
off <- read.csv("../input_data/Opt_off.csv")
index_set <- c('duration', 'distance')
#on 
on <- on[on$start %in% (-30):(299), ]
on_list <- split(on, f = as.factor(on$genotype))
on_plot.list <- lapply(index_set, function(x) plot_wrapper(on_list,"on", x))

#read the data
off <- off[off$start %in% (-30):(299), ]
off_list <- split(off, f = as.factor(off$genotype))
off_plot.list <- lapply(index_set, function(x) plot_wrapper(off_list,"off", x))

#take L1 Norm seperately
on_ave <- do.call("rbind", lapply(on_list, l1norm_preprocess))
on.plot <- L1Norm_plot(on_ave, "on")
#save the plot
out_file <- paste("output/on_ave_L1Norm.pdf", sep = "")
pdf(out_file)
plot(on.plot)
dev.off()

#take L1 Norm seperately
off_ave <- do.call("rbind", lapply(off_list, l1norm_preprocess))
off.plot <- L1Norm_plot(off_ave, "off")
#save the plot
out_file <- paste("output/off_ave_L1Norm.pdf", sep = "")
pdf(out_file)
plot(off.plot)
dev.off()

#take L2 Norm seperately
on_ave <- do.call("rbind", lapply(on_list, l2norm_preprocess))
on.plot <- L2Norm_plot(on_ave, "on")
#save the plot
out_file <- paste("output/on_ave_L2Norm.pdf", sep = "")
pdf(out_file)
plot(on.plot)
dev.off()

#take L2 Norm seperately
off_ave <- do.call("rbind", lapply(off_list, l2norm_preprocess))
off.plot <- L2Norm_plot(off_ave, "off")
#save the plot
out_file <- paste("output/off_ave_L2Norm.pdf", sep = "")
pdf(out_file)
plot(off.plot)
dev.off()
