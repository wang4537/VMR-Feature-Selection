library(rTensor)
library(ggplot2)
library(gridExtra)
#library(parallel)
library(doParallel)
library(foreach)# multi-core processing

setwd('Y:/data/Feature_selection_framework/')
output_folder <- "tensor_decomposition/output_off/"
dir.create(output_folder)#create the output folder if not existing
#import data
off_data <- read.csv('input_data/Opt_off.csv')
#sort by genotype and time
off_data <- off_data[order(off_data$genotype, off_data$start, off_data$location),]
off_data <- off_data[off_data$start >= 0, ]
#transform long to 3d array
off_list <- split(off_data, f = as.factor(off_data$start))
#keep only start and functional variables
off_list_sub <- lapply(off_list, function(x) as.matrix(x[,6:14]))
off_array <- simplify2array(off_list_sub)

#select the best rank
cl <- makeCluster(6)#doParallel, setup parallel processing
registerDoParallel(cl)
rank_list_cp <- foreach(i = 1:50, .packages = c("rTensor")) %dopar% {cp(as.tensor(off_array), num_components = i, max_iter = 100)}
rank_error <- unlist(lapply(rank_list_cp, function(x) x$fnorm_resid))
pdf(paste(output_folder, "error_plot_off.pdf"))
plot(rank_error, type = "b")
dev.off()

#use rank 1:3
cp_decomp_opt <- cp(as.tensor(off_array), num_components = 3, max_iter = 500)
rank_num <- 1:3
#plot by samples
sample_plot <- function(cp_obj, n_rank)
{
  plot_dat <- data.frame("score" = cp_obj$U[[1]][,n_rank])
  plot_dat <- cbind(plot_dat, "sample" = 1:nrow(plot_dat))
  plot_dat$genotype <- sort(rep(c("Q344X", "WT"), nrow(plot_dat)/2))
  dat.type <- plot_dat$genotype
  plot_color <- c("Red", "Black")

  p <- ggplot(data = plot_dat, aes_string(x = "sample", y = "score", color = "genotype")) +
    geom_point() +
    scale_color_manual(limits = unique(dat.type), values = plot_color,
                       labels = unique(dat.type)) +
    geom_vline(xintercept = nrow(plot_dat)/2 + 0.5, color = "blue", size = 0.2) +
    labs(title = paste("Rank = ", n_rank, sep = "")) +
    theme(plot.title = element_text(size=20),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))

  return(p)
}
sample_plot_list <- lapply(rank_num, function(x) sample_plot(cp_decomp_opt, x))
pdf(paste(output_folder, "sample_tensor.pdf"), onefile = TRUE)
do.call('grid.arrange', sample_plot_list)
dev.off()
#plot by behaviors
behavior_plot <- function(cp_obj, n_rank)
{
  plot_dat <- data.frame("score" = cp_obj$U[[2]][,n_rank])
  plot_dat <- cbind(plot_dat, "behavior" = 1:nrow(plot_dat))
  plot_dat$behavior <- colnames(off_data)[6:14]
  p <- ggplot(data = plot_dat, aes_string(x = "behavior", y = "score")) +
    geom_bar(stat = "identity", width = 0.5) +
    labs(title = paste("Rank = ", n_rank, sep = "")) +
    theme(plot.title = element_text(size=20),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))
  return(p)
}
behavior_plot_list <- lapply(rank_num, function(x) behavior_plot(cp_decomp_opt, x))
pdf(paste(output_folder, "behavior_tensor.pdf"), onefile = TRUE)
do.call('grid.arrange', behavior_plot_list)
dev.off()

#plot by samples
time_plot <- function(cp_obj, n_rank)
{
  plot_dat <- data.frame("score" = cp_obj$U[[3]][,n_rank])
  plot_dat <- cbind(plot_dat, "time" = 1:nrow(plot_dat))
  p <- ggplot(data = plot_dat, aes_string(x = "time", y = "score")) +
    geom_line() +
    labs(title = paste("Rank = ", n_rank, sep = "")) +
    theme(plot.title = element_text(size=20),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))
  return(p)
}
time_plot_list <- lapply(rank_num, function(x) time_plot(cp_decomp_opt, x))
pdf(paste(output_folder, "time_tensor.pdf"), onefile = TRUE)
do.call('grid.arrange', time_plot_list)
dev.off()
