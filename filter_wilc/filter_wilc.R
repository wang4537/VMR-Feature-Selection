t_test_fs <- function(test_tf)
{
  #uid to genotype name
  test_tf$uid <- unlist(lapply(strsplit(test_tf$uid, "_"), function(x) x[3]))
  test_list <- lapply(as.list(test_tf[,-1]), function(x) data.frame("genotype" = test_tf$uid, "feature" = as.numeric(as.character(x))))
  #number of tests
  #n_test <- ncol(test_tf)-1

  ## caculate p-value and fold change
  #p-value
  test_list_p <- lapply(test_list, function(x) wilcox_test(x, feature~genotype))
  test_list_p.df <- do.call('rbind', test_list_p) %>% adjust_pvalue(method = 'fdr')
  test_list_p <- -log10(test_list_p.df$p.adj)
  #fc
  test_list_FC <- unlist(lapply(test_list, function(x)
  {
    y <- aggregate(feature~genotype, x, mean)
    fc <- y[2,2]/y[1,2]
    return(log2(fc))
  }))
  #combine p-value and fc
  test_result <- data.frame("-log10P" = test_list_p,
                            "log2FC" = test_list_FC) #bind the result
  test_result$feature <- names(test_list) #consistent feature name as in RF

  return(test_result) #return the file for all
}

volcano <- function(test_result, file_name, fc_cutoff)
{
  ## prepare result for plotting ##
  #add no, up and down
  test_result$SigDiff <- "NO"
  test_result$SigDiff[test_result$X.log10P > -log10(0.01) & test_result$log2FC > log2(fc_cutoff)] <-"UP"
  test_result$SigDiff[test_result$X.log10P > -log10(0.01) & test_result$log2FC < -log2(fc_cutoff)] <-"DOWN"
  #calculate the distance to 0
  test_result$norm <- sqrt(test_result$X.log10P^2+ test_result$log2FC^2)
  top_10 <- test_result[test_result$SigDiff != "NO",] %>% arrange(desc(norm)) %>% top_n(30) 
  
  #save the list to files
  write.csv(test_result[test_result$SigDiff != "NO",],
            gsub(".pdf", ".csv", file_name))

  #label if sig
  test_result$label <- NA
  test_result$label[test_result$feature %in% top_10$feature] <- test_result$feature[test_result$feature %in% top_10$feature]
  #volcano plot
  sig_color <- c("blue", "red", "black")
  names(sig_color) <- c("DOWN", "UP", "NO")
  p <- ggplot(test_result, aes_string(x = "log2FC", y = "X.log10P", col = "SigDiff", label="label")) +
    geom_point() +
    scale_color_manual(values = sig_color) +
    geom_vline(xintercept = c(-log2(fc_cutoff),log2(fc_cutoff)), col = "black", linetype = "longdash") +
    geom_hline(yintercept = -log10(0.01), col = "black") +
    geom_text_repel(max.overlaps = 50) +
    ylab("-log10(p)") +
    xlab("log2(FC)") +
    theme(axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))

  pdf(file_name)
  plot(p)
  dev.off()
}

#main
library(ggplot2)
library(ggrepel)
library(rstatix)
library(dplyr)
library(parallel) # multi-core processing

setwd("Y:/data/Feature_selection_framework/filter_wilc/")
#light on
on_files <- list.files("../transformed_data", pattern = "*on.csv", full.names = TRUE, recursive = TRUE)
mat_names <- gsub("../transformed_data/|/.*$", "", on_files) #keep only the mid dir name
on_list <- lapply(on_files, read.csv) #read the files
names(on_list) <- mat_names #rename the list
#create the output directories
sapply(mat_names, function(x) dir.create(paste("./",x, "/", sep = ""))) #create folds for outputs
#read the yellow page
on_yellow_page <- read.csv("../int_output/Opt_on_kFolds_yellow_page.csv", header = TRUE)
#index for training data
on_train.list <- unique(unlist(lapply(on_yellow_page[,3:ncol(on_yellow_page)], function(x) which(x != "Test"))))
#subset the dfs in the list
cl <- makeCluster(4)#doParallel, setup parallel processing
registerDoParallel(cl)
on_list.fs <- foreach(i = 1:length(on_list), .packages = c("dplyr", "rstatix")) %dopar% {t_test_fs(on_list[[i]][on_train.list,])}
names(on_list.fs) <- names(on_list)

#on_list.fs <- lapply(on_list, function(x)
#{
#  x_train <- x[on_train.list,] #subset the training data
#  x_fs <- t_test_fs(x_train) #p-value and fold change

#  return(x_fs)
#})
#plot and save volcano
for(i in 1:length(on_list.fs))
{
  #write the filter values to a file
  full_out <- paste("./", names(on_list.fs)[i], "/", names(on_list.fs)[i], "_Opt_on_Full.csv", sep = "")
  write.csv(on_list.fs[[i]], full_out, quote = FALSE, row.names = FALSE)
  #volcano plot and filtered results
  out_file <- paste("./", names(on_list.fs)[i], "/", names(on_list.fs)[i], "_Opt_on_vol.pdf", sep = "")
  volcano(on_list.fs[[i]], out_file, 1.3)
  rm(out_file)
}

#light off
off_files <- list.files("../transformed_data", pattern = "*off.csv", full.names = TRUE, recursive = TRUE)
mat_names <- gsub("../transformed_data/|/.*$", "", off_files) #keep only the mid dir name
off_list <- lapply(off_files, read.csv) #read the files
names(off_list) <- mat_names #rename the list
#read the yellow page
off_yellow_page <- read.csv("../int_output/Opt_off_kFolds_yellow_page.csv", header = TRUE)
#index for training data
off_train.list <- unlist(lapply(off_yellow_page[,3:ncol(off_yellow_page)], function(x) which(x != "Test")))

cl <- makeCluster(4)#doParallel, setup parallel processing
registerDoParallel(cl)
off_list.fs <- foreach(i = 1:length(off_list), .packages = c("dplyr", "rstatix")) %dopar% {t_test_fs(off_list[[i]][off_train.list,])}
names(off_list.fs) <- names(off_list)
#subset the dfs in the list
#off_list.fs <- lapply(off_list, function(x)
#{
#  x_train <- x[off_train.list,] #subset the training data
#  x_fs <- t_test_fs(x_train) #p-value and fold change
#  return(x_fs)
#})
#plot and save volcano
for(i in 1:length(off_list.fs))
{
  #write the filter values to a file
  full_out <- paste("./", names(on_list.fs)[i], "/", names(on_list.fs)[i], "_Opt_off_Full.csv", sep = "")
  write.csv(on_list.fs[[i]], full_out, quote = FALSE, row.names = FALSE)

  #volcano plot and filtered results
  out_file <- paste("./", names(off_list.fs)[i], "/", names(off_list.fs)[i], "_Opt_off_vol.pdf", sep = "")
  volcano(off_list.fs[[i]], out_file, 1.5)
  rm(out_file)
}
