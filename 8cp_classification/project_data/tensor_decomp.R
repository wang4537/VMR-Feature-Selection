#FUNCTION: kneedle algorithm to detect elbow point
kneedle <- function(errors, range) #x = error, range = rank list
{
  x_norm <- (range - min(range)) / (max(range) - min(range)) #normalized x
  y_norm <- (errors - min(errors)) / (max(errors) - min(errors)) #normalized y
  distances <- x_norm - y_norm  #negative: y decreases larger x
  # Find knees (peaks in distances)
  optimal_ranks <- range[which.max(distances)] #max value
  return(optimal_ranks)
}
#FUNCTION: calulate cp decomposition and prepare training and testing data
cp_train_test <- function(f_name)
{
  #### prepare file format #####
  #import data
  train_file <- paste('../../2preprocess_sep/', f_name, "_train_norm.csv", sep = "")
  train_data <- read.csv(train_file)
  #sort by genotype and time
  train_data <- train_data[order(train_data$genotype, train_data$start, train_data$location),]
  train_data <- train_data[train_data$start >= 0, ] #data after light change
  #import data
  test_file <- paste('../../2preprocess_sep/', f_name, "_test_norm.csv", sep = "")
  test_data <- read.csv(test_file)
  #sort by genotype and time
  test_data <- test_data[order(test_data$genotype, test_data$start, test_data$location),]
  test_data <- test_data[test_data$start >= 0, ] #data after light change
  
  #### Training set: CP decomposition ####
  #spit by time points
  train_list <- split(train_data, f = as.factor(train_data$start))
  #keep only start and functional variables
  train_list_sub <- lapply(train_list, function(x) as.matrix(x[,6:14]))
  train_array <- simplify2array(train_list_sub)

  #run CP through a range of ranks
  rank_range <- 1:9
  cl <- makeCluster(3)#doParallel, setup parallel processing
  registerDoParallel(cl)
  rank_list_cp <- foreach(i = rank_range, .packages = c("rTensor")) %dopar% {cp(as.tensor(train_array), num_components = i, max_iter = 100)} #CP at each rank number
  rank_error <- unlist(lapply(rank_list_cp, function(x) x$fnorm_resid)) #get the error
  opt_rank <- kneedle(rank_error, rank_range) #rund the kneedle algorithm
  cat("Max rank number:", opt_rank, "\n") #print the rank
  cat("Optimal rank:", opt_rank, "\n") #print the rank

  #run the CP with optimal rank
  cp_decomp_opt <- cp(as.tensor(train_array), num_components = opt_rank, max_iter = 500)
  #extract the projected data from the training set, bind with the original uid and genotype
  train_genotype <- unlist(lapply(strsplit(unique(train_data$uid), split = "_"), function(x) x[3]))
  train_cp <- cbind(unique(train_data$uid), train_genotype, cp_decomp_opt$U[[1]])
  proj_col_names <- c("uid", "genotype", paste("Rank", 1:ncol(cp_decomp_opt$U[[1]]),sep = ""))#rename the columns
  colnames(train_cp) <- proj_col_names
  #save the training data to a csv file
  csv_file <- paste(f_name, "cp_train_data.csv", sep = "_")
  write.csv(train_cp, csv_file, quote = FALSE, row.names = FALSE)

  #### Project Testing Data ####
  #convert the test data to 3-D tensor
  #spit by time points
  test_list <- split(test_data, f = as.factor(test_data$start))
  #keep only start and functional variables
  test_list_sub <- lapply(test_list, function(x) as.matrix(x[,6:14]))
  test_array <- simplify2array(test_list_sub)

  #projection using Moore–Penrose pseudoinverse
  A_train = cp_decomp_opt$U[[1]]
  B = cp_decomp_opt$U[[2]]
  C = cp_decomp_opt$U[[3]]
  lambda_matrix_inv <- diag(1/cp_decomp_opt$lambdas) #inverse of diagonal matrix
  mode1_unfolded <- k_unfold(as.tensor(test_array), m = 1)@data #mode-1 unfolding of test tensor
  test_cp <- mode1_unfolded %*% khatri_rao(C, B) %*% ginv( (t(C)%*%C) * (t(B)%*%B ) ) %*% lambda_matrix_inv #calculate projection

  test_genotype <- unlist(lapply(strsplit(unique(test_data$uid), split = "_"), function(x) x[3])) #test genotype
  test_cp <- cbind(unique(test_data$uid), test_genotype, test_cp) #add uid and genotype columns
  colnames(test_cp) <- proj_col_names #rename columns
  #save the data to a file
  csv_file <- paste(f_name, "cp_test_data.csv", sep = "_")
  write.csv(test_cp, csv_file, quote = FALSE, row.names = FALSE)

  return(list("CP" = cp_decomp_opt, "proj_train" = train_cp, "proj_test" = test_cp)) #retrun a list of CP object, projected train and test set
}

library(rTensor)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)# multi-core processing
library(caret) # cross validation and data partition
library(MASS) # calculate ginv, Moore–Penrose pseudoinverse
setwd('/depot/yleung/data/Feature_selection_framework/8cp_classification/project_data/')
#output_folder <- "tensor_decomposition/output_on/"
#dir.create(output_folder)#create the output folder if not existing
on_cp_list <- cp_train_test("Opt_on")
saveRDS(on_cp_list, "on_cp_obj_data.rds")
off_cp_list <- cp_train_test("Opt_off")
saveRDS(off_cp_list, "off_cp_obj_data.rds")
