#cross validation function for NB
NB_k_acc <- function(train, valid, NB_k)
{
  auc_score <- apply(NB_k, 1, function(i)
  {
    X <- as.matrix(train[,-1])
    Y <- train[,1]
    #build the naive bayes classifier
    nb_classifier <- naive_bayes(X,
                                 Y,
                                 usekernel = i[1],
                                 laplace = i[2])
    valid_pred <- predict(nb_classifier, as.matrix(valid[,-1]))
    #calculate the confusion matrix
    cm <- confusionMatrix(valid_pred, as.factor(valid[,1]), positive = "2")
    #get performance metrics
    auroc <- auc(valid[,1], as.numeric(valid_pred)) #area under ROC curve
    return(as.numeric(auroc))
  })
  return(auc_score)
}
#wrapper for CV on mat method and feature sets
MAT_CV <- function(i) # i = matricization method name
{
  input <- input_list[[i]] #GLOBAL VAR: input_list
  feature_set <- feature.list[[i]] #GLOBAL VAR: feature sets
  feature_set$full <- colnames(input[,-1]) #full feature set
  feature_set <- feature_set[lapply(feature_set, length)>0] #remove empty vectors
  #replace the uid by genotype
  input$uid <- as.factor(as.numeric(factor(yellow_page[,2]))) #GLOBAL VAR: yellow_page

  ## cross validation
  ave_cv_feature.list <- lapply(feature_set, function(y) #loop through feature sets
  {
    CV_list <- lapply(fold.list, function(k) #loop through folds
    {
      input <- input[,c("uid",y)]# parse the feature set
      #split the traning and testing set
      train_fold <- input[-c(k,out_test),]
      valid_fold <- input[k,]
      #calculate the CV accuracy
      score_df <- as.data.frame(NB_k_acc(train_fold, valid_fold, NB_parameters))
      return (score_df) #return the mean accuracy of 10 folds in each parameter
    }) # multicore processing

    ave_CV_df <- cbind(NB_parameters, 'AUROC' = as.numeric(rowMeans(as.data.frame(CV_list)))) #add HPs to the df
    ave_CV_df_max <- ave_CV_df[which.max(ave_CV_df[,'AUROC'])[1],] #take the first row with the maximum ave AUROC
    return(ave_CV_df_max)
  })
  ave_cv_feature.df <- do.call("rbind", ave_cv_feature.list) #bind the dfs by rows
  ave_cv_feature.df <- cbind('mat' = rep(i, nrow(ave_cv_feature.df)), #matricization method
                             'feature' = row.names(ave_cv_feature.df), #feature set
                             ave_cv_feature.df) #add a column to indicate the feature set
  return(ave_cv_feature.df)
}

################ MAIN ################
library(naivebayes) # naive bayes
library(caret) # data splitting and grid search
#library(dplyr)
library(pROC) # auroc
library(doParallel)
library(foreach)# multi-core processing
setwd("/depot/yleung/data/Feature_selection_framework/CV_wilc/")
#create the output directories
class_name <- 'NB'
dir.create(class_name)
#list of hyperparameters
NB_parameters <- expand.grid(usekernel = c(TRUE, FALSE), laplace = 0:3)

################ ON ################
#read the data
on_files <- list.files("../transformed_data", pattern = "*on.csv", full.names = TRUE, recursive = TRUE)
input_list <- lapply(on_files, read.csv)
mat_names <- gsub("../transformed_data/|/.*$", "", on_files) #keep only the mid dir name
names(input_list) <- mat_names #rename the list

#### FILTER #####
#read filter lists
on_filter_files <- list.files("../filter_wilc", pattern = "*on_vol.csv", full.names = TRUE, recursive = TRUE)
on_filter.list <- lapply(on_filter_files, function(x) read.csv(x)$feature) #read the files and subset feature names
names(on_filter.list) <- mat_names#rename the list

#### EMBEDDED ####
#read embedded lists
on_embedded_files <- list.files("../embedded_wilc", pattern = "*on_sigVar_pass_list.csv", full.names = TRUE, recursive = TRUE)[-c(1,2)]
on_embedded.list <- lapply(on_embedded_files, function(x) read.csv(x)[,1]) #read the files and subset feature names
names(on_embedded.list) <- mat_names #rename the list

#### FEATURE LIST #####
#create a list all feature sets
feature.list <- list()
for(i in 1:length(on_filter.list))
{
  feature.list[[i]] <- list('full' = '',
                               'filter' = on_filter.list[[i]],
                               'embedded' = on_embedded.list[[i]],
                               'intersection' = intersect(on_filter.list[[i]], on_embedded.list[[i]]),
                               'union' = union(on_filter.list[[i]], on_embedded.list[[i]]))
}
names(feature.list) <- mat_names


#### Set Up Cross Validation and Yellow Page ####
f_name <- "Opt_on"
yellow_page_file_name <- paste("../int_output/", f_name, "_kFolds_yellow_page.csv", sep = "")
yellow_page <- read.csv(yellow_page_file_name, header = TRUE)
fold.list <- lapply(yellow_page[,3:ncol(yellow_page)], function(x) which(x == "Valid"))
out_test <- unlist(lapply(yellow_page[,3:ncol(yellow_page)], function(x) which(x == "Test")))

#one data, loop through feature sets and hyperparameters
cl <- makeCluster(4)#doParallel, setup parallel processing
registerDoParallel(cl)
mat_cv_list <- foreach(i = 1:length(mat_names), .packages = c("naivebayes", "pROC", "caret")) %dopar% {MAT_CV(mat_names[i])}
#mat_cv_list <- lapply(mat_names, MAT_CV)
names(mat_cv_list) <- mat_names
mat_cv_df <- do.call('rbind', mat_cv_list)
mat_cv_df <- cbind('classfier' = rep(class_name, nrow(mat_cv_df)),
                   'light' = rep(f_name, nrow(mat_cv_df)),
                   mat_cv_df)
#save the CV result
out_file <- paste("./", class_name, "/", class_name, "_", f_name, "_CV.csv",sep = '')
write.csv(mat_cv_df, out_file, quote = FALSE, row.names = FALSE)

################ OFF ################
#read the data
off_files <- list.files("../transformed_data", pattern = "*off.csv", full.names = TRUE, recursive = TRUE)
input_list <- lapply(off_files, read.csv)
mat_names <- gsub("../transformed_data/|/.*$", "", off_files) #keep offly the mid dir name
input_list <- lapply(off_files, read.csv) #read the files
names(input_list) <- mat_names #rename the list

#### FILTER #####
#read filter lists
off_filter_files <- list.files("../filter_wilc", pattern = "*off_vol.csv", full.names = TRUE, recursive = TRUE)
off_filter.list <- lapply(off_filter_files, function(x) read.csv(x)$feature) #read the files and subset feature names
names(off_filter.list) <- mat_names#rename the list

#### EMBEDDED ####
#read embedded lists
off_embedded_files <- list.files("../embedded_wilc", pattern = "*off_sigVar_pass_list.csv", full.names = TRUE, recursive = TRUE)[-c(1,2)]
off_embedded.list <- lapply(off_embedded_files, function(x) read.csv(x)[,1]) #read the files and subset feature names
names(off_embedded.list) <- mat_names #rename the list

#### FEATURE LIST #####
#create a list all feature sets
feature.list <- list()
for(i in 1:length(off_filter.list))
{
  feature.list[[i]] <- list('full' = '',
                            'filter' = off_filter.list[[i]],
                            'embedded' = off_embedded.list[[i]],
                            'intersectioff' = intersect(off_filter.list[[i]], off_embedded.list[[i]]),
                            'union' = union(off_filter.list[[i]], off_embedded.list[[i]]))
}
names(feature.list) <- mat_names


#### Set Up Cross Validatioff and Yellow Page ####
f_name <- "Opt_off"
yellow_page_file_name <- paste("../int_output/", f_name, "_kFolds_yellow_page.csv", sep = "")
yellow_page <- read.csv(yellow_page_file_name, header = TRUE)
fold.list <- lapply(yellow_page[,3:ncol(yellow_page)], function(x) which(x == "Valid"))
out_test <- unlist(lapply(yellow_page[,3:ncol(yellow_page)], function(x) which(x == "Test")))

#offe data, loop through feature sets and hyperparameters
cl <- makeCluster(4)#doParallel, setup parallel processing
registerDoParallel(cl)
mat_cv_list <- foreach(i = 1:length(mat_names), .packages = c("naivebayes", "pROC", "caret")) %dopar% {MAT_CV(mat_names[i])}
mat_cv_list <- lapply(mat_names, MAT_CV)
names(mat_cv_list) <- mat_names
#mat_cv_list <- mclapply(mat_names, MAT_CV, mc.cores = 7)
mat_cv_df <- do.call('rbind', mat_cv_list)
mat_cv_df <- cbind('classfier' = rep(class_name, nrow(mat_cv_df)),
                   'light' = rep(f_name, nrow(mat_cv_df)),
                   mat_cv_df)
#save the CV result
out_file <- paste("./", class_name, "/", class_name, "_", f_name, "_CV.csv",sep = '')
write.csv(mat_cv_df, out_file, quote = FALSE, row.names = FALSE)
