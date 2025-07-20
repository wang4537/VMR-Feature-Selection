#cross validation function for 3nn
XGB_acc <- function(train, valid, i)
{
  #construct the input for the xgb model
  dtrain <- list(data = data.matrix(train[,-1]),
                 label = train[,1])
  dvalid <- list(data = data.matrix(valid[,-1]),
                 label = valid[,1])

  #build the model
  classifier <- xgboost(data = dtrain$data,
                        label = dtrain$label,
                        eta = i[1],
                        gamma = i[2],
                        max_depth = i[3],
                        nrounds = 10,
                        objective = "binary:hinge"
  )

  #predict the test labels
  valid_pred <- predict(classifier, dvalid$data)
  #calculate the confusion matrix
  cm <- confusionMatrix(as.factor(valid_pred), as.factor(dvalid$label), positive = "1")
  #get performance metrics
  auroc <- auc(valid[,1], as.numeric(valid_pred)) #area under ROC curve
  #get performance metrics
  acc <- cm$overall[1] #accuracy
  sns <- cm$byClass[1] #sensitivity
  sps <- cm$byClass[2] #specificity
  prc <- cm$byClass[5] #precision
  kappa <- cm$overall[2] #Cohen's kappa
  #output a vector of metrics
  performance <- c(i[1],
                   i[2],
                   i[3],
                   "AUROC" = auroc,
                   acc,
                   sns,
                   sps,
                   prc,
                   kappa)
   #calculate the CV accuracy
   score_df <- data.frame(performance)
   return(score_df)
}

#wrapper for CV on mat method and feature sets
MAT_CV <- function(i) # i = matricization method name
{
  input_train <- input_list[[i]] #GLOBAL VAR: input_list
  input_test <- test_list[[i]]
  feature_set <- feature.list[[i]] #GLOBAL VAR: feature sets
  feature_set$full <- colnames(input_train[,-1]) #full feature set
  feature_set <- feature_set[lapply(feature_set, length)>0] #remove empty vectors
  #replace the uid by genotype
  train_gen <- sapply(strsplit(input_train$uid, "_"), "[", 3)
  input_train$uid <- as.factor(as.numeric(factor(train_gen))) 
  test_gen <- sapply(strsplit(input_test$uid, "_"), "[", 3)
  input_test$uid <- as.factor(as.numeric(factor(test_gen)))
  
  #encode the labels
  levels(input_train$uid) <- c(0,1) #1 is WT, 0 is Q344X
  levels(input_test$uid) <- c(0,1) #1 is WT, 0 is Q344X
  input_train$uid <- as.numeric(as.character(input_train$uid))
  input_test$uid <- as.numeric(as.character(input_test$uid))
  
  #subset the cv output by mat name
  cv_out_sub <- cv_out[cv_out$mat == i,] #GLOBAL VAR: cv_out

  ## cross validation
  per_feature.list <- lapply(names(feature_set), function(y) #loop through feature sets
  {
    train_fold <- input_train[,c("uid",feature_set[[y]])]# parse the feature set
    test_fold <- input_test[,c("uid",feature_set[[y]])]
    #hyperparameters
    XGB_parameters <- cv_out_sub[cv_out_sub$feature == y , 5:7]

    #calculate the CV accuracy
    score_df <- as.data.frame(XGB_acc(train_fold, test_fold, XGB_parameters))
    return (score_df) #return the mean accuracy of 10 folds in each parameter
  }) # multicore processing

  names(per_feature.list) <- names(feature_set)#rename the list
  per_feature.df <- do.call("rbind", per_feature.list) #bind the dfs by rows
  per_feature.df <- cbind('mat' = rep(i, nrow(per_feature.df)), #matricization method
                          'feature' = row.names(per_feature.df), #feature set
                          per_feature.df) #add a column to indicate the feature set
  return(per_feature.df)
}

library(xgboost) # XGB
library(caret) # data splitting and feature selection grid
library(pROC) # auc
library(parallel) # multi-core processing
setwd("/depot/yleung/data/Feature_selection_framework/6HV_preprocess_sep/")
#create the output directories
class_name <- 'XGB'
dir.create(class_name)

################ ON ################
#read the data
on_train_files <- list.files("../3transformed_data_sep", pattern = "*on_train.csv", full.names = TRUE, recursive = TRUE)
input_list <- lapply(on_train_files, read.csv)
on_test_files <- list.files("../3transformed_data_sep", pattern = "*on_test.csv", full.names = TRUE, recursive = TRUE)
test_list <- lapply(on_test_files, read.csv)
mat_names <- gsub("../3transformed_data_sep/|/.*$", "", on_train_files) #keep only the mid dir name
names(input_list) <- mat_names #rename the list
names(test_list) <- mat_names

#### FILTER #####
#read filter lists
on_filter_files <- list.files("../4feature_selection_sep/filter_wilc", pattern = "*on_vol.csv", full.names = TRUE, recursive = TRUE)
on_filter.list <- lapply(on_filter_files, function(x) read.csv(x)$feature) #read the files and subset feature names
names(on_filter.list) <- mat_names#rename the list

#### EMBEDDED ####
#read embedded lists
on_embedded_files <- list.files("../4feature_selection_sep/embedded_wilc/", pattern = "*on_sigVar_pass_list.csv", full.names = TRUE, recursive = TRUE)
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
#read the output from CV
cv_out_file <- paste("../5CV_sep/", class_name, "/", class_name, "_", f_name, "_CV.csv", sep = "")
cv_out <- read.csv(cv_out_file)

#one data, loop through feature sets and hyperparameters
mat_cv_list <- lapply(mat_names, MAT_CV)
mat_cv_df <- do.call('rbind', mat_cv_list)
mat_cv_df <- cbind('classfier' = rep(class_name, nrow(mat_cv_df)),
                   'light' = rep(f_name, nrow(mat_cv_df)),
                   mat_cv_df)
#save the CV result
out_file <- paste("./", class_name, "/", class_name, "_", f_name, "_test_wilc.csv",sep = '')
write.csv(mat_cv_df, out_file, quote = FALSE, row.names = FALSE)

################ OFF ################
#read the data
off_train_files <- list.files("../3transformed_data_sep", pattern = "*off_train.csv", full.names = TRUE, recursive = TRUE)
input_list <- lapply(off_train_files, read.csv)
off_test_files <- list.files("../3transformed_data_sep", pattern = "*off_test.csv", full.names = TRUE, recursive = TRUE)
test_list <- lapply(off_test_files, read.csv)
#keep offly the mid dir name
mat_names <- gsub("../3transformed_data_sep/|/.*$", "", off_train_files) #keep only the mid dir name
names(input_list) <- mat_names #rename the list
names(test_list) <- mat_names

#### FILTER #####
#read filter lists
off_filter_files <- list.files("../4feature_selection_sep", pattern = "*off_vol.csv", full.names = TRUE, recursive = TRUE)
off_filter.list <- lapply(off_filter_files, function(x) read.csv(x)$feature) #read the files and subset feature names
names(off_filter.list) <- mat_names#rename the list

#### EMBEDDED ####
#read embedded lists
off_embedded_files <- list.files("../4feature_selection_sep", pattern = "*off_sigVar_pass_list.csv", full.names = TRUE, recursive = TRUE)
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
#read the output from CV
cv_out_file <- paste("../5CV_sep/", class_name, "/", class_name, "_", f_name, "_CV.csv", sep = "")
cv_out <- read.csv(cv_out_file)

#offe data, loop through feature sets and hyperparameters
mat_cv_list <- lapply(mat_names, MAT_CV)
mat_cv_df <- do.call('rbind', mat_cv_list)
mat_cv_df <- cbind('classfier' = rep(class_name, nrow(mat_cv_df)),
                   'light' = rep(f_name, nrow(mat_cv_df)),
                   mat_cv_df)
#save the CV result
out_file <- paste("./", class_name, "/", class_name, "_", f_name, "_test_wilc.csv",sep = '')
write.csv(mat_cv_df, out_file, quote = FALSE, row.names = FALSE)
