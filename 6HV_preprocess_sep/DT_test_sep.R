#cross validation function for 3nn
DT_k_acc <- function(train, valid, j)
{
  #DT
  classifier <- rpart(uid~.,
                      train,
                      method = "class",
                      parms = list(split = j))
  #calculate the accracy
  valid_pred <- predict(classifier, valid[,-1, drop=FALSE], type = 'class')
  #calculate the confusion matrix
  cm <- confusionMatrix(valid_pred, valid[,1], positive = "2")
  #get performance metrics
  auroc <- auc(valid[,1], as.numeric(valid_pred)) #area under ROC curve
  #get performance metrics
  acc <- cm$overall[1] #accuracy
  sns <- cm$byClass[1] #sensitivity
  sps <- cm$byClass[2] #specificity
  prc <- cm$byClass[5] #precision
  kappa <- cm$overall[2] #Cohen's kappa
  #output a vector of metrics
  performance <- c(j[1],
                   "AUROC" = auroc,
                   acc,
                   sns,
                   sps,
                   prc,
                   kappa)

  #calculate the CV accuracy
  score_df <- t(data.frame(performance))

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
  
  #subset the cv output by mat name
  cv_out_sub <- cv_out[cv_out$mat == i,] #GLOBAL VAR: cv_out

  ## cross validation
  per_feature.list <- lapply(names(feature_set), function(y) #loop through feature sets
  {
    train_fold <- input_train[,c("uid",feature_set[[y]])]# parse the feature set
    test_fold <- input_test[,c("uid",feature_set[[y]])]
    #hyperparameters
    DT_parameters <- cv_out_sub[cv_out_sub$feature == y , 5]

    #calculate performances
    score_df <- as.data.frame(DT_k_acc(train_fold, test_fold, DT_parameters))
    return (score_df) #return the mean accuracy of 10 folds in each parameter
  })

  names(per_feature.list) <- names(feature_set)#rename the list
  per_feature.df <- do.call("rbind", per_feature.list) #bind the dfs by rows
  per_feature.df <- cbind('mat' = rep(i, nrow(per_feature.df)), #matricization method
                          'feature' = row.names(per_feature.df), #feature set
                          per_feature.df) #add a column to indicate the feature set
  return(per_feature.df)
}

library(rpart) # decision tree
library(caret) # data splitting and grid search
#library(dplyr)
library(pROC) # auroc
library(parallel) # multi-core processing
setwd("/depot/yleung/data/Feature_selection_framework/6HV_preprocess_sep/")
#create the output directories
class_name <- 'DT'
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
rm(i)
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
rm(i)
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
