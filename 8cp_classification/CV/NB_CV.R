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

################ MAIN ################
library(naivebayes) # naive bayes
library(caret) # data splitting and grid search
#library(dplyr)
library(pROC) # auroc
library(doParallel)
library(foreach)# multi-core processing
setwd("/depot/yleung/data/Feature_selection_framework/8cp_classification/CV/")
#create the output directories
class_name <- 'NB'
dir.create(class_name)
#list of hyperparameters
NB_parameters <- expand.grid(usekernel = c(TRUE, FALSE), laplace = 0:3)

################# ON ################
f_name <- "Opt_on"
#input file name
input_file <- paste("../project_data/", f_name, "_cp_train_data.csv", sep = "")
#read the data
input <- read.csv(input_file)

#### Set Up Cross Validation and Yellow Page ####
yellow_page_file_name <- paste("../../2preprocess_sep/", f_name, "_kFolds_yellow_page.csv", sep = "")
yellow_page <- read.csv(yellow_page_file_name, header = TRUE)
#list of validation set per fold
fold_uid_list <- lapply(yellow_page[,3:ncol(yellow_page)], function(x) yellow_page[which(x == "Valid"), 'uid'])
#list of test set
out_test <- unique(unlist(lapply(yellow_page[,3:ncol(yellow_page)], function(x) which(x == "Test"))))
test_uid <- yellow_page[out_test,"uid"]

#encode the genotypes
input$genotype <- as.factor(as.numeric(factor(input$genotype)))

#one data, loop through feature sets and hyperparameters
cl <- makeCluster(3)#doParallel, setup parallel processing
registerDoParallel(cl)
CV_list <- foreach(i = 1:length(fold_uid_list), .packages = c("naivebayes", "pROC", "caret")) %dopar% {
  #split the traning and testing set
  train_fold <- input[!input$uid %in% c(fold_uid_list[[i]],out_test),]
  valid_fold <- input[input$uid %in% fold_uid_list[[i]],]
  #calculate the CV accuracy
  score_df <- as.data.frame(NB_k_acc(train_fold[,-1], valid_fold[-1], NB_parameters))
  return (score_df) #return the mean accuracy of 10 folds in each parameter
}

ave_CV_df <- cbind(NB_parameters, 'AUROC' = as.numeric(rowMeans(as.data.frame(CV_list)))) #add HPs to the df
ave_CV_df_max <- as.data.frame(ave_CV_df[which.max(ave_CV_df[,'AUROC'])[1],]) #take the first row with the maximum ave AUROC
ave_cv_feature.df <- cbind('mat' = "CP", #matricization method
                           'feature' = "default", #feature set
                           ave_CV_df_max) #add a column to indicate the feature set

#add addtional info the the output
mat_cv_df <- cbind('classfier' = rep(class_name, nrow(ave_cv_feature.df)),
                   'light' = rep(f_name, nrow(ave_cv_feature.df)),
                   ave_cv_feature.df)
#save the CV result
out_file <- paste("./", class_name, "/", class_name, "_", f_name, "_CV_cp.csv",sep = '')
write.csv(mat_cv_df, out_file, quote = FALSE, row.names = FALSE)

################ OFF ################
f_name <- "Opt_off"
#input file name
input_file <- paste("../project_data/", f_name, "_cp_train_data.csv", sep = "")
#read the data
input <- read.csv(input_file)

#### Set Up Cross Validation and Yellow Page ####
yellow_page_file_name <- paste("../../2preprocess_sep/", f_name, "_kFolds_yellow_page.csv", sep = "")
yellow_page <- read.csv(yellow_page_file_name, header = TRUE)
#list of validation set per fold
fold_uid_list <- lapply(yellow_page[,3:ncol(yellow_page)], function(x) yellow_page[which(x == "Valid"), 'uid'])
#list of test setf_name <- "Opt_on"
#input file name
input_file <- paste("../project_data/", f_name, "cp_train_data.csv", sep = "_")
#read the data
input <- read.csv("../project_data/Opt_on_cp_train_data.csv")

out_test <- unique(unlist(lapply(yellow_page[,3:ncol(yellow_page)], function(x) which(x == "Test"))))
test_uid <- yellow_page[out_test,"uid"]

#encode the genotypes
input$genotype <- as.factor(as.numeric(factor(input$genotype)))

#one data, loop through feature sets and hyperparameters
cl <- makeCluster(2)#doParallel, setup parallel processing
registerDoParallel(cl)
CV_list <- foreach(i = 1:length(fold_uid_list), .packages = c("naivebayes", "pROC", "caret")) %dopar% {
  #split the traning and testing set
  train_fold <- input[!input$uid %in% c(fold_uid_list[[i]],out_test),]
  valid_fold <- input[input$uid %in% fold_uid_list[[i]],]
  #calculate the CV accuracy
  score_df <- as.data.frame(NB_k_acc(train_fold[,-1], valid_fold[-1], NB_parameters))
  return (score_df) #return the mean accuracy of 10 folds in each parameter
}

ave_CV_df <- cbind(NB_parameters, 'AUROC' = as.numeric(rowMeans(as.data.frame(CV_list)))) #add HPs to the df
ave_CV_df_max <- as.data.frame(ave_CV_df[which.max(ave_CV_df[,'AUROC'])[1],]) #take the first row with the maximum ave AUROC
ave_cv_feature.df <- cbind('mat' = "CP", #matricization method
                           'feature' = "default", #feature set
                           ave_CV_df_max) #add a column to indicate the feature set

#add addtional info the the output
mat_cv_df <- cbind('classfier' = rep(class_name, nrow(ave_cv_feature.df)),
                   'light' = rep(f_name, nrow(ave_cv_feature.df)),
                   ave_cv_feature.df)
#save the CV result
out_file <- paste("./", class_name, "/", class_name, "_", f_name, "_CV_cp.csv",sep = '')
write.csv(mat_cv_df, out_file, quote = FALSE, row.names = FALSE)
