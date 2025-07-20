#linear
SVM_linear_acc <- function(train, valid, par_list)
{
  #build the SVM model
  classifier <- svm(formula = genotype~.,
                    data = train,
                    type = "C-classification",
                    kernel = "linear")

  #calculate the accracy
  valid_pred <- predict(classifier, valid[,-1])
  #calculate the confusion matrix
  cm <- confusionMatrix(valid_pred, as.factor(valid[,1]), positive = "2")
  #get performance metrics
  auroc <- auc(valid[,1], as.numeric(valid_pred)) #area under ROC curve
  return(auroc)
}

#polynomial
SVM_poly_acc <- function(train, valid, par_df)
{
  auc_score <- apply(par_df, 1, function(i)
  {
    #build the SVM model
    classifier <- svm(formula = genotype~.,
                      data = train,
                      type = "C-classification",
                      kernel = i[1],
                      gamma = i[2],
                      degree = i[3],
                      coef0 = i[4])
    #calculate the accuracy
    valid_pred <- predict(classifier, valid[,-1])
    #calculate the confusion matrix
    cm <- confusionMatrix(valid_pred, as.factor(valid[,1]), positive = "2")
    #get performance metrics
    auroc <- auc(valid[,1], as.numeric(valid_pred)) #area under ROC curve
    return(auroc)
  })
  return(auc_score)
}

#radial
SVM_rad_acc <- function(train, valid, par_df)
{
  auc_score <- apply(par_df, 1, function(i)
  {
    #build the SVM model
    classifier <- svm(formula = genotype~.,
                      data = train,
                      type = "C-classification",
                      kernel = i[1],
                      gamma = i[2])
    #calculate the accracy
    valid_pred <- predict(classifier, valid[,-1])
    #calculate the confusion matrix
    cm <- confusionMatrix(valid_pred, as.factor(valid[,1]), positive = "2")
    #get performance metrics
    auroc <- auc(valid[,1], as.numeric(valid_pred)) #area under ROC curve
    return(auroc)
  })
  return(auc_score)
}

#sigmoid
SVM_sig_acc <- function(train, valid, par_df)
{
  auc_score <- apply(par_df, 1, function(i)
  {
    #build the SVM model
    classifier <- svm(formula = genotype~.,
                      data = train,
                      type = "C-classification",
                      kernel = i[1],
                      gamma = i[2],
                      coef0 = i[3])
    #calculate the accracy
    valid_pred <- predict(classifier, valid[,-1])
    #calculate the confusion matrix
    cm <- confusionMatrix(valid_pred, as.factor(valid[,1]), positive = "2")
    #get performance metrics
    auroc <- auc(valid[,1], as.numeric(valid_pred)) #area under ROC curve
    return(auroc)
  })
  return(auc_score)
}

CV_SVM <- function(kernel_names)
{
  #### linear ####
  if(kernel_names == 'linear')
  {
    cv_list <- lapply(fold_uid_list, function(i)
    {
      #split the traning and testing set
      train_fold <- input[!input$uid %in% c(i,out_test),-1] #remove the uid column
      valid_fold <- input[input$uid %in% i,-1] #remove the uid column
      #score per fold
      score_df <- as.data.frame(SVM_linear_acc(train_fold, valid_fold, SVM_parameters[[kernel_names]]))
      return(score_df)
    })
    #average score across folds
    ave_CV_df <- cbind(SVM_parameters[[kernel_names]], 'AUROC' = as.numeric(rowMeans(as.data.frame(cv_list)))) #add HPs to the df
    #highest average score and best hyperparameter
    ave_CV_df_max <- ave_CV_df[which.max(ave_CV_df[,'AUROC'])[1],] #take the first row with the maximum ave AUR
    return(ave_CV_df_max)
  }
  #### poly ####
  if(kernel_names == 'polynomial')
  {
    cv_list <- lapply(fold_uid_list, function(i)
    {
      #split the traning and testing set
      train_fold <- input[!input$uid %in% c(i,out_test),-1] #remove the uid column
      valid_fold <- input[input$uid %in% i,-1] #remove the uid column
      #score per fold
      score_df <- as.data.frame(SVM_poly_acc(train_fold, valid_fold, SVM_parameters[[kernel_names]]))
      return(score_df)
    })
    #average score across folds
    ave_CV_df <- cbind(SVM_parameters[[kernel_names]], 'AUROC' = as.numeric(rowMeans(as.data.frame(cv_list)))) #add HPs to the df
    #highest average score and best hyperparameter
    ave_CV_df_max <- ave_CV_df[which.max(ave_CV_df[,'AUROC'])[1],] #take the first row with the maximum ave AUR
    return(ave_CV_df_max)
  }
  #### rad ####
  if(kernel_names == "radial")
  {
    cv_list <- lapply(fold_uid_list, function(i)
    {
      #split the traning and testing set
      train_fold <- input[!input$uid %in% c(i,out_test),-1] #remove the uid column
      valid_fold <- input[input$uid %in% i,-1] #remove the uid column
      #score per fold
      score_df <- as.data.frame(SVM_rad_acc(train_fold, valid_fold, SVM_parameters[[kernel_names]]))
      return(score_df)
    })
    #average score across folds
    ave_CV_df <- cbind(SVM_parameters[[kernel_names]], 'AUROC' = as.numeric(rowMeans(as.data.frame(cv_list)))) #add HPs to the df
    #highest average score and best hyperparameter
    ave_CV_df_max <- ave_CV_df[which.max(ave_CV_df[,'AUROC'])[1],] #take the first row with the maximum ave AUR
    return(ave_CV_df_max)
  }
  #### sig ####
  if(kernel_names == "sigmoid")
  {
    cv_list <- lapply(fold_uid_list, function(i)
    {
      #split the traning and testing set
      train_fold <- input[!input$uid %in% c(i,out_test),-1]
      valid_fold <- input[input$uid %in% i,-1]
      #score per fold
      score_df <- as.data.frame(SVM_sig_acc(train_fold, valid_fold, SVM_parameters[[kernel_names]]))
      return(score_df)
    })
    #average score across folds
    ave_CV_df <- cbind(SVM_parameters[[kernel_names]], 'AUROC' = as.numeric(rowMeans(as.data.frame(cv_list)))) #add HPs to the df
    #highest average score and best hyperparameter
    ave_CV_df_max <- ave_CV_df[which.max(ave_CV_df[,'AUROC'])[1],] #take the first row with the maximum ave AUR
    return(ave_CV_df_max)
  }

  #aggregate results
  ave_cv_feature.df <- do.call("rbind", cv_list) #bind the dfs by rows
  ave_cv_feature.df <- cbind('mat' = rep("CP", nrow(ave_cv_feature.df)), #matricization method
                             'feature' = rep("default", nrow(ave_cv_feature.df)), #feature set
                             ave_cv_feature.df) #add a column to indicate the feature set
  return(ave_cv_feature.df)
}

library(caTools)
library(e1071)
library(caret)
library(pROC)
#library(parallel)
library(doParallel)
library(foreach)# multi-core processing

setwd("/depot/yleung/data/Feature_selection_framework/8cp_classification/CV/")
class_name <- "SVM"
dir.create(class_name)
#hyper-parameters
SVM_parameters <- list("linear" = expand.grid(kernel = "linear"),
                       "polynomial" = expand.grid(kernel = "polynomial", gamma = 1:3, degree = c(2,3), coef0 = 0:5),
                       "radial" = expand.grid(kernel = "radial", gamma = 1:3),
                       "sigmoid" = expand.grid(kernel = "sigmoid", gamma = 1:3,  coef0 = 0:5))

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
#cl <- makeCluster(3)#doParallel, setup parallel processing
#registerDoParallel(cl)
SVM_names <- c("linear", "polynomial", "radial", "sigmoid")
CV_SVM_results <- lapply(SVM_names, function (z)
{
  mat_cv_df <- CV_SVM(kernel_names = z) #run CV per kernel
  mat_cv_df <- cbind('classfier' = rep(class_name, nrow(mat_cv_df)), #classifier name
                     'light' = rep(f_name, nrow(mat_cv_df)), #light condition
                     'mat' = rep("CP", nrow(mat_cv_df)), #matricization method
                     'feature' = rep("default", nrow(mat_cv_df)), #feature set
                     mat_cv_df)
  #save the CV result
  out_file <- paste("./SVM/SVM_", z, "_", f_name, "_CP_CV.csv",sep = '')
  write.csv(mat_cv_df, out_file, quote = FALSE, row.names = FALSE)

  return(mat_cv_df)
})

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
#list of test set
out_test <- unique(unlist(lapply(yellow_page[,3:ncol(yellow_page)], function(x) which(x == "Test"))))
test_uid <- yellow_page[out_test,"uid"]

#encode the genotypes
input$genotype <- as.factor(as.numeric(factor(input$genotype)))

#one data, loop through feature sets and hyperparameters
#cl <- makeCluster(3)#doParallel, setup parallel processing
#registerDoParallel(cl)
SVM_names <- c("linear", "polynomial", "radial", "sigmoid")
CV_SVM_results <- lapply(SVM_names, function (z)
{
  mat_cv_df <- CV_SVM(kernel_names = z) #run CV per kernel
  mat_cv_df <- cbind('classfier' = rep(class_name, nrow(mat_cv_df)), #classifier name
                     'light' = rep(f_name, nrow(mat_cv_df)), #light condition
                     'mat' = rep("CP", nrow(mat_cv_df)), #matricization method
                     'feature' = rep("default", nrow(mat_cv_df)), #feature set
                     mat_cv_df)
  #save the CV result
  out_file <- paste("./SVM/SVM_", z, "_", f_name, "_CP_CV.csv",sep = '')
  write.csv(mat_cv_df, out_file, quote = FALSE, row.names = FALSE)

  return(mat_cv_df)
})
