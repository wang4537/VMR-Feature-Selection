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
  #get performance metrics
  acc <- cm$overall[1] #accuracy
  sns <- cm$byClass[1] #sensitivity
  sps <- cm$byClass[2] #specificity
  prc <- cm$byClass[5] #precision
  kappa <- cm$overall[2] #Cohen's kappa
  #output a vector of metrics
  performance <- c(par_list[1],
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

#polynomial
SVM_poly_acc <- function(train, valid, j)
{
  #build the SVM model
  classifier <- svm(formula = genotype~.,
                    data = train,
                    type = "C-classification",
                    kernel = j[1],
                    gamma = j[2],
                    degree = j[3],
                    coef0 = j[4])
  #calculate the accuracy
  valid_pred <- predict(classifier, valid[,-1])
  #calculate the confusion matrix
  cm <- confusionMatrix(valid_pred, as.factor(valid[,1]), positive = "2")
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
                   j[2],
                   j[3],
                   j[4],
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

#radial
SVM_rad_acc <- function(train, valid, j)
{
  #build the SVM model
  classifier <- svm(formula = genotype~.,
                    data = train,
                    type = "C-classification",
                    kernel = j[1],
                    gamma = j[2])
  #calculate the accracy
  valid_pred <- predict(classifier, valid[,-1])
  #calculate the confusion matrix
  cm <- confusionMatrix(valid_pred, as.factor(valid[,1]), positive = "2")
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
                   j[2],
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

#sigmoid
SVM_sig_acc <- function(train, valid, j)
{
  #build the SVM model
  classifier <- svm(formula = genotype~.,
                    data = train,
                    type = "C-classification",
                    kernel = j[1],
                    gamma = j[2],
                    coef0 = j[3])
  #calculate the accracy
  valid_pred <- predict(classifier, valid[,-1])
  #calculate the confusion matrix
  cm <- confusionMatrix(valid_pred, as.factor(valid[,1]), positive = "2")
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
                   j[2],
                   j[3],
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
SVM_test <- function(kernel_name, parm) # i = matricization method name; kernel_name = kernel for SMV
{
  ## cross validation
  #linear
  if(kernel_name == 'linear')
  {
    #calculate the CV accuracy
    score_df <- as.data.frame(SVM_linear_acc(train_input, test_input, parm))
    return(score_df)
  }
  #poly
  if(kernel_name == 'polynomial')
  {
    #calculate the CV accuracy
    score_df <- as.data.frame(SVM_poly_acc(train_input, test_input, parm))
    return(score_df)
  }


  #rad
  if(kernel_name == "radial")
  {
    #calculate the CV accuracy
    score_df <- as.data.frame(SVM_rad_acc(train_input, test_input, parm))
    return(score_df)
  }


  #sig
  if(kernel_name == "sigmoid")
  {
    #calculate the CV accuracy
    score_df <- as.data.frame(SVM_sig_acc(train_input, test_input, parm))
    return(score_df)
  }
}

library(caTools)
library(e1071)
library(caret)
library(pROC)
library(parallel)
setwd("/depot/yleung/data/Feature_selection_framework/8cp_classification/final/")
class_name <- "SVM"
dir.create(class_name)

################ ON ################
f_name <- "Opt_on"
#input file name
train_file <- paste("../project_data/", f_name, "_cp_train_data.csv", sep = "")
test_file <- paste("../project_data/", f_name, "_cp_test_data.csv", sep = "")
#read the data, remoe the uid column, keep the genotype column
train_input <- read.csv(train_file)[,-1]
test_input <- read.csv(test_file)[,-1]
#encode the genotype column
train_input$genotype <- as.factor(as.numeric(factor(train_input$genotype)))
test_input$genotype <- as.factor(as.numeric(factor(test_input$genotype)))
#get the optimal parameter
SVM_names <- c("linear", "polynomial", "radial", "sigmoid") #kernel names
SVM_test.list <- lapply(SVM_names, function(s)
{
  #read the output from CV
  cv_out_file <- paste("../CV/", class_name, "/", class_name, "_", s, "_", f_name, "_CP_CV.csv", sep = "") #cv out file name
  cv_out <- read.csv(cv_out_file) #name 
  best_parm <- cv_out[,-c(1:4, ncol(cv_out)), drop = FALSE] #first 4 column = meta info, last columns = AUROC
  
  #run the SVMs
  mat_cv <- SVM_test(s, best_parm) 
  #add columns
  per_feature.df <- cbind('classfier' = rep(class_name, nrow(mat_cv)),
                          'light' = rep(f_name, nrow(mat_cv)),
                          'mat' = rep("CP", nrow(mat_cv)), #matricization method
                          'feature' = rep("default", nrow(mat_cv)), #feature set
                          mat_cv) #add a column to indicate the feature set
  
  #save the CV result
  out_file <- paste("./SVM/SVM_", s, "_", f_name, "_test_CP.csv",sep = '')
  write.csv(per_feature.df, out_file, quote = FALSE, row.names = FALSE)
  return(per_feature.df)
})



################ OFF ################
f_name <- "Opt_off"
#input file name
train_file <- paste("../project_data/", f_name, "_cp_train_data.csv", sep = "")
test_file <- paste("../project_data/", f_name, "_cp_test_data.csv", sep = "")
#read the data, remoe the uid column, keep the genotype column
train_input <- read.csv(train_file)[,-1]
test_input <- read.csv(test_file)[,-1]
#encode the genotype column
train_input$genotype <- as.factor(as.numeric(factor(train_input$genotype)))
test_input$genotype <- as.factor(as.numeric(factor(test_input$genotype)))
#get the optimal parameter
SVM_names <- c("linear", "polynomial", "radial", "sigmoid") #kernel names
SVM_test.list <- lapply(SVM_names, function(s)
{
  #read the output from CV
  cv_out_file <- paste("../CV/", class_name, "/", class_name, "_", s, "_", f_name, "_CP_CV.csv", sep = "") #cv out file name
  cv_out <- read.csv(cv_out_file) #name 
  best_parm <- cv_out[,-c(1:4, ncol(cv_out)), drop = FALSE] #first 4 column = meta info, last columns = AUROC
  
  #run the SVMs
  mat_cv <- SVM_test(s, best_parm) 
  #add columns
  per_feature.df <- cbind('classfier' = rep(class_name, nrow(mat_cv)),
                          'light' = rep(f_name, nrow(mat_cv)),
                          'mat' = rep("CP", nrow(mat_cv)), #matricization method
                          'feature' = rep("default", nrow(mat_cv)), #feature set
                          mat_cv) #add a column to indicate the feature set
  
  #save the CV result
  out_file <- paste("./SVM/SVM_", s, "_", f_name, "_test_CP.csv",sep = '')
  write.csv(per_feature.df, out_file, quote = FALSE, row.names = FALSE)
  return(per_feature.df)
})

