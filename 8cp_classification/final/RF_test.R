#cross validation function for 3nn
RF_k_acc <- function(train, valid, i)
{
  #build the svm classifier
  classifier <- randomForest(genotype~.,
                             data = train,
                             proximity = TRUE,
                             ntree = as.numeric(i[1]),
                             mtry = as.numeric(i[2]))

  #calculate the accracy
  valid_pred <- predict(classifier, valid[,-1,drop = FALSE])
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
  performance <- c(i[1],
                   i[2],
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

library(randomForest) # random forest
library(caret) # data splitting and grid search
#library(dplyr)
library(pROC) # auroc
library(parallel) # multi-core processing
setwd("/depot/yleung/data/Feature_selection_framework/8cp_classification/final/")
#create the output directories
class_name <- 'RF'
dir.create(class_name)

################ ON ################
f_name <- "Opt_on"
#input file name
train_file <- paste("../project_data/", f_name, "_cp_train_data.csv", sep = "")
test_file <- paste("../project_data/", f_name, "_cp_test_data.csv", sep = "")
#read the data, remoe the uid column, keep the genotype column
train_input <- read.csv(train_file)[,-1]
test_input <- read.csv(test_file)[,-1]
#get the optimal parameter
cv_out_file <- paste("../CV/", class_name, "/", class_name, "_", f_name, "_CV_cp.csv", sep = "")
cv_out <- read.csv(cv_out_file)
best_parm <- cv_out[,-c(1:4, ncol(cv_out)), drop = FALSE] #first 4 column = meta info, last columns = AUROC

#encode the genotype column
train_input$genotype <- as.factor(as.numeric(factor(train_input$genotype)))
test_input$genotype <- as.factor(as.numeric(factor(test_input$genotype)))

#run the classifier
mat_cv_df <- RF_k_acc(train_input,test_input, best_parm)

#save the output
per_feature.df <- cbind('classfier' = rep(class_name, nrow(mat_cv_df)),
                        'light' = rep(f_name, nrow(mat_cv_df)),
                        'mat' = rep("CP", nrow(mat_cv_df)), #matricization method
                        'feature' = rep("default", nrow(mat_cv_df)), #feature set
                        mat_cv_df) #add a column to indicate the feature set
#save the CV result
out_file <- paste("./", class_name, "/", class_name, "_", f_name, "_CP_test.csv",sep = '')
write.csv(per_feature.df, out_file, quote = FALSE, row.names = FALSE)

################ OFF ################
f_name <- "Opt_off"
#input file name
train_file <- paste("../project_data/", f_name, "_cp_train_data.csv", sep = "")
test_file <- paste("../project_data/", f_name, "_cp_test_data.csv", sep = "")
#read the data, remoe the uid column, keep the genotype column
train_input <- read.csv(train_file)[,-1]
test_input <- read.csv(test_file)[,-1]
#get the optimal parameter
cv_out_file <- paste("../CV/", class_name, "/", class_name, "_", f_name, "_CV_cp.csv", sep = "")
cv_out <- read.csv(cv_out_file)
best_parm <- cv_out[,-c(1:4, ncol(cv_out)), drop = FALSE] #first 4 column = meta info, last columns = AUROC

#encode the genotype column
train_input$genotype <- as.factor(as.numeric(factor(train_input$genotype)))
test_input$genotype <- as.factor(as.numeric(factor(test_input$genotype)))

#run the classifier
mat_cv_df <- RF_k_acc(train_input,test_input, best_parm)

#save the output
per_feature.df <- cbind('classfier' = rep(class_name, nrow(mat_cv_df)),
                        'light' = rep(f_name, nrow(mat_cv_df)),
                        'mat' = rep("CP", nrow(mat_cv_df)), #matricization method
                        'feature' = rep("default", nrow(mat_cv_df)), #feature set
                        mat_cv_df) #add a column to indicate the feature set
#save the CV result
out_file <- paste("./", class_name, "/", class_name, "_", f_name, "_CP_test.csv",sep = '')
write.csv(per_feature.df, out_file, quote = FALSE, row.names = FALSE)
