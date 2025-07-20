#cross validation function for 3nn
kNN_k_acc <- function(train, valid, i)
{
  #build the svm classifier
  classifier <- knn(train[,-1,drop = FALSE],
                    valid[,-1,drop = FALSE],
                    cl = train[,1],
                    k = i[1],
                    prob = TRUE)
  #calculate the confusion matrix
  cm <- confusionMatrix(classifier, valid[,1], positive = "2")
  #get performance metrics
  acc <- cm$overall[1] #accuracy
  sns <- cm$byClass[1] #sensitivity
  sps <- cm$byClass[2] #specificity
  prc <- cm$byClass[5] #precision
  kappa <- cm$overall[2] #Cohen's kappa
  auroc <- auc(valid[,1], as.numeric(classifier)) #area under ROC curve
  #output a vector of metrics
  performance <- c('k' = i[1],
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

################ MAIN ################
library(class)
library(caret)
library(pROC)
library(parallel)
setwd("/depot/yleung/data/Feature_selection_framework/8cp_classification/final/")
#create the output directories
class_name <- 'KNN'
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
mat_cv_df <- kNN_k_acc(train_input,test_input, best_parm)

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
mat_cv_df <- kNN_k_acc(train_input,test_input, best_parm)

#save the output
per_feature.df <- cbind('classfier' = rep(class_name, nrow(mat_cv_df)),
                        'light' = rep(f_name, nrow(mat_cv_df)),
                        'mat' = rep("CP", nrow(mat_cv_df)), #matricization method
                        'feature' = rep("default", nrow(mat_cv_df)), #feature set
                        mat_cv_df) #add a column to indicate the feature set
#save the CV result
out_file <- paste("./", class_name, "/", class_name, "_", f_name, "_CP_test.csv",sep = '')
write.csv(per_feature.df, out_file, quote = FALSE, row.names = FALSE)
