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


library(xgboost) # XGB
library(caret) # data splitting and feature selection grid
library(pROC) # auc
library(parallel) # multi-core processing
setwd("/depot/yleung/data/Feature_selection_framework/8cp_classification/final/")
#create the output directories
class_name <- 'XGB'
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
levels(train_input$genotype) <- c(0,1) #1 is WT, 0 is Q344X
train_input$genotype <- as.numeric(as.character(train_input$genotype))

test_input$genotype <- as.factor(as.numeric(factor(test_input$genotype)))
levels(test_input$genotype) <- c(0,1) #1 is WT, 0 is Q344X
test_input$genotype <- as.numeric(as.character(test_input$genotype))

#run the classifier
mat_cv_df <- XGB_acc(train_input,test_input, best_parm)

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
levels(train_input$genotype) <- c(0,1) #1 is WT, 0 is Q344X
train_input$genotype <- as.numeric(as.character(train_input$genotype))

test_input$genotype <- as.factor(as.numeric(factor(test_input$genotype)))
levels(test_input$genotype) <- c(0,1) #1 is WT, 0 is Q344X
test_input$genotype <- as.numeric(as.character(test_input$genotype))

#run the classifier
mat_cv_df <- XGB_acc(train_input,test_input, best_parm)

#save the output
per_feature.df <- cbind('classfier' = rep(class_name, nrow(mat_cv_df)),
                        'light' = rep(f_name, nrow(mat_cv_df)),
                        'mat' = rep("CP", nrow(mat_cv_df)), #matricization method
                        'feature' = rep("default", nrow(mat_cv_df)), #feature set
                        mat_cv_df) #add a column to indicate the feature set
#save the CV result
out_file <- paste("./", class_name, "/", class_name, "_", f_name, "_CP_test.csv",sep = '')
write.csv(per_feature.df, out_file, quote = FALSE, row.names = FALSE)
