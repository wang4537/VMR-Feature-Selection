#### Main ####
library(readxl) #read excel
library(caret) #data partition
library(pbapply) #progress indicator
library(parallel) #parallel processing
setwd("/depot/yleung/data/Feature_selection_framework/2preprocess_sep/")
#get f_names for on and off data
f_name <- c("Opt_on", "Opt_off")

for(i in f_name)
{
  yellow_page_file_name <- paste("./", i, "_kFolds_yellow_page.csv", sep = "") #create yellow page file name
  #read the input data
  input_file <- paste("./", i, ".csv", sep = "")
  input <- read.csv(input_file)
  #create the yellow page
  yellow_page <- cbind(input$uid,input$genotype)
  colnames(yellow_page) <- c("uid", "genotype")
  
  #keep the unique rows
  yellow_page <- yellow_page[!duplicated(yellow_page), ]
  #get the row numbers for each row
  set.seed(222)
  out_test <- unlist(createDataPartition(yellow_page[,2], times = 1, p = 0.2))
  set.seed(333)
  fold.list <- createFolds(yellow_page[-out_test,2], k = 10) # create the folds based the the #mark the outer testing set
  
  #mark each sample with train or test, based on folds
  yellow_page_fold <- lapply(fold.list, function(x)
  {
    y <- ifelse(1:nrow(yellow_page) %in% out_test,
                "Test", #if in out_test, then "outer test set"
                ifelse(1:nrow(yellow_page) %in% x,
                       "Valid", #if in fold list, then "inner test"
                       "Train"))
    return(y)
  })
  #join all markings by column
  yellow_page_fold.df <- do.call("cbind", yellow_page_fold)
  #join the identifiers with the markings
  yellow_page.full <- cbind(yellow_page, yellow_page_fold.df) 
  #save the yellow page to a csv file
  write.csv(yellow_page.full, yellow_page_file_name, row.names = FALSE, quote = FALSE, col.names = TRUE)
}