#function: RF and tuneRF #
run_RF_tuneRF <- function(TRAIN, LIGHT, F_NAME)
{
  #### RF ####
  TRAIN[,1] <- as.factor(TRAIN[,1])

  #### tuneRF ####
  set.seed(222)
  #mtry refers to the number of variables
  tune_out <- tuneRF(TRAIN[,-1],
              TRAIN[,1],
              stepFactor =  0.5,
              plot = TRUE,
              ntreeTry = 500,
              trace = TRUE,
              improve = 0.05)
  print(tune_out)

  #rebuild the model using the optimal mtry
  mtry.opt <- tune_out[which.min(tune_out[,2]),1]
  set.seed(222)
  rf.opt <- randomForest(uid~., data = TRAIN, proximity = TRUE, mtry = mtry.opt) #default mtry  = sqrt(p), ntree = 500

  #importatn variables after optimization
  #importance of the variables
  var.imp <- importance(rf.opt)
  cut_off <- ifelse(LIGHT == "on",
                    on_cutoff[on_cutoff$mat == F_NAME,2],
                    off_cutoff[off_cutoff$mat == F_NAME,2])
  sig <- var.imp[var.imp > cut_off,,drop=FALSE]

  #output the full list
  out_dir <- paste("./", F_NAME, "/", sep = "")
  if(!dir.exists(out_dir)){dir.create(out_dir)}
  out_file <- paste(out_dir, F_NAME, "_Opt_", LIGHT, "_RF_Gini_list.csv", sep = "")
  write.csv(var.imp, out_file)

  #output a list of the sig vars passing the threshold
  out_file <- paste(out_dir, F_NAME, "_Opt_", LIGHT, "_sigVar_pass_list.csv", sep = "")
  write.csv(sig, out_file)

  #output the plot
  out_file <- paste(F_NAME, "_Opt_", LIGHT, "_RF_model_sig_var.pdf", sep = "")
  pdf(paste(out_dir, out_file, sep = ""))
  varImpPlot(rf.opt, sort = TRUE, n.var = 30, cex.lab=2.4, cex.axis=1.8, cex.main=2.4, cex.sub=1.5, main = "First 30 Features")
  abline(v = cut_off, col = "black", lwd = 1.5, lty = 2)
  dev.off()

  return(rf.opt)
}

#Main function
main <- function(x, light, MAT_NAME) #x = input file name;light = light on or off;MAT_NAME = matricization method
{
  #### Set Up Cross Validation and Yellow Page ####
  #create the uid by batch number, location and genotype
  yellow_page_file_name <- paste("../int_output/Opt_", light, "_kFolds_yellow_page.csv", sep = "")
  if(file.exists(yellow_page_file_name))
  {
    yellow_page <- read.csv(yellow_page_file_name, header = TRUE)
    out_test <- unlist(lapply(yellow_page[,3:ncol(yellow_page)], function(x) which(x == "Test")))
  }else
  {
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
    {y <- ifelse(1:nrow(yellow_page) %in% out_test,
                 "Test", #if in out_test, then "outer test set"
                 ifelse(1:nrow(yellow_page) %in% x,
                        "Valid", #if in fold list, then "inner test"
                        "Train")
    )
    return(y)
    })
    #join all markings by column
    yellow_page_fold.df <- do.call("cbind", yellow_page_fold)
    #join the identifiers with the markings
    yellow_page.full <- cbind(yellow_page, yellow_page_fold.df)
    #save the yellow page to a csv file
    write.csv(yellow_page.full, yellow_page_file_name, row.names = FALSE, quote = FALSE, col.names = TRUE)
  }

  #### Prepare Train and Test
  #split the training and testing set
  x$uid <- as.factor(as.numeric(factor(yellow_page[,2])))
  levels(x$uid) <- c(0,1) #1 is WT, 0 is Q344X
  x_input <- x[-out_test,] #use the training data

  #run the RF function
  best_RF <- run_RF_tuneRF(x_input,  LIGHT = light, F_NAME = MAT_NAME)

  return(best_RF)
}

#load the library
library(aod)
library(broom)
library(caret)
library(randomForest)
#set the working directory
setwd("/depot/yleung/data/Feature_selection_framework/embedded_wilc/")
#read the data
on_files <- list.files("../transformed_data", pattern = "*on.csv", full.names = TRUE, recursive = TRUE)
mat_names <- gsub("../transformed_data/|/.*$", "", on_files) #keep only the mid dir name
on_list <- lapply(on_files, read.csv) #read the files
names(on_list) <- mat_names #rename the list
#create the output directories
sapply(mat_names, function(x) dir.create(paste("./",x, "/", sep = ""))) #create folds for outputs
#cutoff for mean decrease gini
on_cutoff <- data.frame("mat" = mat_names,
                        "cutoff" = c(0.65, 0.85, 0.75, 2, 1.5, 5, 3))
#run the RF wrapper
on_models <- list()
for(i in 1:length(on_list))
{
  on_models[[i]] <- main(on_list[[i]], "on", names(on_list)[i])
}
names(on_models) <- names(on_list)

#read the off data
off_files <- list.files("../transformed_data", pattern = "*off.csv", full.names = TRUE, recursive = TRUE)
mat_names <- gsub("../transformed_data/|/.*$", "", on_files) #keep only the mid dir name
off_list <- lapply(off_files, read.csv) #read the files
names(off_list) <- mat_names
#cutoff for mean decrease gini
off_cutoff <- data.frame("mat" = mat_names,
                         "cutoff" = c(0.8, 0.75, 1.6, 1.2, 1.5, 2, 2))
#run the RF wrapper
off_models <- list()
for(i in 1:length(off_list))
{
  off_models[[i]] <- main(off_list[[i]], "off", names(off_list)[i])
}
names(off_models) <- names(off_list)
