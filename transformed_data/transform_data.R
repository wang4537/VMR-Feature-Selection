#long to wide transformation
transform_mydata <- function(test, method, light_con) #method = matricization method; light=
{
  #total distance 
  if(method == 'total_dist')
  {
    test$distance <- apply(test[,c(8,11,14)], 1, function(x) as.numeric(norm(as.matrix(x),type = '1'))) #inadist, smldist, lardist
    #subset only the uid, time(end) and distance
    test.sub <- test[,c(17,5,18)]
  }
  #total duration
  if(method == 'total_dur')
  {
    test$duration <- apply(test[,c(7,10,13)], 1, function(x) as.numeric(norm(as.matrix(x),type = '1')))#inadur, smldur, lardur
    #subset only the uid, time(end) and duration
    test.sub <- test[,c(17,5,18)]
  }
  #L1 Norm
  if(method == 'L1-Norm')
  {
    test.list <- split(test, f = as.factor(test$start)) #split by second
    test.list <- lapply(test.list, function(x)
    {
      x <- cbind(x[,c(1:5, 15:17)], scale(x[,c(6:14)])) #info columns and numeric columns
      x$L1Norm <- apply(x[,9:17], 1, function(x) as.numeric(norm(as.matrix(x),type = '1'))) #l1-norm of a vector
      return(x)
    })  
    #subset only the uid, time(end) and L1-norm
    test.sub <- do.call('rbind', test.list)[,c(8,5,18)]
  }
  if(method == 'L2-Norm')
  {
    test.list <- split(test, f = as.factor(test$start)) #split by second
    test.list <- lapply(test.list, function(x)
    {
      x <- cbind(x[,c(1:5, 15:17)], scale(x[,c(6:14)])) #info columns and numeric columns
      x$L2Norm <- apply(x[,9:17], 1, function(x) as.numeric(norm(as.matrix(x),type = '2'))) #l1-norm of a vector
      return(x)
    })  
    #subset only the uid, time(end) and L2-norm
    test.sub <- do.call('rbind', test.list)[,c(8,5,18)]
  }
  #concatenated distance
  if(method == 'concat_dist')
  {
    #subset only the uid, time(end) and all distance variables
    test.sub <- test[,c(8,11,14,5,17)]
  }
  #concatenated duration
  if(method == 'concat_dur')
  {
    #subset only the uid, time(end) and all duration variables
    test.sub <- test[,c(7,10,13,5,17)]
  }
  #concatenated duration
  if(method == 'concat_9')
  {
    #subset only the uid, time(end) and all variables
    test.sub <- test[,c(6:14,5,17)]
  }
  
  #reshape the long to wide format
  test.sub.df <- reshape(test.sub, idvar = "uid", timevar = "end", direction = "wide")
  #make the variables legal
  names(test.sub.df) <- make.names(names(test.sub.df))
  
  #specify the out dir
  out_file <- paste("./", method, "/", method, "_Opt_", light_con, ".csv", sep = "") 
  #write the data to the corresponding dir
  write.csv(test.sub.df, out_file, quote = FALSE, row.names = FALSE)
  
  return(test.sub.df)
}

#main
library(ggplot2)
library(ggrepel)
setwd("/depot/yleung/data/Feature_selection_framework/transformed_data/")

#read the data
lightOn <- read.csv("../input_data/Opt_on.csv")
lightOff <- read.csv("../input_data/Opt_off.csv")
#subset the right time window
lightOn <- lightOn[lightOn$start %in% 0:300,]
lightOff <- lightOff[lightOff$start %in% 0:300,]
#transform the data
method_list <- c('total_dist', 'total_dur', 'L1-Norm', 'L2-Norm', 'concat_dist', 'concat_dur', 'concat_9')
#create dirs for output
sapply(method_list, function(x) dir.create(paste("./",x, "/", sep = ""))) #create folds for outputs

on_tf_list <- lapply(method_list, function(x) transform_mydata(lightOn, method = x, light_con = unique(lightOn$light)))
off_tf_list <- lapply(method_list, function(x) transform_mydata(lightOff, method = x, light_con = unique(lightOff$light)))
