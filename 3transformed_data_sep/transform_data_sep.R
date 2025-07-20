#long to wide transformation
transform_mydata <- function(input, method, light_con, partition) #method = matricization method; light_con = on or off; partiion = train or test
{
  #total distance 
  if(method == 'total_dist')
  {
    input$distance <- apply(input[,c(8,11,14)], 1, function(x) as.numeric(norm(as.matrix(x),type = '1'))) #inadist, smldist, lardist
    #subset only the uid, time(end) and distance
    input.sub <- input[,c(17,5,18)]
  }
  
  #L1 Norm
  if(method == 'L1-Norm')
  {
    input.list <- split(input, f = as.factor(input$start)) #split by second
    input.list <- lapply(input.list, function(x)
    {
      x <- cbind(x[,c(1:5, 15:17)], scale(x[,c(6:14)])) #info columns and numeric columns
      x$L1Norm <- apply(x[,9:17], 1, function(x) as.numeric(norm(as.matrix(x),type = '1'))) #l1-norm of a vector
      return(x)
    })  
    #subset only the uid, time(end) and L1-norm
    input.sub <- do.call('rbind', input.list)[,c(8,5,18)]
  }
  if(method == 'L2-Norm')
  {
    input.list <- split(input, f = as.factor(input$start)) #split by second
    input.list <- lapply(input.list, function(x)
    {
      x <- cbind(x[,c(1:5, 15:17)], scale(x[,c(6:14)])) #info columns and numeric columns
      x$L2Norm <- apply(x[,9:17], 1, function(x) as.numeric(norm(as.matrix(x),type = '2'))) #l1-norm of a vector
      return(x)
    })  
    #subset only the uid, time(end) and L2-norm
    input.sub <- do.call('rbind', input.list)[,c(8,5,18)]
  }
  #concatenated distance
  if(method == 'concat_dist')
  {
    #subset only the uid, time(end) and all distance variables
    input.sub <- input[,c(8,11,14,5,17)]
  }
  #concatenated duration
  if(method == 'concat_9')
  {
    #subset only the uid, time(end) and all variables
    input.sub <- input[,c(6:14,5,17)]
  }
  
  #reshape the long to wide format
  input.sub.df <- reshape(input.sub, idvar = "uid", timevar = "end", direction = "wide")
  #make the variables legal
  names(input.sub.df) <- make.names(names(input.sub.df))
  
  #specify the out dir
  out_file <- paste("./", method, "/", method, "_Opt_", light_con, "_", partition,".csv", sep = "") 
  #write the data to the corresponding dir
  write.csv(input.sub.df, out_file, quote = FALSE, row.names = FALSE)
  
  return(input.sub.df)
}

#### main ####
library(ggplot2)

setwd("/depot/yleung/data/Feature_selection_framework/3transformed_data_sep/")
source("../2preprocess_sep/FUNCTION_confounding_removal.R") #read the normalization functions
#list of transformation
method_list <- c('total_dist', 'L1-Norm', 'L2-Norm', 'concat_dist', 'concat_9')
#create dirs for output
sapply(method_list, function(x) dir.create(paste("./",x, "/", sep = ""))) #create folds for outputs

### ON ###
lightOn <- read.csv("../2preprocess_sep/Opt_on.csv") #read the data
yellow_on <- read.csv("../2preprocess_sep/Opt_on_kFolds_yellow_page.csv") #read the yellow page
out_on <- unique(unlist(lapply(yellow_on[,3:ncol(yellow_on)], function(x) which(x == "Test")))) #row number for test labels
#create a list to seperate train and test
lightOn_list <- list("train" = lightOn[lightOn$uid %in% yellow_on[-out_on,'uid'], ],
                "test" = lightOn[lightOn$uid %in% yellow_on[out_on,'uid'], ])
#normalize and transform the data
transform_list_on <- list()
for (i in names(lightOn_list))
{
  norm_data <- normalization(lightOn_list[[i]]) #normalize data
  #save the normalized data
  out_name <- paste("../2preprocess_sep/Opt_on_", i,"_norm.csv", sep = "")
  write.csv(norm_data, out_name, row.names = FALSE)
  #remove the baseline period
  norm_data <- norm_data[norm_data$start >= 0, ]
  #transform data
  transform_list_on[i] <- lapply(method_list, function(x) return(transform_mydata(norm_data, x, 'on', i)))
  rm(i,norm_data)
}

### OFF ###
lightOff <- read.csv("../2preprocess_sep/Opt_off.csv") #read the data
yellow_off <- read.csv("../2preprocess_sep/Opt_off_kFolds_yellow_page.csv") #read the yellow page
out_off <- unique(unlist(lapply(yellow_off[,3:ncol(yellow_off)], function(x) which(x == "Test")))) #row number for test labels
#create a list to seperate train and test
lightOff_list <- list("train" = lightOff[lightOff$uid %in% yellow_off[-out_off,'uid'], ],
                     "test" = lightOff[lightOff$uid %in% yellow_off[out_off,'uid'], ])
#normalize and transform the data
transform_list_off <- list()
for (i in names(lightOff_list))
{
  norm_data <- normalization(lightOff_list[[i]]) #normalize data
  #save the normalized data
  out_name <- paste("../2preprocess_sep/Opt_off_", i,"_norm.csv", sep = "")
  write.csv(norm_data, out_name, row.names = FALSE)
  #remove the baseline period
  norm_data <- norm_data[norm_data$start >= 0, ]
  #transform data
  transform_list_off[i] <- lapply(method_list, function(x) return(transform_mydata(norm_data, x, 'off', i)))
  rm(i)
}