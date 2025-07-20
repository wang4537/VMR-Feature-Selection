#require library to read xlsx files
library(readxl) #read excel
library(dplyr) #pipe

#### light intensity normalization ####
light_norm <- function(x) #x = input data
{
  #read the file 
  intensity_file <- "/depot/yleung/data/Feature_selection_framework/2preprocess_sep/96 Well Intensities.xlsx"
  light_intensity <- as.data.frame(read_xlsx(intensity_file,
                                             sheet = "Sheet1",
                                             range = "A2:C97",
                                             col_names = c("Location", "Well", "Intensity")))
  
  #get the light intensity for corresponding well
  x$lt_int <- light_intensity[match(x$location, light_intensity$Location),3]
  
  #wrap a function of linear regression
  lm_func <- function (y, data) lm(y ~ lt_int, data)
  #conduct linear regression to the for each of the parameters
  lm.list <- pblapply(x[,6:14], lm_func, x)
  #get the summary of the linear regression result
  lm.sum <- lapply(lm.list, summary)
  #make a data frame with the residuals from each list
  lm.residuals <- list()
  lm.fitted_value <- list()
  for(i in 1:length(lm.sum))
  {
    if(lm.sum[[i]]$coefficients[2,4] < 0.05) #when the coeffficient is signficant
    {
      lm.fitted_value[[i]] <- lm.list[[i]]$fitted.values
      lm.residuals[[i]] <- lm.list[[i]]$residuals
    }
    else #when the coeffcient is not signficant
    {
      lm.fitted_value[[i]] <- rep(lm.sum[[i]]$coefficients[1,1], nrow(x)) #use only the intercept for coefficient
      lm.residuals[[i]] <- x[,i+4] - lm.sum[[i]]$coefficients[1,1] #residuals = actual value - intercept
      
    }
  }
  #transform the list to a data frame
  lm.residuals <- as.data.frame(lm.residuals)
  colnames(lm.residuals) <- colnames(x)[6:14] #renmae the columns
  names(lm.fitted_value) <- colnames(x)[6:14] #rename the list names
  #get off set values
  mul <- as.data.frame(lapply(lm.residuals, function(y) abs(min(y))))
  lm.nm <- apply(lm.residuals, 1, function(m) m+mul)
  lm.nm <- do.call("rbind", lm.nm)
  #join two data frames to generate light intensity normalized output.sub
  x_lt_norm <- cbind(x[,1:5], lm.nm, x[,15:17])
  
  return(x_lt_norm)
}

#### Batch Effect Normalization ####
batch_norm <- function(x) #x = input data
{
  ## batch normalization ##
  #batch normalization function
  lm_func <- function(Y,X) lm(Y ~ X) # X = bahavior, Y = Batch
  #run the linear fit through every parameters
  lm.list <- pblapply(x[,6:14], lm_func, as.factor(x$batch))
  #output the summary for the linear regression result
  lm.sum <- pblapply(lm.list, summary)
  #make a data frame with the residuals from each list
  lm.residuals <- list()
  for(j in 1:length(lm.sum))
  {
    lm.residuals[[j]] <- lm.list[[j]]$residuals
  }
  #transform the list to a data frame
  lm.residuals <- as.data.frame(lm.residuals)
  colnames(lm.residuals) <- colnames(x)[6:14] #rename the columns
  #get off set values
  mul <- t(apply(lm.residuals, 2, function(y) abs(min(y))))
  lm.nm <- t(apply(lm.residuals, 1, function(m) m+mul))
  colnames(lm.nm) <- colnames(x)[6:14]
  #join two data frames to generate normalized output.sub
  x_batch_norm <- cbind(x[,1:5], lm.nm, x[,15:17])
  
  return(x_batch_norm)
}

#### Baseline Normalization ####
baseline_norm <- function(x) #x = input data
{
  #define basline period based on data set
  baseline_T <- -30:-1
  #parse the baseline period
  dk_adpt <- x[x$start %in% baseline_T,]
  #take the mean for all measurement
  beta <- colMeans(dk_adpt[,6:14])
  #split the dk_adpt to each of the genotypes
  dk_adpt.list <- list()
  group.name <- unique(dk_adpt$genotype)
  #assign genotype names as the dataframe name, i.e. parse the genotype dataframes
  for(i in 1:length(group.name))
  {
    dk_adpt.list[[i]] <- dk_adpt[dk_adpt$genotype == group.name[i],]
    names(dk_adpt.list)[i] <- as.character(group.name[i])
  }
  #calculate the mean for each genotype
  dk_adpt.mean <- pblapply(dk_adpt.list, function(y) colMeans(y[,6:14]))
  #calculate predicated baseline activity
  dk_adpt.prd <- pblapply(dk_adpt.mean, function(y) y-beta)
  #convert it to dataframe
  dk_adpt.prd.df <- data.frame(matrix(unlist(dk_adpt.prd), nrow = length(dk_adpt.prd), byrow = TRUE))
  rownames(dk_adpt.prd.df) <- names(dk_adpt.prd)
  colnames(dk_adpt.prd.df) <- names(beta)
  #set offset values, min for each column
  mub <- pbapply(dk_adpt.prd.df, 2, function(y) abs(min(y)))
  
  treatment_list <- list()
  group_type <- unique(x$genotype)
  for(i in 1:length(group_type))
  {
    treatment_list[[i]] <- x[x$genotype == group_type[i],]
    names(treatment_list)[i] <- as.character(group_type[i])
  }
  #make an empty list to store normlalized data
  treatment_list.blNorm <- list()
  #loop through each genotype
  for(i in 1:length(treatment_list))
  {
    #normalized activity = actual activity - predicated activity + offset
    #loop through each row of data frame in the list
    genorm <- pbapply(treatment_list[[i]][,6:14], 1, function(x) x- dk_adpt.prd.df[i,]+mub+0.0322) #offset to be non-negative
    #transform the list to a data frame
    genorm <- data.frame(matrix(unlist(genorm),
                                nrow = length(genorm),
                                byrow = TRUE))
    #join two data frames
    genorm <- cbind(treatment_list[[i]][,1:5], genorm, treatment_list[[i]][,15:17])
    #add the normalized data frame to the list
    treatment_list.blNorm[[i]] <- genorm
    #change the columns to the measurement names
    colnames(treatment_list.blNorm[[i]])[6:14] <- colnames(treatment_list[[i]])[6:14]
    rm(genorm)
  }
  #change the names of the list to the genotype
  names(treatment_list.blNorm) <- names(treatment_list)
  #save the baseline normalized data
  bl_df <- do.call("rbind", treatment_list.blNorm)
  
  return(bl_df)
}

#### Master Wrap ####
normalization <- function(z)#z = master data
{
  return(light_norm(z) %>% batch_norm %>% baseline_norm)
}



