Prep_LightNorm<- function(x) #x, filename of the data file
{
  ## preprocessing ##
  output <- read.csv(x, header = TRUE, sep = "\t") 
  out.columns <- colnames(output)[c(1,4,5,7,8,12,13,14,15,16,17,18,19,20)]#subset the output dataframe by the columns
  output.sub <- as.data.frame(output[out.columns])
  output.sub <- output.sub[output.sub$an == 2,]
  #keep only the desired time segment
  time_segment <- (3600-30):3899
  output.sub <- output.sub[output.sub$start %in% time_segment, ] 
  
  ## light normalization ##
  output.sub$lt_int <- rep(light_intensity$Intensity, nrow(output.sub)/96)
  #keep only the WT and Q344X
  output.sub <- output.sub[output.sub$genotype %in% c("WT", "Q344X", "WT+DMSO", "Q344X+DMSO"), ]
  #wrap a function of linear regression
  lm_func <- function (y, data) lm(y ~ lt_int, data)
  #conduct linear regression to the for each of the parameters
  lm.list <- pblapply(output.sub[,6:14], lm_func, output.sub)
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
      lm.fitted_value[[i]] <- rep(lm.sum[[i]]$coefficients[1,1], nrow(output.sub)) #use only the intercept for coefficient
      lm.residuals[[i]] <- output.sub[,i+4] - lm.sum[[i]]$coefficients[1,1] #residuals = actual value - intercept
      
    }
  }
  #transform the list to a data frame
  lm.residuals <- as.data.frame(lm.residuals)
  colnames(lm.residuals) <- colnames(output.sub)[6:14]
  names(lm.fitted_value) <- colnames(output.sub)[6:14]
  #get off set values
  mul <- as.data.frame(lapply(lm.residuals, function(y) abs(min(y))))
  lm.nm <- apply(lm.residuals, 1, function(x) x+mul)
  lm.nm <- do.call("rbind", lm.nm)
  #join two data frames to generate light intensity normalized output.sub
  output.sub_lt_norm <- cbind(output.sub[,1:5], lm.nm)
}

#### Main ####
library(readxl) #read excel
library(pbapply) #progress indicator
library(parallel) #parallel processing
setwd("/depot/yleung/data/Feature_selection_framework/preprocessing/")
#read the light intensity
light_intensity <- read_xlsx("96 Well Intensities.xlsx",
                             sheet = "Sheet1",
                             range = "A2:B97",
                             col_names = c("Well", "Intensity"))
#read the training data
data_files <- c(list.files("../../beichen temp/optimized DMSO Q344X/light-on/", full.names = TRUE),
                list.files("../../beichen temp/optimized untreated Q344X/light-on/", full.names = TRUE))

## preprocess and light normalization ##
training_list <- mclapply(data_files, Prep_LightNorm, mc.cores = 8)
#add the batch number 
for(i in 1:length(training_list))
{
  training_list[[i]]$batch = i
}
#save the intermediate output
saveRDS(training_list, "int_output/light_normalized_List_opt_all_lightOn.Rds")

## batch normalization ##
training_df <- do.call("rbind", training_list) #rbind all dfs in the list
rm(training_list)
#batch normalization function
lm_func <- function(y,x) lm(y ~ x)
#run the linear fit through every parameters
lm.list <- pblapply(training_df[,6:14], lm_func, as.factor(training_df$batch))
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
colnames(lm.residuals) <- colnames(training_df)[6:14]
#get off set values
mul <- t(apply(lm.residuals, 2, function(y) abs(min(y))))
lm.nm <- t(apply(lm.residuals, 1, function(x) x+mul))
colnames(lm.nm) <- colnames(training_df)[6:14]
#join two data frames to generate normalized output.sub
training_df.norm <- cbind(training_df[,1:5], lm.nm, training_df$batch)
colnames(training_df.norm)[ncol(training_df.norm)] <- "batch"
rm(lm_func, lm.list, lm.sum, lm.residuals, mul, lm.nm)
#save the batch corrrected data
saveRDS(training_df.norm, "int_output/batch_corrected_opt_Untreated_all_lightOn.Rds")
rm(training_df)

training_df.norm <- readRDS("int_output/batch_corrected_opt_Untreated_all_lightOn.Rds")
## baseline normalization ##
baseline_T <- (3600-30):(3600-1)
#parse the baseline period
dk_adpt <- training_df.norm[training_df.norm$start %in% baseline_T,]
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
rm(i)
rm(group.name)
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
#remove unnecessary variables
rm(list = c("dk_adpt", "beta", "dk_adpt.mean", "dk_adpt.prd", "dk_adpt.list"))
#subset the dataframe for start time between 5340 (-60 before light-offset) and 5540 (120 after light-offset)
#output.sub <- training_df.norm[training_df.norm$start %in% off_t,]
#assign genotype names as the dataframe name, i.e. parse the genotype dataframes
treatment_list <- list()
group_type <- unique(training_df.norm$genotype)
for(i in 1:length(group_type))
{
  treatment_list[[i]] <- training_df.norm[training_df.norm$genotype == group_type[i],]
  names(treatment_list)[i] <- as.character(group_type[i])
}
rm(i)
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
  genorm <- cbind(treatment_list[[i]][,1:5], genorm, treatment_list[[i]]$batch)
  #add the normalized data frame to the list
  treatment_list.blNorm[[i]] <- genorm
  #change the columns to the measurement names
  colnames(treatment_list.blNorm[[i]])[6:14] <- colnames(treatment_list[[i]])[6:14]
  rm(genorm)
}
#change the names of the list to the genotype
names(treatment_list.blNorm) <- names(treatment_list)
#save the baseline normalized data
Training_bl_df <- do.call("rbind", treatment_list.blNorm)
names(Training_bl_df)[15] <- "batch" 
saveRDS(Training_bl_df, "int_output/blNorm_opt_all_lightOn.Rds")
