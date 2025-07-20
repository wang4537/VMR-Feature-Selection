Prep <- function(x) #x, filename of the data file
{
  ## preprocessing ##
  output <- read.csv(x, header = TRUE, sep = "\t") 
  out.columns <- colnames(output)[c(1,4,5,7,8,12,13,14,15,16,17,18,19,20)]#subset the output dataframe by the columns
  output.sub <- as.data.frame(output[out.columns])
  output.sub <- output.sub[output.sub$an == 2,] #center of the frame
  output.sub <- output.sub[output.sub$genotype %in% c("WT", "Q344X", "WT+DMSO", "Q344X+DMSO"), ]
  #keep only the desired time segment
  time_segment <- (3600-30):(7200+299)
  output.sub <- output.sub[output.sub$start %in% time_segment, ] 

  return (output.sub)
}

#### Main ####
library(readxl) #read excel
library(pbapply) #progress indicator
library(parallel) #parallel processing
setwd("/depot/yleung/data/Feature_selection_framework/2preprocess_sep/")

#### ON ####
#read files
data_files <- c(list.files("../../beichen temp/optimized DMSO Q344X/light-on/", full.names = TRUE),
                list.files("../../beichen temp/optimized untreated Q344X/light-on/", full.names = TRUE))

data_list <- mclapply(data_files, Prep, mc.cores = 4) #read files into a list
#assign batch number 
for(i in 1:length(data_list))
{
  data_list[[i]]$batch <- i
}
#bind together
on_df <- do.call('rbind', data_list)
#subset the time frames
on.subset <- on_df[on_df$start %in% (3600-30):(3600+299),]
#rename the start from 0
on.subset$start <- on.subset$start-3600
on.subset$end <-on.subset$start+1 
#add the column for light-on label
on.subset$light <- rep("on", nrow(on.subset))
#replace the genotype names
on.subset$genotype <- gsub("*\\+DMSO", "", on.subset$genotype)
#fix location issue in batch 1 to 3
on_4_loc <- on.subset[on.subset$genotype == "Q344X" & on.subset$batch == 4,"location"]
on_fix <- on.subset[on.subset$batch %in% 1:3,]
on_fix_list <- split(on_fix, as.factor(on_fix$batch))
on_fix_list <- lapply(on_fix_list, function(x)
{
  on_fix_W <- x[x$genotype == "WT",]
  on_fix_Q <- x[x$genotype == "Q344X",]
  on_fix_Q$location <-on_4_loc
  new_on <- rbind(on_fix_W, on_fix_Q)
  return(new_on)
}
)
#paste the new on_1 back to the main
on.subset <- rbind(do.call("rbind", on_fix_list),
                   on.subset[on.subset$batch %in% c(4:9),])
#add uid
on.subset$uid <- paste(on.subset$batch, on.subset$location, on.subset$genotype, sep = "_")
#save the master on data
write.csv(on.subset, "./Opt_on.csv", quote = FALSE, row.names = FALSE, col.names = TRUE)

#### OFF ####
#read files
#read the training data
data_files <- c(list.files("../../beichen temp/optimized DMSO Q344X/light-off", full.names = TRUE),
                list.files("../../beichen temp/optimized untreated Q344X/light-off", full.names = TRUE))[1:10]
data_list <- lapply(data_files, Prep) #read files into a list
#assign batch number 
for(i in 1:length(data_list))
{
  data_list[[i]]$batch <- i
}
#bind together
off_df <- do.call('rbind', data_list)
#subset the time frames
off.subset <- off_df[off_df$start %in% (7200-30):(7200+299),]
#rename the start from 0
off.subset$start <- off.subset$start-7200
off.subset$end <- off.subset$start+1 
#add the column for light-off label
off.subset$light <- rep("off", nrow(off.subset))
#replace the genotype names
off.subset$genotype <- gsub("*\\+DMSO", "", off.subset$genotype)
#light-off batch 1 to 5 different location convention
off_6_loc <- off.subset[off.subset$genotype == "Q344X" & off.subset$batch == 6,"location"]
off_fix <- off.subset[off.subset$batch %in% 1:5,]
off_fix_list <- split(off_fix, as.factor(off_fix$batch))
off_fix_list <- lapply(off_fix_list, function(x)
{
  off_fix_W <- x[x$genotype == "WT",]
  off_fix_Q <- x[x$genotype == "Q344X",]
  off_fix_Q$location <-off_6_loc
  new_off <- rbind(off_fix_W, off_fix_Q)
  return(new_off)
}
)
off.subset <- rbind(do.call("rbind", off_fix_list),
                    off.subset[off.subset$batch %in% c(6:9),])
#add uid
off.subset$uid <- paste(off.subset$batch, off.subset$location, off.subset$genotype, sep = "_")
#save the master on data
write.csv(off.subset, "./Opt_off.csv", quote = FALSE, row.names = FALSE, col.names = TRUE)
