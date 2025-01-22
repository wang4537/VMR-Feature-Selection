setwd("Y:/data/Feature_selection_framework/preprocessing/")
on_df <- readRDS("int_output/blNorm_opt_all_lightOn.Rds")
off_df <- readRDS("int_output/blNorm_opt_all_lightOff.Rds")

#### transoform the data ####
#subset the time frames
on.subset <- on_df[on_df$start %in% (3600-30):(3600+299),]
#rename the start from 0
on.subset$start <- on.subset$start-3600
on.subset$end <-on.subset$start+1 
#add the column for light-on label
on.subset$light <- rep("on", nrow(on.subset))
#replace the genotype names
on.subset$genotype <- gsub("*\\+DMSO", "", on.subset$genotype)
#light-on batch 1 to 3 different location convention
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

#subset the time frames
off.subset <- off_df[off_df$start %in% (7200-30):(7200+299),]
#rename the start from 3600
off.subset$start <- off.subset$start-7200
off.subset$end <- off.subset$start+1
#add the column for light-on label
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
on.subset$uid <- paste(on.subset$batch, on.subset$location, on.subset$genotype, sep = "_")
off.subset$uid <- paste(off.subset$batch, off.subset$location, off.subset$genotype, sep = "_")
#saved RDS
write.csv(on.subset, "../input_data/Opt_on.csv", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.csv(off.subset, "../input_data/Opt_off.csv", quote = FALSE, row.names = FALSE, col.names = TRUE)
#write the data to a .csv file for python input
#write.csv(on.subset, "py_input_data/Opt_on.csv", quote = FALSE, row.names = FALSE)
#write.csv(off.subset, "py_input_data/Opt_off.csv", quote = FALSE, row.names = FALSE)
