feature_count <- function(MAT_NAME, F_NAME)
{
  Var_file_filter <- paste("../filter_wilc/", MAT_NAME, "/", MAT_NAME, "_", F_NAME , "_vol.csv", sep = "")
  Var_list_filter <- read.csv(Var_file_filter, row.names = 1)
  
  Var_file_embeded <- paste("../embedded_wilc/", MAT_NAME, "/", MAT_NAME, "_", F_NAME, "_sigVar_pass_list.csv", sep = "")
  Var_list_embeded <- read.csv(Var_file_embeded)
  
  intersect_var_names <- intersect(Var_list_filter$feature, Var_list_embeded[,1])
  union_var_names <- union(Var_list_filter$feature, Var_list_embeded[,1])
  
  return(list("filter" = Var_list_filter$feature,
              "embeded" = Var_list_embeded[,1],
              "intersect" = intersect_var_names,
              "union" = union_var_names))
}

#Main
setwd("Y:/data/Feature_selection_framework/survey/")
mat_names <- list.dirs("../transformed_data/", full.names = FALSE)[-1]
f_name <- c("Opt_on", "Opt_off")
feature_list <- lapply(mat_names, function(x) lapply(f_name, function(y) feature_count(x, y)))
names(feature_list) <- mat_names
