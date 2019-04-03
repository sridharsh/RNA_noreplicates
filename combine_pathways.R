library(reshape)
library(data.table)

filenames<- list.files("/GSEA_results/", pattern = "GO*.xls")

##Create list of data frame names without the ".csv" part 
names <- gsub(".xls","",filenames)
list_dataframes <- list()
adam_list_dataframes <- list()
names_list <- list()
up_down <- list()

###Load all files
for(i in 1:length(names)){
  filepath <- file.path("/GSEA_results/",paste(names[i],".xls",sep=""))
  list_dataframes[[i]] <- read.table(filepath, sep = "\t", skip = 1)
  list_dataframes[[i]]$pathway <- names[[i]]
  list_dataframes[[i]]$V1 <- list_dataframes[[i]]$V4 <- list_dataframes[[i]]$V5 <- list_dataframes[[i]]$V6 <- list_dataframes[[i]]$V3 <- list_dataframes[[i]]$V9 <- NULL
}
df <-rbindlist(list_dataframes)
write.table(df,"pathway_count.txt", quote = FALSE, row.names = FALSE, sep = "\t", col.names = FALSE)
