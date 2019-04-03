library(gtools)
library(ggplot2)
library(reshape2)
library(ggrepel)

### run combine counts python script to get results for the below part of the code ####
#load python rpkm data#

updown_data<- read.csv("rpkm.csv", row.names = "GeneID")
updown_data <- updown_data[which(updown_data$WT_counts != updown_data$KO_counts ),]
updown_data <- updown_data[which(updown_data$WT_counts > 0 & updown_data$KO_counts > 0),]
updown_data$foldchange <- foldchange2logratio(foldchange(updown_data$WT_counts,updown_data$KO_counts))
updown_data<- updown_data[order(-updown_data$foldchange),]
updown_data <- updown_data[which(updown_data$WT_RPKM >= 5 | updown_data$KO_RPKM >= 5),]
updown_data$GeneID <- row.names(updown_data)

p <- ggplot(data = updown_data, aes(x=log(WT_counts), y=log(KO_counts))) +
  geom_point(aes(colour =foldchange), alpha =0.3) + scale_colour_gradient(low = "blue", high = "red")+ 
  geom_text(aes(label=ifelse(abs(foldchange) > 4,as.character(GeneID),'')),
  hjust=0,vjust=0, size=2, fontface = "bold", 
  position=position_jitter(width=1,height=1)) + 
  scale_x_continuous(name ="WT", limits = c(0,16)) + 
  scale_y_continuous(name ="KO", limits = c(0,16))
pdf("scatter_DEgenes.pdf")
p 
dev.off()

x <- ggplot(data = updown_data, aes(x=log(WT_counts), y=log(ADAMTSL2_KO_counts))) +
  geom_point(aes(colour =foldchange), alpha =0.3) + scale_colour_gradient(low = "blue", high = "red")+ 
  geom_text(aes(label=ifelse(GeneID == "<insert gene of interest>",as.character(GeneID),'')), hjust=0, vjust=0, size=2, 
  fontface = "bold", ,position=position_jitter(width=1,height=1)) + 
  scale_x_continuous(name ="WT", limits = c(0,16)) + 
  scale_y_continuous(name ="KO", limits = c(0,16))
x

DEgenes <- updown_data
DEgenes$WT_RPKM <- DEgenes$KO_RPKM <- DEgenes$WT_counts <- DEgenes$KO_counts <- NULL

write.table(DEgenes,"DEgenes.rnk", quote = FALSE, row.names = FALSE, sep = "\t", col.names = FALSE)
write.table(DEgenes,"DEgenes.csv", quote = FALSE, row.names = FALSE, sep = ",", col.names = TRUE)

### run pathway parse python script to get results for the below part of the code ####

pathway_data <- read.table("output_updownpathways.tsv", header = TRUE, sep = "\t")
pathway_data <- pathway_data[which((pathway_data$upregulated + pathway_data$downregulated) > 50),]
pathway_data$downregulated <- -(pathway_data$downregulated)
pathway_data <- melt(pathway_data)
q<- ggplot(pathway_data, aes(fill=variable, y=value, x=Pathway)) + 
  geom_bar( stat="identity", alpha = 0.5, color="black") + coord_flip() + 
  theme(axis.text.y =element_text(face="bold", size = 5)) + 
  labs(title = "GSEA Results: Enriched Pathways", x = "Gene Ontology Pathways", y = "Number of enriched genes") +
  scale_fill_manual("Regulation",values=c("Red","Blue"))
pdf("pathway_updown.pdf")
q
dev.off()
