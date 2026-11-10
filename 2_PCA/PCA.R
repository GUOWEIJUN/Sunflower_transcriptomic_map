#PCA analysis
install.packages("FactoMineR")
install.packages("factoextra")
library("FactoMineR")
library("factoextra")

#change folders
setwd("./")

#read file of gene expression
count <- read.table("FPKM_28sample.txt",header = T,sep ="\t",row.names = 1)

count1 <- as.data.frame(t(count)) 
#1.PCA
count.pca <- PCA(count1, graph = FALSE)
PCA(count1)


group <- read.table("group.txt",header = T,sep ="\t",row.names = 1)


count1 <- cbind(count1,group)
#2.
pdf("PCA_group.pdf",width = 6, height = 4)
fviz_pca_ind(count.pca,
             # Only display points without showing the text. By default, both are displayed.
             geom.ind = "point",
             # Set classification categories
             col.ind = count1$group,
             # set color
             palette = c("#00AFBB", "#E7B800","#FC4E07"),
             #  Concentration ellipses
             #palette = "jco",
             addEllipses = FALSE,
             legend.title = "Groups",repel = TRUE # Avoid text overlapping
)
dev.off()

png("PCA_group.png",width = 600, height = 600)
fviz_pca_ind(count.pca,
             # Only display points without showing the text. By default, both are displayed.
             geom.ind = "point",
             # Set classification categories
             col.ind = count1$group,
             # set color
             palette = c("#00AFBB", "#E7B800","#FC4E07"),
             # Concentration ellipses
             #palette = "jco",
             addEllipses = FALSE,
             legend.title = "Groups",repel = TRUE # Avoid text overlapping
)
dev.off()

#####
#tree

tree=hclust(dist(t(count)),method = 'average')
pdf(file='sample clustering.pdf',w=10,h=7)
plot(tree,xlab="", sub="", cex = 0.7,main="Sample Clustering") 
dev.off()

tree=hclust(dist(t(count)),method = 'average')
png(file='sample clustering.png',w=1000,h=400)
plot(tree,xlab="", sub="", cex = 0.7,main="Sample Clustering") 
dev.off()
