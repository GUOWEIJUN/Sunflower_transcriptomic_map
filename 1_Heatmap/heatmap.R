#install.packages('devtools')
#library(devtools)
#install_github("raivokolde/pheatmap")

library(RColorBrewer)
library(pheatmap)


# change folder
setwd('./')
# set color
cc = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))


#####
# read file
data=read.delim(header=T,file='FPKM_28sample.txt',row.names = 1,sep="\t")

#remove row that (SD=0)：
#dat1 <- log10(data+1)
# matrix <- data[apply(data, 1, function(x) sd(x)!=0),]


p=cor(data,data,method='pearson')

write.table(p,file = "cor_pearson.txt",quote = F, sep = "\t")

#log10(fpkm+1)


# plot heatmap                                   
h=pheatmap(p,main = "pearson heatmap",
           na_col = "grey",color=cc(100),
           cluster_rows =F,cluster_cols = F,
           show_colnames=T,show_rownames = F,
           #         cellwidth = 20, cellheight = 20,
           display_numbers=F,legend=T,angle_col = 45,treeheight_row = 15,treeheight_col = 15,border_color = "white") 

pdf("heatmap.pdf", height = 6, width = 6)
h
dev.off()

png("heatmap.png", height = 600, width = 600)
h
dev.off()

# memb <- cutree(h$tree_row, k = 10)
# write.table(memb,file='gene.group.xls',sep="\t",quote=F)

