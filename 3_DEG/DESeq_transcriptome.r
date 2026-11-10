#====设置安装源为清华大学TUNA镜像站点=====================
# options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
#=========================================================
#安装DESeq2
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

#加载DESeq2
rm(list=ls())
#setwd('./')
library(DESeq2)
library("ggplot2")
data <- read.delim('FPKM_28sample_rounding.txt', header=T,row.names=1,sep="\t",check.names = F)

#####
##### 
countMatrix <- data[,c(1,2,15,16)]    #1
countMatrix <- data[,c(3,4,17,18)]    #2
countMatrix <- data[,c(5,6,19,20)]    #3
countMatrix <- data[,c(7,8,21,22)]    #4
countMatrix <- data[,c(9,10,23,24)]   #5
countMatrix <- data[,c(11,12,25,26)]  #6
countMatrix <- data[,c(13,14,27,28)]  #7

c <- colnames(countMatrix)
a <- c(rep("WT",2), rep("MT",2))
group <- as.data.frame(cbind(c,a)) 
colnames(group) <- c("Sample","Group")
group <- group[c(3,4,1,2),]
#group=read.delim('group_NIP.txt',header=T,sep="\t")  #读取分组文件

dds <- DESeqDataSetFromMatrix(countMatrix, colData=group, ~ Group)
#dds <- DESeqDataSetFromMatrix(countMatrix, colData=coldata, ~ condition)

#过滤掉所有基因中表达量为0的基因
dds <- dds[ rowSums(counts(dds)) > 1, ]

dds <- DESeq(dds)

res <- results(dds)
out=cbind(res$log2FoldChange,res$pvalue,res$padj)
rownames(out)=row.names(res)
colnames(out)=c('log2(FoldChange)','Pvalue','FDR')
#summary(res) #总结
#===================FDR过滤=============
# FCcut=2
# FDRcut=0.05
# diff=out[(!is.na(res$padj) & res$padj<FDRcut) & abs(res$log2FoldChange)>abs(log2(FCcut)),]
#=======================================
##2.火山图
# 设置p_value和logFC的阈值
cut_off_pvalue = 0.05  #统计显著性
cut_off_logFC = 1     #差异倍数值

plot_df <- data.frame(res)
plot_df <- cbind(geneid=row.names(plot_df),
                 plot_df)
#plot_df <- plot_df[!is.na(plot_df$pvalue), ]


plot_df$change = ifelse(plot_df$pvalue < cut_off_pvalue & abs(plot_df$log2FoldChange) >= cut_off_logFC, 
                        ifelse(plot_df$log2FoldChange > cut_off_logFC ,'Up','Down'),
                        'No diff')
p <- ggplot(
  # 数据、映射、颜色
  plot_df, aes(x = log2FoldChange, y = -log10(pvalue), colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

pdf("valcano.pdf",width = 6,height = 4)
p
dev.off()
png("valcano.png",width = 600,height = 400)
p
dev.off()

#==============P值过滤==================
#与FDR过滤任选一个，首选FDR
FCcut=2
Pcut=0.05
diff=out[(!is.na(res$pvalue) & res$pvalue<Pcut) & abs(res$log2FoldChange)>abs(log2(FCcut)),]

#上下调数据统计
res_sub <- subset(res, abs(log2FoldChange) > 1 & pvalue <Pcut)
##测试是否准确
gene_counts <- counts(dds, normalized = TRUE)
gene_counts <- gene_counts[row.names(res_sub), ]
head(gene_counts)
##

gene_up <- row.names(res_sub[res_sub$log2FoldChange > 0, ])
up <- length(gene_up) # 1024
gene_down <- row.names(res_sub[res_sub$log2FoldChange < 0, ])
down <- length(gene_down) # 702
DEG_number <- data.frame(c(down,up),c("gene_down","gene_up")) 

colnames(DEG_number) <- c("Number","regulated")
write.table(DEG_number,file = "DEG.number.txt",quote = FALSE)

p1 <- ggplot(data=DEG_number, mapping=aes(x=regulated,y=Number,fill=regulated))+  geom_bar(stat="identity",width = 0.3)+
  scale_color_manual(values=c("#546de5", "#ff4757"))+
  geom_text(aes(label=Number), vjust=1.6, color="white", size=3.5)+  theme_bw()#theme_minimal()
pdf("DEG.number.barplot.pdf",width = 6,height = 4)
p1
dev.off()
png("DEG.number.barplot.png",width = 600,height = 400)
p1
dev.off()

#
#"#546de5", "#ff4757"
#"#56B4E9","#999999"
#ggsave("DEG.number.barplot.pdf")
#ggsave("DEG.number.barplot.png")
#=======================================

all.re <- merge(countMatrix,diff, by ="row.names")

write.table(out,"DESeq2_AllResult.txt", sep = "\t", quote=F,row.names = T,col.names=T)
write.table(all.re,"DESeq2_DEG.txt", sep = "\t", row.names = TRUE,quote=F)
write.table(rownames(diff),"DESeq2_DEGlist.txt", sep = "\t",quote=F,row.names=F,col.names=F)

