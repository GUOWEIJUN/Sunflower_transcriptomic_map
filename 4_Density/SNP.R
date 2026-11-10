#install.packages("CMplot")
library(CMplot)
mydata<-read.table("AZB_363_snp.maf5mis10.pos.txt",header=TRUE,'\t')
head(mydata)
# snp         chr       pos
# snp1_1    1        2041
# snp1_2    1        2062
# snp1_3    1        2190
CMplot(mydata,plot.type="d",bin.size=1e6,
       col=c("darkgreen","yellow", "red"),
       file="pdf",dpi=300)
CMplot(mydata,plot.type="d",bin.size=1e6,
       col=c("darkgreen","yellow", "red"),
       file="jpg",dpi=300)

#indel
mydata<-read.table("AZB_363_indel.maf5mis10.pos.txt",header=TRUE,'\t')
head(mydata)
# snp         chr       pos
# snp1_1    1        2041
# snp1_2    1        2062
# snp1_3    1        2190
CMplot(mydata,plot.type="d",bin.size=1e6,
       col=c("darkgreen","yellow", "red"),
       file="pdf",dpi=300)
CMplot(mydata,plot.type="d",bin.size=1e6,
       col=c("darkgreen","yellow", "red"),
       file="jpg",dpi=300)

