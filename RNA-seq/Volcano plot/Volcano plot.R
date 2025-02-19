##set up a work directory##????????????????????????,??????????????????????????????????????????????????????,??????????????????????????????##
setwd("E:/L3vs0hTr5/L3no1vs0hTr5no3 flybaseID")


##volcano plot##


install.packages("IDPmisc")

library(IDPmisc)
##need to read the data from exported table, cannot use Deseq matrix or result matrix directly##
library(ggplot2)
library(ggrepel)
# f<-read.csv("0hrvsL3 for volcano.csv")
#f<-read.csv("InRDN_L3vsL3 Geneid for volca.csv")
f<-read.csv("5hS-L3-T5vsL3-T5 genename for volca.csv")
f<-na.omit(f) 

f$threshold = factor(ifelse(f$padj < 0.05 & abs(f$log2FoldChange) >= 2, 
                            ifelse(f$log2FoldChange>= 2 ,'Up','Down'),'N.S.'),
                     levels=c('Up','Down','N.S.'))



ggobj <- ggplot(f,aes(x=log2FoldChange,y=-log10(padj),color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#CC0000","#2f5688","#BBBBBB"))+#??????????????????
  geom_text_repel(
    data = f[f$padj<0.05&abs(f$log2FoldChange)>=2,],
    aes(label = Label),
    size = 6,max.overlaps = 10000,
    col="black",
    segment.color = "black", show.legend = FALSE )+#??????????????????????????????
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+#??????????????????
  theme(
    legend.title = element_blank()#?????????????????????
  )+
  ylab('-log10 (p-adj)')+#??????y?????????
  xlab('log2 (FoldChange)')+#??????x?????????
  geom_vline(xintercept=c(-2,2),lty=3,col="black",lwd=0.5) +#????????????|FoldChange|>=2
  geom_hline(yintercept =c(0,1.3),lty=3,col="black",lwd=0.5)#????????????padj<0.05 

# ggsave("0hrvsL3 for volcano.tiff", device = tiff, height = 8, width = 8,plot = ggobj)
#ggsave("InRDN_L3vsL3 Geneid for volca.tiff", device = tiff, height = 8, width = 8,plot = ggobj)
ggsave("5hS-L3-T5vsL3-T5 genename for volca.tiff", device = tiff, height = 8, width = 8,plot = ggobj)
