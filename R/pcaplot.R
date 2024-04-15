#' Title draw 2D pca plot
#'
#' @param input csv file//example: data/metabolitesg.csv
#' @param color color of the group; the number of colors need to match the number of groups
#'
#' @return the PCA plot
#' @export
#' @import cowplot
#' @import patchwork
#' @import tidyverse
#' @import aplot
#' @import ggplotify
#'
#' @examples
#'      setwd("D:/R_work/meTool")
#'      a<-pcaplot("data/metabolitesg.csv",color= c("CK"="violet","L"="lightblue2","H"="lightgreen"))
pcaplot<-function(input, color = c("CK"="violet","L"="lightblue2","H"="lightgreen")){
  data<- read.table(input,header=T,row.names=1,check.names = F,sep=",")
  #scale. = TRUE表示分析前对数据进行归一化；#行为样本，列为特征
  com1 <- prcomp(data[,1:dim(data)[2]-1], center = TRUE,scale. = TRUE)
  summary(com1)
  #提取PC score；
  df1<-com1$x
  #将metabolites数据集的第24列类别合并进来；
  df1<-data.frame(df1,data[,dim(data)[2]])
  summ<-summary(com1)
  #设置横纵坐标
  xlab<-paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")#round保留小数点后两位
  ylab<-paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
  p3<-ggplot(data = df1,aes(x=PC1,y=PC2,color=data$Group))+
    stat_ellipse(aes(fill=data$Group),#添加置信区间;以Group为标准
                 type = "t",level = 0.95, na.rm = FALSE, geom ="polygon",alpha=0.25,color=NA, show.legend=F,inherit.aes = T)+
    geom_point()+labs(x=xlab,y=ylab,color="")+
    guides(fill="none")+
    scale_fill_manual(values = color)+ #为三个种类手动添加颜色
    scale_colour_manual(values = color)+
    theme(axis.title.y = element_text(size = 12, family = "serif",face="bold"),
          axis.title.x = element_text(size = 12, family = "serif",face="bold"),
          axis.text.x = element_text(size = 10,  family = "serif",face="bold"),
          axis.text.y = element_text(size = 10,  family = "serif",face="bold"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  p3
  #pc1的密度图
  p1<- ggplot(df1, aes(df1$PC1, fill = df1$data...dim.data..2..,color = df1$data...dim.data..2.. )) +
    geom_density(alpha = 0.15,show.legend = FALSE) +
    scale_fill_manual(values = color) + scale_colour_manual(values = color) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x=NULL,y=NULL)+#去掉横纵坐标名称
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),#去掉横坐标刻度
          axis.text.y = element_text(size = 8,family = "serif", face="italic"),
          panel.background = element_blank())#去掉背景
  p1

  #pc2的密度图
  p2<- ggplot(df1, aes(df1$PC2, fill = df1$data...dim.data..2..,color= df1$data...dim.data..2..)) +
    geom_density(alpha = 0.15, show.legend = FALSE) +
    scale_fill_manual(values = color) + scale_colour_manual(values = color) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x=NULL,y=NULL)+#去掉横纵坐标名称
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),#去掉横坐标刻度
          axis.text.x = element_text(size = 8, family = "serif", face="italic"),
          panel.background = element_blank())+#去掉背景
    coord_flip()
  p2
  #拼图
  plot_layout(ncol = 2,nrow = 3)
  p <-p3%>%
    insert_top(p1,height=0.25)%>%
    insert_right(p2,width=0.25)%>%
    as.ggplot()
  return(p)
}
