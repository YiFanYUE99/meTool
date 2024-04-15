#' Title draw pls-da plot
#'
#' @param input csv file//example: data/metabolitesg.csv
#' @param color color of the group; the number of colors need to match the number of groups
#'
#' @return the PLS-DA plot
#' @export
#' @import cowplot
#' @import caret
#' @import patchwork
#' @import tidyverse
#' @import aplot
#' @import ggplotify
#'
#' @examples
#'      #先设置要展示的组的行数
#'      #G=c("CK","L")
#'      #再设置各个组的颜色
#'      #color= c("CK"="violet","L"="lightblue2")
#'      #b<-plsdaplot("data/metabolitesg.csv",G=G,color= color)
#'      #b
#'      #output<-"D:/R_work/master-s-thesis-project/pic/plsdaplot.png"
#'      #ggsave(output, plot = b, width = 20, height = 20,dpi = 300,units = "cm")

plsdaplot<-function(input, G=c("CK","L"),color = c("CK"="violet","L"="lightblue2")){
  data<-read.table(input,header=T,row.names=1,check.names = F,sep=",")
  data1<-data[data$Group==G[1]|data$Group==G[2],2:dim(data)[2]-1]
  group<-as.factor(data[data$Group==G[1]|data$Group==G[2],dim(data)[2]])
  data1<-scale(data1,center = T,scale = T)#中心化且归一化数据
  a<-plsda(x = data1, group, ncomp =11)

  ascore<-matrix(nrow= dim(a$scores)[1],ncol =2)
  for (i in 1:dim(ascore)[1]){
    for (j in 1:dim(ascore)[2]){
      ascore[i,j]<-a$scores[i,j]
    }
  }

  colnames(ascore) <- c("Comp1","Comp2")
  ascore<-as.data.frame(ascore)
  ascore<-cbind(ascore,group)

  xlab<-paste0("Component1(",round(a$Xvar[1]/10,2),"%)")#round保留小数点后两位
  ylab<-paste0("Component2(",round(a$Xvar[2]/10,2),"%)")
  p3<-ggplot(data = ascore,aes(x=Comp1,y=Comp2,color = group))+
    stat_ellipse(aes(fill = group),#添加置信区间;以Group为标准
                 type = "t",level = 0.95, na.rm = FALSE, geom ="polygon",alpha=0.25,color=NA, show.legend=F,inherit.aes = T)+
    geom_point()+labs(x=xlab,y=ylab)+
    guides(fill=F)+
    scale_fill_manual(values = color)+ #为2个种类手动添加颜色
    scale_colour_manual(values = color)+
    theme(axis.title.y = element_text(size = 12, family = "serif",face="bold"),
          axis.title.x = element_text(size = 12, family = "serif",face="bold"),
          axis.text.x = element_text(size = 10, family = "serif", face="bold"),
          axis.text.y = element_text(size = 10, family = "serif", face="bold"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  p3
  #c1的密度图
  p1<- ggplot(ascore, aes(Comp1, fill = group, color=group)) +
    scale_fill_manual(values = color)+ #为2个种类手动添加颜色
    scale_colour_manual(values = color)+
    geom_density(alpha = 0.15,show.legend = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x=NULL,y=NULL)+#去掉横纵坐标名称
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),#去掉横坐标刻度
          axis.text.y = element_text(size = 8, family = "serif", face="italic"),
          panel.background = element_blank())#去掉背景
  p1
  #pc2的密度图
  p2<- ggplot(ascore, aes(Comp2, fill = group,color =group)) +
    scale_fill_manual(values = color)+ #为2个种类手动添加颜色
    scale_colour_manual(values = color)+
    geom_density(alpha = 0.15, show.legend = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x=NULL,y=NULL)+
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
  p
  return(p)
}
