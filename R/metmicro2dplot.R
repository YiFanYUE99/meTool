#' Title draw pls-da plot
#' @param metabolites metabolites csv file
#' @param microbiome microbiome csv file
#' @param sub the metabolites or microbiome name
#' @param color color of the group; the number of colors need to match the number of groups
#'
#' @return 2d relative Abundance plot
#' @export
#' @import cowplot
#' @import patchwork
#' @import tidyverse
#' @import aplot
#' @import ggplotify
#' @examples
#' #metabolites<-"data/metabolitesg.csv"
#' #microbiome<-"data/microbiome.csv"
#' #sub=c("glutamate","urea")
#' #color<-c("CK"="violet","L"="lightblue2","H"="lightgreen")
#' #p1<-metmicro2dplot(metabolites,microbiome,sub,color)
#' #ggsave("pic/glutamate_urea.png", plot = p1, width = 20, height = 20,dpi = 300,units = "cm")

metmicro2dplot<-function(metabolites,microbiome,sub=c("glutamate","urea"),color=c("CK"="violet","L"="lightblue2","H"="lightgreen")){
  meta<-read.table(metabolites,header=T,row.names=1,check.names = F,sep=",")
  met<-meta[,-dim(meta)[2]]
  met<-as.data.frame(scale(met,center = T,scale = T))#中心化且归一化数据
  met$Group<-as.factor(meta$Group)
  mic<-read.table(microbiome,header=T,row.names=1,check.names = F,sep=",")
  micro<-mic[,-dim(mic)[2]]
  miname<-as.data.frame(matrix(ncol=1,nrow=1))
  for (k in 1:dim(micro)[2]) {
    # 找到最后一个逗号的位置
    microname <- max(gregexpr(";", colnames(micro)[k])[[1]])
    # 提取最后一个逗号后的字符
    extracted <- substr(colnames(micro)[k],
                        start = microname + 1,
                        stop = nchar(colnames(micro)[k]))
    miname<-rbind(miname,extracted)
  }
  miname<-miname[-1,]
  colnames(micro)<-miname
  micro<-as.data.frame(scale(micro,center = T,scale = T))#中心化且归一化数据
  micro$Group=mic$Group
  #如果两个变量都是代谢物
  if(any(colnames(met) == sub[1], na.rm = TRUE)&any(colnames(met) == sub[2], na.rm = TRUE)){
    indices1 <- which(colnames(met) == sub[1], arr.ind = TRUE)
    indices2 <- which(colnames(met) == sub[2], arr.ind = TRUE)
    plot2d<-as.data.frame(cbind(met[,indices1],met[,indices2]))
    colnames(plot2d)<-c(colnames(met)[indices1],colnames(met)[indices2])
    plot2d$Group<-met$Group
    rownames(plot2d)<-rownames(met)
    #作图
    xlab<-sub[1]
    ylab<-sub[2]
    p3<-ggplot(data = plot2d,aes(x=plot2d[,1],y=plot2d[,2],color=Group))+
      stat_ellipse(segments =  nlevels(plot2d$Group)+1,aes(fill=Group),#添加置信区间;以Group为标准
                   type = "t",level = 0.95, na.rm = FALSE, geom ="polygon",alpha=0.25,color=NA, show.legend=F,inherit.aes = T)+
      geom_point()+labs(x=xlab,y=ylab,color="")+
      guides(fill="none")+
      scale_fill_manual(values = color)+ #为三个种类手动添加颜色
      scale_colour_manual(values = color)+
      theme(axis.title.y = element_text(size = 12, family = "serif",face="bold"),
            axis.title.x = element_text(size = 12, family = "serif",face="bold"),
            axis.text.x = element_text(size = 10,  family = "serif",face="bold"),
            axis.text.y = element_text(size = 10,  family = "serif",face="bold"),
            panel.background = element_rect(colour = "#FFF5CC",fill = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    p3
    p1 <- ggplot(data=plot2d, aes(x = plot2d[,1],fill = Group)) +
      geom_density(color="transparent",alpha = 0.25, show.legend = FALSE) +
      guides(fill="none")+
      labs(title = "Density Plot", x = "X", y = "Density")+
      scale_fill_manual(values = color)+ #为2个种类手动添加颜色
      scale_colour_manual(values = color)+
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      labs(x=NULL,y=NULL)+
      theme(plot.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),#去掉横坐标刻度
            axis.text.x = element_text(size = 8, family = "serif", face="italic"),
            panel.background = element_blank())
    p1
    p2 <- ggplot(data=plot2d, aes(x = plot2d[,2],fill =Group)) +
      geom_density(color="transparent",alpha = 0.25, show.legend = FALSE) +
      guides(fill="none")+
      labs(title = "Density Plot", x = "X", y = "Density")+
      scale_fill_manual(values = color)+ #为2个种类手动添加颜色
      scale_colour_manual(values = color)+
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      labs(x=NULL,y=NULL)+
      theme(plot.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),#去掉横坐标刻度
            axis.text.x = element_text(size = 8, family = "serif", face="italic"),
            panel.background = element_blank())+#去掉背景
      coord_flip()
    p2
    plot_layout(ncol = 2,nrow = 3)
    p <-p3%>%
      insert_top(p1,height=0.25)%>%
      insert_right(p2,width=0.25)%>%
      as.ggplot()
    p<-p+plot_annotation(title = "Metabolites")
  }else{
    indices1 <- which(colnames(micro) == sub[1], arr.ind = TRUE)
    indices2 <- which(colnames(micro) == sub[2], arr.ind = TRUE)
    plot2d<-as.data.frame(cbind(micro[,indices1],micro[,indices2]))
    colnames(plot2d)<-c(colnames(micro)[indices1],colnames(micro)[indices2])
    plot2d$Group<-micro$Group
    rownames(plot2d)<-rownames(met)
    #作图
    xlab<-sub[1]
    ylab<-sub[2]
    p3<-ggplot(data = plot2d,aes(x=plot2d[,1],y=plot2d[,2],color=Group))+
      stat_ellipse(segments =  nlevels(plot2d$Group)+1,aes(fill=Group),#添加置信区间;以Group为标准
                   type = "t",level = 0.95, na.rm = FALSE, geom ="polygon",alpha=0.25,color=NA, show.legend=F,inherit.aes = T)+
      geom_point()+labs(x=xlab,y=ylab,color="")+
      guides(fill="none")+
      scale_fill_manual(values = color)+ #为三个种类手动添加颜色
      scale_colour_manual(values = color)+
      theme(axis.title.y = element_text(size = 12, family = "serif",face="bold"),
            axis.title.x = element_text(size = 12, family = "serif",face="bold"),
            axis.text.x = element_text(size = 10,  family = "serif",face="bold"),
            axis.text.y = element_text(size = 10,  family = "serif",face="bold"),
            panel.background = element_rect(colour = "#FFF5CC",fill = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    p3
    p1 <- ggplot(data=plot2d, aes(x = plot2d[,1],fill = Group)) +
      geom_density(color="transparent",alpha = 0.25, show.legend = FALSE) +
      guides(fill="none")+
      labs(title = "Density Plot", x = "X", y = "Density")+
      scale_fill_manual(values = color)+ #为2个种类手动添加颜色
      scale_colour_manual(values = color)+
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      labs(x=NULL,y=NULL)+
      theme(plot.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),#去掉横坐标刻度
            axis.text.x = element_text(size = 8, family = "serif", face="italic"),
            panel.background = element_blank())
    p1
    p2 <- ggplot(data=plot2d, aes(x = plot2d[,2],fill =Group)) +
      geom_density(color="transparent",alpha = 0.25, show.legend = FALSE) +
      guides(fill="none")+
      labs(title = "Density Plot", x = "X", y = "Density")+
      scale_fill_manual(values = color)+ #为2个种类手动添加颜色
      scale_colour_manual(values = color)+
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      labs(x=NULL,y=NULL)+
      theme(plot.title = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),#去掉横坐标刻度
            axis.text.x = element_text(size = 8, family = "serif", face="italic"),
            panel.background = element_blank())+#去掉背景
      coord_flip()
    p2
    plot_layout(ncol = 2,nrow = 3)
    p <-p3%>%
      insert_top(p1,height=0.25)%>%
      insert_right(p2,width=0.25)%>%
      as.ggplot()
    p<-p+plot_annotation(title = "Microbiome")
  }
  return(p)

}
