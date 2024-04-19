#' Title draw 2d Relative Abundance plot
#' @param metabolites metabolites csv file
#' @param microbiome microbiome csv file
#' @param sub the metabolites or microbiome name
#' @return 2d relative Abundance plot
#' @export
#' @import tidyverse
#' @import RColorBrewer
#' @examples
#' #metabolites<-"data/metabolitesg.csv"
#' #microbiome<-"data/microbiome.csv"
#' #sub=c("glutamate","urea")
#' #p1<-metmicro2dplot(metabolites,microbiome,sub)
#' #ggsave("pic/glutamate_urea.png", plot = p1, width = 20, height = 20,dpi = 300,units = "cm")
metmicro2dplot<-function(metabolites,microbiome,sub=c("glutamate","urea")){
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
  if(any(colnames(met) == sub[1], na.rm = TRUE)){
    if(any(colnames(met) == sub[2], na.rm = TRUE)){
      color<-brewer.pal(nlevels(met$Group),"Accent")
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
        ggtitle("Metabolites")+
        theme(plot.title = element_text(size=16,family = "serif",face = "bold",hjust = 0.5),
              axis.title.y = element_text(size = 12, family = "serif",face="bold"),
              axis.title.x = element_text(size = 12, family = "serif",face="bold"),
              axis.text.x = element_text(size = 10,  family = "serif",face="bold"),
              axis.text.y = element_text(size = 10,  family = "serif",face="bold"),
              panel.background = element_rect(colour = "#FFF5CC",fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    }else{
      color<-brewer.pal(nlevels(met$Group),"Set2")
      indices1 <- which(colnames(met) == sub[1], arr.ind = TRUE)
      indices2 <- which(colnames(micro) == sub[2], arr.ind = TRUE)
      plot2d<-as.data.frame(cbind(met[,indices1],micro[,indices2]))
      colnames(plot2d)<-c(colnames(met)[indices1],colnames(micro)[indices2])
      plot2d$Group<-met$Group
      rownames(plot2d)<-rownames(met)
      #作图
      xlab<-sub[1]
      ylab<-sub[2]
      p3<-ggplot(data = plot2d,aes(x=plot2d[,1],y=plot2d[,2],color=Group))+
        stat_ellipse(segments =  nlevels(plot2d$Group)+1,aes(fill=Group),#添加置信区间;以Group为标准
                     type = "t",level = 0.95, na.rm = FALSE, geom ="polygon",alpha=0.25,color=NA, show.legend=F,inherit.aes = T)+
        geom_point()+
        labs(x=xlab,y=ylab,color="")+
        guides(fill="none")+
        scale_fill_manual(values = color)+ #为三个种类手动添加颜色
        scale_colour_manual(values = color)+
        ggtitle("Metabolites & Microbiome")+
        theme(plot.title = element_text(size=16,family = "serif",face = "bold",hjust = 0.5),
              axis.title.y = element_text(size = 12, family = "serif",face="bold"),
              axis.title.x = element_text(size = 12, family = "serif",face="bold"),
              axis.text.x = element_text(size = 10,  family = "serif",face="bold"),
              axis.text.y = element_text(size = 10,  family = "serif",face="bold"),
              panel.background = element_rect(colour = "#FFF5CC",fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    }

  }else{
    if(any(colnames(micro) == sub[2], na.rm = TRUE)){
      color<-brewer.pal(nlevels(met$Group),"Dark2")
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
        ggtitle("Microbiome")+
        theme(plot.title = element_text(size=16,family = "serif",face = "bold",hjust = 0.5),
              axis.title.y = element_text(size = 12, family = "serif",face="bold"),
              axis.title.x = element_text(size = 12, family = "serif",face="bold"),
              axis.text.x = element_text(size = 10,  family = "serif",face="bold"),
              axis.text.y = element_text(size = 10,  family = "serif",face="bold"),
              panel.background = element_rect(colour = "#FFF5CC",fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    }else{
      color<-brewer.pal(nlevels(met$Group),"Set3")
      indices1 <- which(colnames(micro) == sub[1], arr.ind = TRUE)
      indices2 <- which(colnames(met) == sub[2], arr.ind = TRUE)
      plot2d<-as.data.frame(cbind(micro[,indices1],met[,indices2]))
      colnames(plot2d)<-c(colnames(micro)[indices1],colnames(met)[indices2])
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
        ggtitle("Microbiome & Metabolites")+
        theme(plot.title = element_text(size=16,family = "serif",face = "bold",hjust = 0.5),
              axis.title.y = element_text(size = 12, family = "serif",face="bold"),
              axis.title.x = element_text(size = 12, family = "serif",face="bold"),
              axis.text.x = element_text(size = 10,  family = "serif",face="bold"),
              axis.text.y = element_text(size = 10,  family = "serif",face="bold"),
              panel.background = element_rect(colour = "#FFF5CC",fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    }
  }
  return(p3)

}
