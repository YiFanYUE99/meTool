#' Title Get the correlation picture of microbiome
#'
#' @param input the anaylzing filename
#' @import ggplot2
#' @return a picture of the microbiome correlation map
#' @export
#' @examples
#' #plot<-cor_met("data/microbiome.csv")
#' #plot
#' #ggsave(plot=plot,"pic/correlation_microbiome.png",width=20,height = 18,dpi=300,units="cm")
cor_micro<-function(input="data/microbiome.csv"){
  micro<-read.table(input,header=T,row.names=1,check.names = F,sep=",")
  micro<-micro[,1:(dim(micro)[2]-1)]
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


  cormat<-as.data.frame(matrix(ncol=5,nrow=1))
  colnames(cormat)<-c("Microbiome1","Microbiome2","r-value","p-value","Sig")
  corcell<-as.data.frame(matrix(ncol=5,nrow=1))
  colnames(corcell)<-c("Microbiome1","Microbiome2","r-value","p-value","Sig")

  for (i in 1:dim(micro)[2]) {
    for (j in 1:dim(micro)[2]) {
      corcell[1,1]<-colnames(micro)[i]
      corcell[1,2]<-colnames(micro)[j]
      corcell[1,3]<-ifelse(i==j,NA,cor(micro[,i],micro[,j]))
      corcell[1,4]<-ifelse(i==j,NA,cor.test(micro[,i], micro[,j], method = "pearson")$p.value)
      corcell[1,5]<-ifelse(i==j,NA,getSig(corcell[1,4]))
      cormat<-rbind(cormat,corcell)
    }
  }
  cormat<-cormat[-1,]

  allsize<-12
  textsize<-2
  titlesize<-20
  xtisize<-14
  ytisize<-14
  xtesize<-8
  ytesize<-8
  legtesize<-8
  legtisize<-12
  color<-c("#00E6E6","white","#FF7A00")
  p<-ggplot(data=cormat,aes(x=Microbiome1,y=Microbiome2))+
    geom_tile(width=1,height =1,color = "lightgrey", fill="white")+
    geom_point(aes(size = abs(`r-value`), color = `r-value`))+
    coord_equal() + # 绘制正方形，长宽相等
    scale_color_gradientn(colors = color)+#legend上下限
    guides(size = "none",  #隐藏size图例
           color = guide_colorbar(title = "Pearson's r"))+#修改图例标签
    ggtitle("Microbiome Correlation")+
    geom_text(aes(label = Sig),
              color = "#FFFFB3",
              size = textsize,
              family = "mono",
              fontface="bold")+
    xlim(colnames(micro))+
    ylim(colnames(micro))+
    theme(
      text = element_text(family = "serif", size=allsize),
      plot.title = element_text(size = titlesize,hjust = 0.5,family = "serif"),#title 大小
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(face = "bold", angle = 45, hjust = 1,size = xtesize),#横坐标
      axis.text.y = element_text(face = "bold", angle = 25, hjust = 1,size = ytesize),
      legend.text = element_text(size = legtesize),#legend字体大小
      legend.title = element_text(size = legtisize), #设置legend标签字体大小
      legend.position = "right"
    )
  return(p)
}



