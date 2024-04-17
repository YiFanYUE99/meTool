#' Title Get the correlation map of Metabolites and microbiome
#'
#' @param metabolites the metabolites table
#' @param microbiome the microbiome table
#' @import ggplot2
#' @import ggDoubleHeat
#' @return a picture of the metabolites and microbiome correlation map
#' @export
#' @examples
#' #p<-cor_metmicro(metabolites = "data/metabolitesg.csv", microbiome = "data/microbiome.csv")
#' #ggsave(plot=p,"pic/cor_met_micro.png",width=60,height = 20,dpi=300,units="cm")
cor_metmicro<-function(metabolites,microbiome){
  met<-read.table(metabolites,header=T,row.names=1,check.names = F,sep=",")
  micro<-read.table(microbiome,header=T,row.names=1,check.names = F,sep=",")
  meta<-met[,1:(dim(met)[2]-1)]
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
  colnames(cormat)<-c("Metabolite","Microbiome","r-value","p-value","Sig")
  corcell<-as.data.frame(matrix(ncol=5,nrow=1))
  colnames(corcell)<-c("Metabolite","Microbiome","r-value","p-value","Sig")
  for (i in 1:dim(meta)[2]) {
    for(j in 1:dim(micro)[2]){
      corcell[1,1]<-colnames(meta)[i]
      corcell[1,2]<-colnames(micro)[j]
      corcell[1,3]<-cor(meta[,i],micro[,j])
      corcell[1,4]<-cor.test(meta[,i], micro[,j], method = "pearson")$p.value
      corcell[1,5]<-getSig(corcell[1,4])
      cormat<-rbind(cormat,corcell)
    }
  }
  cormat<-cormat[-1,]
  textsize<-1.5
  titlesize<-9
  xtisize<-7
  ytisize<-7
  xtesize<-7
  ytesize<-7
  legtesize<-5.5
  legtisize<-5.5
  p<-ggplot(data=cormat,aes(x=Microbiome,y=Metabolite))+
    geom_heat_grid(outside = `r-value`,
                   outside_colors = c("#00E6E6","white","#FF7A00"),
                   inside =  `p-value`,
                   inside_colors = c("#6E00E6","white"),
                   r=3*sqrt(2))+
    ggtitle("Metabolites & Microbiome Correlation")+
    geom_text(aes(label = Sig),
              color = "#FFFF00",
              size = textsize,
              family = "mono",
              fontface="bold")+
    xlim(colnames(micro))+
    ylim(colnames(meta))+
    theme(
      text = element_text(family = "mono", size=16),
      plot.title = element_text(size = titlesize,hjust = 0.5,family = "serif"),#title 大小
      axis.title.x = element_text(face = "bold", angle = 0, size = xtisize),
      axis.title.y = element_text(face = "bold", angle = 90, size = ytisize),
      axis.text.x = element_text(face = "bold", angle = 25, hjust = 1,size = xtesize),#横坐标
      axis.text.y = element_text(face = "bold", angle = 25, hjust = 1,size = ytesize),
      legend.text = element_text(size = legtesize),#legend字体大小
      legend.title = element_text(size = legtisize), #设置legend标签字体大小
      legend.position = "right"
    )
  return(p)

}
