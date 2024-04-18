#' Title Get the correlation picture of metabolites
#'
#' @param input the anaylzing filename
#' @import ggplot2
#' @import ggDoubleHeat
#' @return a picture of the metabolites correlation map
#' @export
#' @examples
#' #plot<-cor_met("data/metabolitesg.csv")
#' #plot
#' #ggsave(plot=plot,"pic/correlation_metabolites.png",width=20,height = 18,dpi=300,units="cm")
cor_met<-function(input="data/metabolitesg.csv"){
  data<-read.table(input,header=T,row.names=1,check.names = F,sep=",")
  data1<-data[,-dim(data)[2]]


  cormat<-as.data.frame(matrix(ncol=5,nrow=1))
  colnames(cormat)<-c("Metabolites1","Metabolites2","r-value","p-value","Sig")
  corcell<-as.data.frame(matrix(ncol=5,nrow=1))
  colnames(corcell)<-c("Metabolites1","Metabolites2","r-value","p-value","Sig")

  for (i in 1:dim(data1)[2]) {
    for (j in 1:dim(data1)[2]) {
      corcell[1,1]<-colnames(data1)[i]
      corcell[1,2]<-colnames(data1)[j]
      corcell[1,3]<-ifelse(i==j,NA,cor(data1[,i],data1[,j]))
      corcell[1,4]<-ifelse(i==j,NA,cor.test(data1[,i], data1[,j], method = "pearson")$p.value)
      corcell[1,5]<-ifelse(i==j,NA,getSig(corcell[1,4]))
      cormat<-rbind(cormat,corcell)
    }
  }
  cormat<-cormat[-1,]

  textsize<-1.5
  titlesize<-9
  xtisize<-7
  ytisize<-7
  xtesize<-6
  ytesize<-6
  legtesize<-5.5
  legtisize<-5.5
  p<-ggplot(data=cormat,aes(x=Metabolites1,y=Metabolites2))+
    geom_heat_grid(outside = `r-value`,
                   outside_colors = c("#86FF86","white","#FF869D"),
                   inside =  `p-value`,
                   inside_colors = c("#0000E6","white"),
                   r=3*sqrt(2))+
    ggtitle("Metabolites Correlation")+
    geom_text(aes(label = Sig),
              color = "#FFFFB3",
              size = textsize,
              family = "mono",
              fontface="bold")+
    xlim(colnames(data1))+
    ylim(colnames(data1))+
    theme(
      text = element_text(family = "mono", size=16),
      plot.title = element_text(size = titlesize,hjust = 0.5,family = "serif"),#title 大小
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(face = "bold", angle = 25, hjust = 1,size = xtesize),#横坐标
      axis.text.y = element_text(face = "bold", angle = 25, hjust = 1,size = ytesize),
      legend.text = element_text(size = legtesize),#legend字体大小
      legend.title = element_text(size = legtisize), #设置legend标签字体大小
      legend.position = "right"
    )
  return(p)
}



