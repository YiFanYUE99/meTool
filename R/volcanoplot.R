#' Title draw volcano plot
#'
#' @param input the metabolites file//example:data/metabolitesg.csv
#' @param CK the Control group rows
#' @param TR the treatment group rows
#' @param legendp the position of legend
#' @param xlab the range of abscissa
#' @param ylab the range of ordinate
#' @param textsize the text size of the point
#'
#' @return the volcano plot
#' @import ggplot2
#' @import tidyverse
#' @export
#'
#' @examples
#'    setwd("D:/R_work/meTool")
#'    a<-volcanoplot("data/metabolitesg.csv",CK=1:6, TR=13:18,legendp=c(0.8,0.1), xlab=c(-2,2),ylab=c(0,5),textsize=0.8)
#'    a
#'
#'

volcanoplot<-function(input, CK=1:6, TR=13:18,legendp=c(0.73,0.8),
                      xlab=c(-2.5,2.5),ylab=c(0,4.6),textsize=0.8){
  metabolites<-read.table(input,header=T,row.names=1,check.names = FALSE,sep = ",")
  met<-metabolites[,1:dim(metabolites)[2]-1]
  #创建火山图所需列表H/CK
  cnames<-c("FC","log2FC","p","log10p")
  rnames<-colnames(met)
  volcano_hck<-matrix(nrow = dim(metabolites)[2]-1 ,ncol=4,dimnames = list(rnames,cnames))
  #计算foldchange
  for (i in 1:dim(met)[2]){
    volcano_hck[i,1]<-mean(met[TR,i])/mean(met[CK,i])
  }
  volcano_hck[,2]<-log2(volcano_hck[,1])
  #计算p值 t.test
  q2<-c(rep(0,length(CK)),rep(1,length(TR)))#t检验第一个变量是数据；第二个是binary数据（就像类别
  for (i in 1:dim(met)[2]){
    q1<-met[c(CK,TR),i]
    volcano_hck[i,3]<-t.test(q1~q2)$p.value
  }
  volcano_hck[,4]<--log10(volcano_hck[,3])

  volcano_hck<-as.data.frame(volcano_hck)
  volcano_hck<-na.omit(volcano_hck)
  #作图
  data<-volcano_hck
  data<-data%>% arrange(desc(FC))#按照FC大小降序排列
  la<-rownames(data)
  data$label<-la


  p<-ggplot(data,aes(log2FC,log10p))+
    #添加水平误差线 虚线 dashed 颜色为灰色#999999
    geom_hline(yintercept = -log10(0.05),linetype="dashed",color="#999999")+
    #添加垂直误差线
    geom_vline(xintercept = c(-0.5,0.5),linetype="dashed",color="#999999")+
    geom_point(aes(size=log10p,color=log10p))+#根据-logp改变颜色
    #scale_color_gradient2(low ="#33429a",mid = "#f9ed36",high = "#b11f24",#gradient 2种颜色,gradient2 3种颜色，gradientn自定义
    #                      midpoint = 3)+
    scale_color_gradientn(values = seq(0,1,0.2),#gradientn自定义
                          colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
    theme_bw()+
    scale_size_continuous(range = c(0.2,8))+#更改点的大小
    theme(panel.grid=element_blank(),
          #标签位置
          legend.position = legendp,
          legend.justification = c(0,0.3)#相对位置
    )+
    guides(col=guide_colorbar(title="-log10(p adjusted)"),
           size="none")+#去掉size的标签
    scale_x_continuous(limits=xlab)+#固定横坐标范围
    scale_y_continuous(limits=ylab)+#固定纵坐标范围
    #添加点的标签
    geom_text(aes(label= label, color=log10p,size=textsize,vjust=1.5,hjust=1))+
    xlab("log2Foldchange")+
    ylab("-log10(FDR)")

  return(p)
}
