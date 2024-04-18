#' Title draw volcano plot
#'
#' @param input the metabolites file//example:data/metabolitesg.csv
#' @param titlename Title name
#' @param G compared groups
#'
#' @return the volcano plot
#' @import ggplot2
#' @import tidyverse
#' @export
#'
#' @examples
#'    #setwd("D:/R_work/meTool")
#'    #G=c("CK","H")
#'    #a<-volcanoplot("data/metabolitesg.csv",G=G)
#'    #a
#'
#'

volcanoplot<-function(input, titlename="volcano plot",G=c("CK","L")){
  metabolites<-read.table(input,header=T,row.names=1,check.names = FALSE,sep = ",")
  metCK<-metabolites[metabolites$Group==G[1],1:dim(metabolites)[2]-1]
  metTR<-metabolites[metabolites$Group==G[2],1:dim(metabolites)[2]-1]
  met<-rbind(metCK,metTR)
  #创建火山图所需列表H/CK
  cnames<-c("FC","log2FC","p","log10p")
  rnames<-colnames(met)
  volcano_hck<-matrix(nrow = dim(metabolites)[2]-1 ,ncol=4,dimnames = list(rnames,cnames))
  #计算foldchange
  for (i in 1:dim(met)[2]){
    volcano_hck[i,1]<-mean(metTR[,i])/mean(metCK[,i])
  }
  volcano_hck[,2]<-log2(volcano_hck[,1])
  #计算p值 t.test
  q2<-c(rep(0,sum(metabolites$Group == G[1])),rep(1,sum(metabolites$Group == G[2])))#t检验第一个变量是数据；第二个是binary数据（就像类别
  for (i in 1:dim(met)[2]){
    q1<-met[,i]
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
    ggtitle(label = titlename)+
    theme(plot.title = element_text(size = 16,hjust = 0.5, family = "serif"),
          legend.title = element_text(size = 10,family = "serif"),
          legend.text = element_text(size = 10,family = "serif"),
          axis.title =  element_text(size = 12,family = "serif"),
          axis.text =  element_text(size = 10,family = "serif"),
          panel.grid=element_blank(),
          #标签位置
          legend.position = "right",
          legend.box = "horizontal"
    )+
    guides(col=guide_colorbar(title="-log10(p adjusted)"),
           size="none")+#去掉size的标签
    scale_x_continuous(limits=c(-max(abs(data$log2FC)),max(abs(data$log2FC))))+#固定横坐标范围
    scale_y_continuous(limits=c(0,max(abs(data$log10p))))+#固定纵坐标范围
    #添加点的标签
    geom_text(aes(label= ifelse(abs(log2FC)>0.5&log10p>-log10(0.05),label,NA),
                  size=1,
                  vjust=2,hjust=1,angle=0),
              color="#333333",family="serif")+
    xlab("log2Foldchange")+
    ylab("-log10(FDR)")

  return(p)
}
