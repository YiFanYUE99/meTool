#' Title pathway plot
#'
#' @param input pathway table
#' @import ggplot2
#' @import readxl
#' @import tidyverse
#' @return pathway fig
#' @export
#' @examples
#' #input<-"data/pathway_CK_H.xlsx"
#' #p<-pathway(input)
#' #ggsave(plot = p,'pic/pathway_CK_H.png',width = 7.5, height = 6, dpi = 300)
pathway<-function(input){
  metabolites<-read_excel(input)
  metabolites$`Match Status`<-as.numeric(metabolites$`Match Status`)
  metabolites$p<-as.numeric(metabolites$p)
  metabolites$`-log(p)`<-as.numeric(metabolites$`-log(p)`)
  metabolites$`Holm p`<-as.numeric(metabolites$`Holm p`)
  metabolites$FDR<-as.numeric(metabolites$FDR)
  metabolites$Impact<-as.numeric(metabolites$Impact)
  #作图
  size=0.04
  allsize<-6
  textsize<-4
  titlesize<-14
  xtisize<-12
  ytisize<-12
  xtesize<-8
  ytesize<-8
  legtesize<-6
  legtisize<-6
  xintercept<-0.1
  yintercept<--log10(0.05)
  p<-ggplot(metabolites,aes(Impact,`-log(p)`))+
    #添加水平误差线 虚线 dashed 颜色为灰色#999999
    geom_hline(yintercept = yintercept,linetype="dashed",color="#999999")+
    #添加垂直误差线
    geom_vline(xintercept = xintercept,linetype="dashed",color="#999999")+
    geom_point(aes(size= Impact,color=`Match Status`))+
    #scale_color_gradient2(low ="#33429a",mid = "#f9ed36",high = "#b11f24",#gradient 2种颜色,gradient2 3种颜色，gradientn自定义
    #                      midpoint = 3)+
    scale_color_gradientn(colors =c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25","#990000"))+
    theme_bw()+
    scale_size_continuous(range = c(0.2,8))+#更改点的大小
    guides(colour = guide_colorbar(title="Match Status",order = 1),#order固定legend的顺序
           size = guide_legend(title = "Pathway Impact",order = 2))+#size的标签
    scale_x_continuous(limits=c(min(metabolites$Impact),max(metabolites$Impact)+0.2))+#固定横坐标范围
    scale_y_continuous(limits=c(min(metabolites$`-log(p)`)-0.02, max(metabolites$`-log(p)`)+0.05))+#固定纵坐标范围
    #添加点的标签
    geom_text(aes(label= ifelse(Impact>xintercept&`-log(p)`> yintercept,`Pathway Name`,NA),
                  vjust=2,hjust=0.3,angle =25,size=size),color="black", family="serif")+
    xlab("Pathway Impact")+
    ylab("-log10(p)")+
    theme(panel.grid=element_blank(),
          text = element_text(family = "serif", size=16),
          plot.title = element_text(size = titlesize,hjust = 0.5,family = "serif"),#title 大小
          axis.title.x = element_text(face = "bold", angle = 0, size = xtisize),
          axis.title.y = element_text(face = "bold", angle = 90, size = ytisize),
          axis.text.x = element_text(face = "bold", angle = 25, hjust = 1,size = xtesize),#横坐标
          axis.text.y = element_text(face = "bold", angle = 25, hjust = 1,size = ytesize),
          legend.text = element_text(size = legtesize),#legend字体大小
          legend.title = element_text(size = legtisize), #设置legend标签字体大小
          plot.margin = unit(c(1, 1, 1, 1), "cm"),#下、左、上、右
          #标签位置
          legend.position = "right"#去除标签
    )
  return(p)

}
