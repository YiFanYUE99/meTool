#' Title Get the VIP value of a pls model
#'
#' @param input the anaylzing filename
#' @param the compared groups
#' @import pls
#' @return a datframe of the vip value of all features
#' @export
#' @examples
#' #input<- "data/metabolitesg.csv"
#' #comgr=c("CK","H")
#' #tables<-metadiff(input,comgr)
#' #vipp_plot(tables,G=comgr)


vipp_plot<-function(table,G=c("CK","H"),color=c("#86FF86","#E0FFE0")){
  title<-paste0(G[1]," ",G[2]," ","VIP & p value")
  bg<-c("#E6C8FF","#F0E0FF")
  gm<-c("#FFE0E0","#FFE0E6")
  pointsize<-0.1
  textsize<-2
  linesize<-3.5
  titlesize<-15
  xtisize<-14
  ytisize<-14
  xtesize<-10
  ytesize<-10
  ggplot(table, aes(x = Metabolites, y = vip)) +
    geom_segment(aes(x= Metabolites, y= vip, color=vip,
                     xend=Metabolites, yend=1),
                 linewidth=linesize,
                 linetype=1)+#使用reorder()排序变量
    geom_point(aes(fill=log10(`p-value`), size=-log10(`p-value`)/pointsize), pch = 21, color=bg[1]) +
    geom_text(aes(label =sprintf("%.4f",`p-value`)),
              color = "black",
              size = textsize,
              family = "mono",
              fontface="bold")+
    scale_fill_continuous(high=color[2],low=color[1])+
    scale_color_continuous(high=color[1],low=color[2])+
    guides(fill ="none",color="none",size="none",linewidth="none")+
    labs(title = title, x = "Meatabolites", y = "vip")+
    xlim(rev(table$Metabolites))+#固定横坐标
    coord_flip()+
    #annotate("text",x=Inf, y=-Inf, label="(A)",hjust=1,vjust=-0.5)+
    theme(
      text = element_text(family = "mono", size=16),
      plot.title = element_text(size = titlesize,hjust = 0.5),#title 大小
      axis.title.x = element_text(size = xtisize),#x轴标题
      axis.title.y = element_text(size = ytisize),
      axis.text.x = element_text(face = "bold", angle = 0, hjust = 1,size = xtesize),#横坐标
      axis.text.y = element_text(face = "bold", angle = 25, hjust = 1,size = ytesize),
      legend.text = element_text(size = 8),#legend字体大小
      legend.title = element_text(size = 8), #设置legend标签字体大小
      panel.grid.major = element_line(color = gm[1]),  # 主要网格线颜色为红色
      panel.grid.minor = element_line(color = gm[2]),  # 次要网格线颜色为蓝色
      panel.background = element_rect(fill = bg[2]),  # 去除背景
      panel.border = element_rect(color = bg[1], fill = NA)#保留边框
    )

}
