#' Title A function to make heatmap with mantel test
#'
#' @param metabolites A csv file contains sample names in the first column, Groups in the last column, and metabolites abundance in the middle.
#' @param micro Run micro<-microdiff("data/microbiome.csv") to find different microbiome.
#' @param num From which columns represent "Class", "Order", "Family", and "Genus"
#' @import devtools
#' @return A picture
#' @export
#' @examples
#' #micro<-microdiff("data/microbiome.csv")
#' #metabolites<-"data/metabolite.csv"
#' #num<-c(4,8,12,22)
#' #output<-"pic/micro_meta.png"
#' #p<-mantelheatmap(metabolites, micro, num)
#' #ggsave(output, p, width = 11, height = 6, dpi = 300)
#'


mantelheatmap<-function(metabolites,micro,num){
  metabo<-read.csv(metabolites, check.names = FALSE, header = 1,row.names = 1)
  metabo1<-metabo[,1:(dim(metabo)[2]-1)]
  Group<-metabo$Group
  micro1<-micro[,1:(dim(micro)[2]-1)]
  mantel <- mantel_test(micro1, metabo1,
                        mantel_fun = 'mantel',
                        spec_select = list(Phylum = 1:(num[1]-1),
                                           Class = num[1]:(num[2]-1),
                                           Order = num[2]:(num[3]-1),
                                           Family = num[3]:(num[4]-1),
                                           Genus = num[4]:dim(micro1)[2])) %>%
    mutate(r1 = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                    labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
           p1 = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                    labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
  cor2 <- correlate(metabo1)
  #绘制相关性热图
  p1<- qcorrplot(cor2,
                 grid_col = "grey50",
                 grid_size = 0.2,
                 type = "upper",
                 diag = FALSE) +
    theme(
      axis.text.x = element_text(color = "black",size = 5),
      axis.text.y = element_text(color = "black",size = 5, angle = 45),
      legend.position = c(-0.1,0.5),
      legend.key = element_blank()
    ) +
    geom_square() +
    scale_fill_gradientn(colours = c("#FF66CC", "#FFCCFF", "#f5f4f4", "#99CCFF","#3399FF"),
                         limits = c(-1, 1))
  p1
  #添加显著性标签：
  p2 <- p1+
    geom_mark(size = 1.7,
              only_mark = T,
              sig_level = c(0.05, 0.01, 0.001),
              sig_thres = 0.05,
              colour = 'black')
  p2
  #加上mantel
  p3<- p2+
    geom_couple(data = mantel,
                aes(colour = p1, size = r1),
                curvature = nice_curvature())+
    scale_size_manual(values = c(0.25, 0.5, 0.75)) + #连线粗细
    scale_colour_manual(values = c("#CFF5D2","#D8CCE6","#F2D2D2")) + #连线配色
    #修改图例：
    guides(size = guide_legend(title = "Mantel r",
                               override.aes = list(colour = "grey35"),
                               order = 2),
           colour = guide_legend(title = "Mantel p",
                                 override.aes = list(size = 3),
                                 order = 1),
           fill = guide_colorbar(title = "Pearson r", order = 3))+
    theme(panel.grid=element_blank(),
          #标签位置
          legend.position = c(0,0.7),
          #legend.justification = c(0,0.3)#相对位置
    )
  p3
  return(p3)
}








