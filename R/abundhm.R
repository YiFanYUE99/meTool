

#' Title
#'
#' @param tables A table that contains differential metabolites
#' @param title The name of the pictures
#'
#' @return A double heatmap
#' @export
#' @import ggDoubleHeat
#' @examples
#' #plot<-abundhm(tables,"CK & H")
#' #output<-"pic/doubleheatmap.png"
#' #ggsave(output, plot, width = 6, height = 6, dpi = 300)
#'
abundhm<-function(tables,title){
  p<-ggplot(data=tables,aes(x=Sample,y=metabolites))+
    geom_heat_grid(outside = Abundance,
                   outside_colors = c("#FFFFE0", "#FFFFB3", "#FFFF86", "#FFFF59", "#FFFF2C","#FFFF00","#E6E600"),
                   inside =  Pvalue,
                   inside_colors = c("#6E00E6","#8200FF","#9628FF", "#AA50FF", "#BE78FF", "#D2A0FF", "#E6C8FF"),
                   r=3*sqrt(2))+
    ggtitle(title)
  return(p)
}




