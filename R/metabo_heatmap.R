#' Title draw differential metabolites
#'
#' @param input csv file//example: data/metabolitesg.csv
#' @param comgr compared group names(just 2)
#'
#' @return the metabolites abundance heatmap
#' @export
#' @import pheatmap
#'
#' @examples
#' #input<- "data/metabolitesg.csv"
#' #comgr=c("CK","L")
#' #p<-metabo_heatmap(input,comgr)
#' #ggsave(plot = p,"pic/metabolites_CK_L.png",width = 20,height = 15,dpi = 300,units = "cm")


metabo_heatmap<-function(input,comgr){
  table<-metadiff(input,comgr)
  data<-read.table(input,header=T,row.names=1,check.names = F,sep=",")
  data1<-data[data$Group==comgr[1]|data$Group==comgr[2],1:dim(data)[2]]
  data2<-data1[,1:(dim(data)[2]-1)]
  data2<-as.data.frame(scale(data2,center = T,scale = T))#中心化且归一化数据
  data2$Group<-data1$Group

  hm<-matrix(ncol = dim(table)[1], nrow = dim(data2)[1])
  hm<-as.data.frame(hm)
  for (i in 1:dim(table)[1]) {
    for (j in 1:(dim(data2)[2]-1)) {
      if(colnames(data2)[j]==table[i,1]){
        hm[,i]=data2[,j]
      }
    }
  }
  colnames(hm)<-table$Metabolites
  rownames(hm)<-rownames(data2)
  p<-pheatmap(
    hm,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("#2CFF2C", "white", "#FF2C63"))(100),
    main = paste0(comgr[1]," ",comgr[2]," ","Differential Metabolite Heatmap"),
    gaps_row = sum(data$Group==comgr[1]),
    fontsize = 8,
    cellwidth = 15,
    cellheight = 15,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontfamily= "serif"
  )
  return(p)
}














