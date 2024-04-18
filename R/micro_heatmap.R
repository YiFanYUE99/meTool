#' Title draw differential microbiome
#'
#' @param input csv file//example: data/microbiome.csv
#' @param comgr compared group names(just 2)
#' @return the microbiome abundance heatmap
#' @export
#' @import pheatmap
#'
#' @examples
#' #input<- "data/microbiome.csv"
#' #comgr=c("CK","L")
#' #p<-metabo_heatmap(input,comgr)
#' #ggsave(plot = p,"pic/metabolites_CK_L.png",width = 20,height = 15,dpi = 300,units = "cm")


micro_heatmap<-function(input,comgr){
  table<-microdiff(input,comgr)
  mi<-read.csv(input,check.names = FALSE,row.names = 1)
  #简化菌群名称
  micro<-mi[,1:(dim(mi)[2]-1)]
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
  micro<-as.data.frame(scale(micro,center = T,scale = T))#中心化且归一化数据
  micro$Group<-mi$Group
  micro<-micro[micro$Group==comgr[1]|micro$Group==comgr[2],-dim(micro)[2]]

  hm<-matrix(ncol = dim(table)[1], nrow = dim(micro)[1])
  hm<-as.data.frame(hm)
  for (i in 1:dim(table)[1]) {
    for (j in 1:(dim(micro)[2]-1)) {
      if(colnames(micro)[j]==table[i,1]){
        hm[,i]=micro[,j]
      }
    }
  }
  colnames(hm)<-table$Microbiome
  rownames(hm)<-rownames(micro)
  p<-pheatmap(
    hm,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("#00E6E6","white","#FF7A00"))(100),
    main = paste0(comgr[1]," ",comgr[2]," ","Differential Microbiome Heatmap"),
    gaps_row = sum(mi$Group==comgr[1]),
    fontsize = 8,
    cellwidth = 15,
    cellheight = 15,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontfamily= "serif"
  )
  return(p)
}
