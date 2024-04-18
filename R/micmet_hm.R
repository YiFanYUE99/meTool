#' Title The abundance heatmap of all microbiome and metabolites
#'
#' @param metabolites the metabolites table
#' @param microbiome the microbiome table
#' @param comgr the compared groups(2 and only 2)
#' @import ComplexHeatmap
#' @import circlize
#' @import dendextend
#' @import dendsort
#' @import gridBase
#' @import ggplot2
#' @return a picture of all the metabolites and microbiome relative abundance
#' @export
#' @examples
#' #micmet(microbiome="data/microbiome.csv", metabolites="data/metabolitesg.csv",comgr=c("CK","L"),output="pic/circle_hp1.png")
#'
micmet<-function(microbiome, metabolites,comgr,output){
  mycol<-colorRamp2(c(-1,0,1),c("#00FFFF","white","#FF0046"))

  mi<-read.csv(microbiome,check.names = FALSE,row.names = 1)

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
  micro<-as.data.frame(scale(micro,center = T,scale = T))
  micro1<-as.data.frame(t(micro[mi$Group==comgr[1],]))
  micro1<-micro1[complete.cases(micro1),]
  micro2<-as.data.frame(t(micro[mi$Group==comgr[2],]))
  micro2<-micro2[complete.cases(micro2),]


  me<-read.csv(metabolites,check.names = FALSE,row.names = 1)
  rownames(me)<-rownames(mi)
  met<-me[,1:dim(me)[2]-1]
  met<-as.data.frame(scale(met,center = T,scale = T))
  met1<-as.data.frame(t(met[me$Group==comgr[1],]))
  met1<-met1[complete.cases(met1),]
  met2<-as.data.frame(t(met[me$Group==comgr[2],]))
  met2<-met2[complete.cases(met2),]

  CK<-rbind(met1,micro1)
  TR<-rbind(met2,micro2)


  ann_row = data.frame(pathway=c(rep("Metabolites",dim(met1)[1]),
                                 rep("Microbiome",dim(micro1)[1])))#对行进行注释，用于后续的热图分裂
  ann_row2 = data.frame(pathway=c(rep("Metabolites",dim(met2)[1]),
                                  rep("Microbiome",dim(micro2)[1])))#对行进行注释，用于后续的热图分裂
  row.names(ann_row) = rownames(CK)
  row.names(ann_row2) = rownames(TR)
  ann_row <- as.matrix(ann_row)#在circlize函数中，需要为matrix
  ann_row2 <- as.matrix(ann_row2)



  png(output, width = 60, height = 60, res=300, unit="cm")
  circos.par(gap.after=c(20,14))#分开处空隙
  circos.heatmap(CK,col =mycol,
                 split=ann_row,
                 track.height=0.2,#圆环粗细
                 rownames.side = "outside",
                 rownames.col = ifelse(ann_row=="Metabolites","#FF8F00","#2C2CFF"),
                 rownames.cex=1,#字体大小
                 rownames.font = 1,#字体粗细
                 show.sector.labels = F,#显示种类
                 cluster=TRUE,#对行进行聚类
                 bg.border="pink",
                 dend.callback=function(dend,m,si) { #dend.callback：用于聚类树的回调，当需要对聚类树进行重新排序，或者添加颜色时使用包含的三个参数：dend：当前扇区的树状图；m：当前扇区对应的子矩阵；si：当前扇区的名称
                   color_branches(dend,k=10,col=1:10) #color_branches():修改聚类树颜色#聚类树颜色改为1，即单色/黑色
                 }

  )
  circos.heatmap(TR,
                 col = mycol,
                 split=ann_row2,
                 dend.side = "inside",#聚类树方向
                 track.height=0.2,#圆环粗细
                 cluster=TRUE,#对行进行聚类
                 dend.track.height = 0.18,#聚类树的高度
                 bg.border="pink",
                 dend.callback=function(dend,m,si) { #dend.callback：用于聚类树的回调，当需要对聚类树进行重新排序，或者添加颜色时使用包含的三个参数：dend：当前扇区的树状图；m：当前扇区对应的子矩阵；si：当前扇区的名称
                   color_branches(dend,k=10,col=1:10) #color_branches():修改聚类树颜色#聚类树颜色改为1，即单色/黑色
                 })
  circos.track(track.index=get.current.track.index(),panel.fun=function(x,y){
    if(CELL_META$sector.numeric.index==1){   #if(CELL_META$sector.numeric.index == 3) { # the last sector
      cn=colnames(CK)
      n=length(cn)
      circos.text(rep(CELL_META$cell.xlim[2],n)+convert_x(1,"mm"),#x坐标
                  (1:n)*1.2+7,#调整y坐标
                  col = "#009E73",
                  cn,cex=1,adj=c(0,1),facing="inside")}
  },bg.border=NA)

  circos.track(track.index=get.current.track.index(),panel.fun=function(x,y){
    if(CELL_META$sector.numeric.index==1){# the last sector
      cn=colnames(TR)
      n=length(cn)
      circos.text(rep(CELL_META$cell.xlim[2],n)+convert_x(1,"mm"),#x坐标
                  (1:n)*1.2+15,#调整y坐标,行距+距离中心距(1:n)*1.2+5,
                  col = "#D55E00",
                  cn,cex=1,adj=c(0,1),facing="inside")
    }

  },bg.border=NA)



  #添加标题
  lg=Legend(title="Abundance",col_fun = mycol,
            direction = c("vertical"))
  grid.draw(lg)
  #draw(lg,x=unit(0.9,"npc"),y=unit(0.6,"npc"),just=c("right","center"))

  dev.off()
  circos.clear()
}
