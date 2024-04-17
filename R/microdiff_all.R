#' Title Find differential bacteria; but not remove those insignificant microbiome
#'
#'
#' @param dataname  The file name of microbiome abundance
#' @param comgr the compared groups
#' @import tidyverse
#' @import rstatix
#' @return bacteria with t-test p values of two groups
#' @export
#'
#' @examples
#' #b<-microdiff_all("data/microbiome.csv",comgr=c("CK","L"))
#'
microdiff_all<-function(dataname,comgr=c("CK","L")){
  mi<-read.csv(dataname,check.names = FALSE,row.names = 1)
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
  micro$Group<-mi$Group
  micro<-micro[micro$Group==comgr[1]|micro$Group==comgr[2],-dim(micro)[2]]
  q2 <- c(rep(0, sum(mi$Group==comgr[1])), rep(1, sum(mi$Group==comgr[2])))
  sel<-matrix(NA, nrow = ncol(micro), ncol = 2)
  sel<-as.data.frame(sel)
  for (i in 1:dim(micro)[2]) {
    q1 <- micro[, i]
    sel[i,1]<-colnames(micro)[i]
    sel[i,2] <- t.test(q1 ~ q2)$p.value
  }
  #删除含nan的行
  sel<-sel[complete.cases(sel),]
  colnames(sel)<-c("Microbiome","p-value")
  sel<-sel[order(sel$`p-value`),]#升序排列
  return(sel)
}
