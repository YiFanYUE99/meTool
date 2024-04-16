#' Title Find differential metabolites through t test, but export them all
#' @param input Input file name and its file place
#' @param comgr The two groups you want to compare
#'
#' @return A table contains differential metabolites abundance and their P values
#' @export
#' @import rstatix
#' @import tidyverse
#' @import pls
#' @examples
#' #input<- "data/metabolitesg.csv"
#' #comgr=c("CK","H")
#' #tables<-metadiff_all(input,comgr)
metadiff_all<-function(input,comgr=c("CK","H")){
  vip<-calc_vip(input,comgr)
  metabolites <- read.table(input, header = T, row.names = 1,
                            check.names = FALSE, sep = ",")
  met<-metabolites[metabolites$Group==comgr[1]|metabolites$Group==comgr[2],1:(dim(metabolites)[2]-1)]
  q2 <- c(rep(0, sum(metabolites$Group==comgr[1])), rep(1, sum(metabolites$Group==comgr[2])))
  sel<-matrix(NA, nrow = ncol(met), ncol = 2)
  sel<-as.data.frame(sel)
  for (i in 1:dim(met)[2]) {
    q1 <- met[, i]
    sel[i,1]<-colnames(met)[i]
    sel[i,2] <- t.test(q1 ~ q2)$p.value
  }
  sel<-cbind(sel,vip)
  colnames(sel)<-c("Metabolites","p-value","vip")
  sel<-sel[order(-sel$vip),]#降序排列
  return(sel)
}



