#' Title Find differential metabolites through t test
#'
#' @param input Input file name and its file place
#' @param comgr The two groups you want to compare
#'
#' @return A table contains differential metabolites abundance and their P values
#' @export
#' @import rstatix
#' @import tidyverse
#' @examples
#' #input<- "data/metabolite.csv"
#' #comgr=c("CK","H")
#' #tables<-metadiff(input,comgr)
metadiff<-function(input,comgr=c("CK","H")){
  metabolites<-read.csv(input,check.names = FALSE,header = 1)
  metabolite<-metabolites[,1:(dim(metabolites)[2]-1)]
  metabo<-cbind(Group=metabolites$Group,metabolite)
  meta<-as_tibble(metabo)%>%
    pivot_longer(cols=3:dim(metabo)[2],
                 names_to = "metabolites",
                 values_to = "Abundance")%>%
    filter(Group %in% comgr)%>%
    group_by(metabolites)

  metatest<-meta%>%
    t_test(Abundance~Group)%>%
    filter(p<0.05)

  mylist<-NULL
  mylist1<-NULL
  #筛选
  metad<-as.data.frame(meta)
  metatestd<-as.data.frame(metatest)
  for(i in 1:dim(metad)[1]){
    for (j in 1:dim(metatestd)[1]) {
      if(metatestd[j,1]==metad[i,3]){
        mylist<-append(mylist,list(metad[i,]))
        mylist1<-append(mylist1,list(metatestd[j,9]))
      }
    }
  }
  tables0<-data.table::rbindlist(mylist)
  tables1<-as.data.frame(mylist1)
  tables1<-t(tables1)
  colnames(tables1)<-"Pvalue"
  rownames(table1)<-NULL
  tables<-cbind(tables0,tables1)
  return(tables)
}



