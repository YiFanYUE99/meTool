#' Title Find differential bacteria
#'
#'
#' @param dataname  The file name of microbiome abundance
#' @import tidyverse
#' @import rstatix
#' @return Differential bacteria
#' @export
#'
#' @examples #a<-microdiff("data/microbiome.csv")
#'
microdiff<-function(dataname){
  mi<-read.csv(dataname,check.names = FALSE,row.names = 1)
  mic<-as.tibble(mi)%>%
    pivot_longer(cols = 1 : (dim(mi)[2]-1), names_to = "variable", values_to = "Abundance") %>%#改格式
    group_by(variable) %>%
    kruskal_test(Abundance~Group) %>%#批量计算kw检验p值
    filter(p<0.05)%>%
    select(variable,p)
  micro<-as.data.frame(mic)
  mylist<-NULL
  mylist1<-NULL
  for (i in 1:(dim(mi)[2]-1)) {
    for (j in 1:dim(micro)[1]) {
      if(colnames(mi)[i]==micro[j,1]){
        mylist<-append(mylist,list(mi[,i]))
        mylist1<-append(mylist1,list(colnames(mi)[i]))
      }
    }
  }
  names<-as.data.frame(mylist1)
  abundance<-as.data.frame(mylist)
  colnames(abundance)<-names
  rownames(abundance)<-row.names(mi)
  all<-cbind(abundance,Group=mi$Group)
  return(all)
}







