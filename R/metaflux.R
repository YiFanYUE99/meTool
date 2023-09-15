#' Title A function to cope with metabolic flux
#' @description Sometimes it may be slow if the dataset is big. Just be patient.
#' all the data used are save in /data/ file
#'
#' @param input the data you want to put in;see data/SXT-neg.csv
#' @param CKrepeat the repeat count of Control group
#' @param EXPrepeat the repeat count of every Control group
#'                  NOTICE: it can be c(3,3) or 3. However, make sure every experiment group has the same repeat count
#' @param EXPgroup  (Optional) the number of experiment groups( except CK ); it can count automatically
#' @param deltaRT the allowed deviation of retention time (default=0.01)
#' @param ddC the allowed deviation of the carbon atom (default =0.0002)
#' @param Cnum the number of carbon atoms that are marked by 13C (default =6)
#' @param deltaC the mass of a carbon atom(default=1.0033)
#'
#' @return a data frame that shows which peak may be the result of 13C isotope labelling
#' @export
#'
#' @examples \donttest{
#' setwd("D:/R_work/meTool")
#' neg25<-metaflux("data/20230728cell_flux_25pos.csv", CKrepeat = 3, EXPrepeat = 3, Cnum = 6)
#'
#' #neg<-metaflux("data/SXT-neg.csv", CKrepeat = 3, EXPrepeat = c(3,3), Cnum = 6)}
#'
metaflux<- function(input,CKrepeat, EXPrepeat, EXPgroup=length(EXPrepeat),deltaRT=0.01,ddC=0.0002,Cnum=6,deltaC=1.0033){
  input<-read.csv(file = input,header=T,check.names = F)
  mylist1<-NULL
  mylist2<-NULL
  a<-0
  EXParea<-vector(length = EXPgroup)
  for(i in 1:dim(input)[1]){
    RTdown=input[i,2]-deltaRT
    RTup=input[i,2]+deltaRT
    Massup=input[i,1]-Cnum*(deltaC-ddC)
    Massdown=input[i,1]-Cnum*(deltaC+ddC)
    #CK在前；experiment组在后
    #计算每行EXP面积
    #EXPrepeat可以为c(3,3)或3,但每个实验组重复个数需要相同
    for(l in 1:EXPgroup){
      EXParea[l]<-rowMeans(input[i,(3+CKrepeat+(l-1)*EXPrepeat[l]):(2+CKrepeat+l*EXPrepeat[l])])
    }
    #计算每行的CK面积
    CKarea<-rowMeans(input[i,3:(2+CKrepeat)])
    for(j in 1:EXPgroup){
      if(EXParea[j]>5000 & CKarea<5000 & (EXParea[j]/CKarea)>10){
        for (k in 1:dim(input)[1]) {
          if((RTdown<=input[k,2])&(input[k,2]<=RTup) & (Massdown<=input[k,1])&(input[k,1]<=Massup)){
            EXParea1=rowMeans(input[k,(3+CKrepeat+(j-1)*EXPrepeat[j]):(2+CKrepeat+j*EXPrepeat[j])])
            CKarea1<-rowMeans(input[k,3:(2+CKrepeat)])
            if(EXParea1>5000 & CKarea1>5000){
              a<-a+1
              mylist1<-append(mylist1,list(input[i,]))
              mylist1<-append(mylist1,list(input[k,]))
              mylist2<-append(mylist2,list(a))
              mylist2<-append(mylist2,list(j))
            }
          }
        }
      }
    }
  }
  tables<-data.table::rbindlist(mylist1)
  no<-t(as.data.frame(mylist2))
  mytable<-cbind(as.data.frame(no),tables)
  return(mytable)
}
