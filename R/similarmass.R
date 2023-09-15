#' Title find similar mass
#'
#' @param input1 the filename of the first csv//  example: data/20230906_C18_5min_QC_MS1_NEG.csv
#' @param input2 the filename of the second csv//  example: data/20230906_C18_5min_QC_MS1_POS.csv
#' @param deltaRT the deviation of RT; If not the same chromatographic column, deltaRT=NULL
#' @param deltamass the deviation of mass; unit:ppm
#'
#' @return a dataframe showing similar mass
#' @export
#' @examples  setwd("D:/R_work/meTool")
#' # a<-similarmass("data/20230906_C18_5min_QC_MS1_NEG.csv","data/20230906_C18_5min_QC_MS1_POS.csv",deltaRT = 0.05,deltamass =500)
#' # b<-similarmass("data/20230906_C18_5min_QC_MS1_NEG.csv","data/20230906_HILIC_5min_QC_MS1_NEG.csv",deltaRT = NULL,deltamass =500)
#' # c<-similarmass("data/20230906_C18_5min_QC_MS1_POS.csv","data/20230906_HILIC_5min_QC_MS1_POS.csv",deltaRT = NULL,deltamass =500)
#' # d<-similarmass("data/20230906_HILIC_5min_QC_MS1_POS.csv","data/20230906_HILIC_5min_QC_MS1_NEG.csv",deltaRT = NULL,deltamass =500)
#' #write.csv(a,file = "output/C18_POS_NEG.csv")
#' #write.csv(b,file = "output/C18_HILIC_NEG.csv")
#' #write.csv(c,file = "output/C18_HILIC_NEG_POS.csv")
#' #write.csv(d,file = "output/HILIC_POS_NEG.csv")

similarmass<-function(input1, input2, deltaRT=0.05, deltamass=1000){
  MS1<-read.table(input1,header=T,,check.names = FALSE,sep = ",")
  MS2<-read.table(input2,header=T,,check.names = FALSE,sep = ",")
  if(is.null(deltaRT)){
    mylist<-NULL
    for (i in 1:dim(MS1)[1]) {
      for (j in 1:dim(MS2)[1]) {
        if(abs(MS1[i,1]-MS2[j,1])/MS2[j,1]*1000000<deltamass){
          mylist<-append(mylist,list(MS1[i,]))
          mylist<-append(mylist,list(MS2[j,]))
        }
      }
    }
  }else{
    mylist<-NULL
    for (i in 1:dim(MS1)[1]) {
      for (j in 1:dim(MS2)[1]) {
        if(abs(MS1[i,1]-MS2[j,1])/MS2[j,1]*1000000<deltamass & abs(MS1[i,2]-MS2[j,2])<deltaRT){
          mylist<-append(mylist,list(MS1[i,]))
          mylist<-append(mylist,list(MS2[j,]))
        }
      }
    }
  }
  tables<-data.table::rbindlist(mylist)
  return(tables)
}

