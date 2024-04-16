#' Title Get the VIP value of a pls model
#'
#' @param input the anaylzing filename
#' @param the compared groups
#' @import pls
#' @return a datframe of the vip value of all features
#' @export
#' @examples
#' #input<- "data/metabolitesg.csv"
#' #comgr=c("CK","H")
#' #vip<-calc_vip(input,comgr)
calc_vip<-function(input,comgr){
  metabolites <- read.table(input, header = T, row.names = 1,
                            check.names = FALSE, sep = ",")
  met<-metabolites[metabolites$Group==comgr[1]|metabolites$Group==comgr[2],1:(dim(metabolites)[2]-1)]
  X<-metabolites[metabolites$Group==comgr[1]|metabolites$Group==comgr[2],1:(dim(metabolites)[2]-1)]
  y<-metabolites[metabolites$Group==comgr[1]|metabolites$Group==comgr[2],dim(metabolites)[2]]
  y<-as.factor(y)
  # 建立 PLS-DA 模型
  X<-as.matrix(X)
  y<-as.numeric(y)
  pls_model <- plsr(y~X, ncomp = 2, scale = TRUE)
  # 提取模型的加载量和变量重要性
  # 提取主成分数量
  m <- pls_model$ncomp
  # 提取主成分得分矩阵
  T <- pls_model$scores
  # 提取载荷矩阵
  W <- pls_model[["loading.weights"]]
  p<-nrow(W)
  # 初始化存储 VIP 值的向量
  vip_values <- numeric(p)
  # 计算每个预测变量的 VIP 值
  for (j in 1:p) {
    numerator<-0
    denominator<-0
    part<-0
    for(k in 1:m){
      numerator <- numerator+T[,k]^2*W[j,k]^2
    }

    for (k in 1:m){
      denominator<-denominator+T[,k]^2
    }
    part<-numerator/denominator+part


    vip_values[j] <-sqrt(part*p)
  }
  vip_values<-as.data.frame(vip_values)
  rownames(vip_values)<-colnames(met)

  return(vip_values)
}

















