#' Title draw all the metabolites' abundance plot
#'
#' @param input csv file//example: data/metabolitesg.csv
#' @param G the group need to be shown
#' @param color color of the group; the number of colors need to match the number of groups
#' @param filepath the output file path; all the pictures will be put in that file
#'
#' @return A batch of pictures
#' @export
#' @import ggplot2
#' @import agricolae
#'
#' @examples
#'    #met_pics(input="data/metabolitesg.csv",G=c("CK","L","H"),color=c("CK"="violet","L"="lightblue2","H"="lightgreen"),filepath="D:/R_work/meTool/pic/")
#'



met_pics<-function(input="data/metabolitesg.csv",G=c("CK","L","H"),
                   color=c("CK"="violet","L"="lightblue2","H"="lightgreen"),
                   filepath="D:/R_work/meTool/pic/"){
  data <- read.table(input, header = T, row.names = 1,
                     check.names = FALSE, sep = ",")
  for (i in 1:(dim(data)[2]-1)) {
    compound<-colnames(data)[i]
    output<-paste(filepath,"plot_compound_dose/",compound,".png",sep = "")
    df<-cbind(data[,dim(data)[2]],data[,i])
    df<-as.data.frame(df)
    colnames(df)<-c("Class","Value")
    df$Class<-as.factor(df$Class)
    df$Value<-as.numeric(df$Value)
    # 计算每个类别的均值和标准误差
    mean_values <- aggregate(Value ~ Class, data = df, FUN = mean)
    se_values <- aggregate(Value ~ Class, data = df, FUN = function(x) sd(x)/ sqrt(length(x)))

    #计算显著性
    # 执行单因子方差分析
    # 执行 Tukey's HSD 测试
    oneway1<-aov(Value~Class,data=df)
    out1<-LSD.test(oneway1,"Class",p.adj="bonferroni")
    mar<-out1$groups

    mean_values_index<- order(mean_values$Value, decreasing = TRUE)
    mv<-mean_values[mean_values_index, ]
    mv<-cbind(mv,mar[,2])
    colnames(mv)<-c("Class","Value","Label")

    p <- ggplot(mv, aes(x = Class, y = Value)) +
      geom_bar(stat = "identity", aes(fill= Class, color=Class), alpha =0.7 ) +
      scale_fill_manual(values = color)+
      scale_color_manual(values = color)+
      # 添加误差线
      geom_errorbar(data = mv, aes(ymin = Value - se_values$Value, ymax = Value + se_values$Value, color = Class),
                    width = 0.2,
                    position = position_dodge(0.9)) +
      # 添加显著性标记
      geom_text(data = mv, aes(label = Label,y=Value + se_values$Value, color = Class),vjust = -0.3, size = 6,family ="serif") +
      # 设置标签和标题
      labs(x = "", y = "Abundance",title = compound)+
      xlim(G)+
      theme(
        text = element_text(family = "serif",face="bold"),#更改字体为times new roman
        axis.text.x.bottom = element_text(size=16,angle = -45,
                                          hjust = 0,vjust= 0.5,face="bold"),#x轴刻度文字大小
        axis.text.y.left = element_text(size = 16,face="bold"), #y轴刻度文字大小
        axis.title.x =element_text(size=18,face= "bold"), # x轴title字体大小
        axis.title.y=element_text(size=18,face="bold"), # y轴title字体大小
        plot.title = element_text(size=20,hjust = 0.5,face="bold"),
        panel.border = element_rect(color = "black", fill = NA),#保留边框
        panel.background = element_blank(),
        legend.position = c(1.35,0.6),
        legend.key.height = unit(4, "cm"),#legend的高
        legend.key.width = unit(0.2, "cm"),#legend的宽
        plot.margin = unit(c(1,1,1,1),"cm")#调整边距t,r,b,l
      )
    ggsave(output, width=22,height=22,units = "cm",dpi=300,create.dir = T)
  }
}
