% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plsdaplot.R
\name{plsdaplot}
\alias{plsdaplot}
\title{Title draw pls-da plot}
\usage{
plsdaplot(input, G = c("CK", "L"), color = c(CK = "violet", L = "lightblue2"))
}
\arguments{
\item{input}{csv file//example: data/metabolitesg.csv}

\item{color}{color of the group; the number of colors need to match the number of groups}
}
\value{
the PLS-DA plot
}
\description{
Title draw pls-da plot
}
\examples{
     #先设置要展示的组的行数
     #G=c("CK","L")
     #再设置各个组的颜色
     #color= c("CK"="violet","L"="lightblue2")
     #b<-plsdaplot("data/metabolitesg.csv",G=G,color= color)
     #b
     #output<-"D:/R_work/master-s-thesis-project/pic/plsdaplot.png"
     #ggsave(output, plot = b, width = 20, height = 20,dpi = 300,units = "cm")
}
