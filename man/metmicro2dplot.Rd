% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metmicro2dplot.R
\name{metmicro2dplot}
\alias{metmicro2dplot}
\title{Title draw 2d Relative Abundance plot}
\usage{
metmicro2dplot(metabolites, microbiome, sub = c("glutamate", "urea"))
}
\arguments{
\item{metabolites}{metabolites csv file}

\item{microbiome}{microbiome csv file}

\item{sub}{the metabolites or microbiome name}
}
\value{
2d relative Abundance plot
}
\description{
Title draw 2d Relative Abundance plot
}
\examples{
#metabolites<-"data/metabolitesg.csv"
#microbiome<-"data/microbiome.csv"
#sub=c("glutamate","urea")
#p1<-metmicro2dplot(metabolites,microbiome,sub)
#ggsave("pic/glutamate_urea.png", plot = p1, width = 20, height = 20,dpi = 300,units = "cm")
}
