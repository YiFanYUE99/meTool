% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/microdiff_all.R
\name{microdiff_all}
\alias{microdiff_all}
\title{Title Find differential bacteria; but not remove those insignificant microbiome}
\usage{
microdiff_all(dataname, comgr = c("CK", "L"))
}
\arguments{
\item{dataname}{The file name of microbiome abundance}

\item{comgr}{the compared groups}
}
\value{
bacteria with t-test p values of two groups
}
\description{
Title Find differential bacteria; but not remove those insignificant microbiome
}
\examples{
#b<-microdiff_all("data/microbiome.csv",comgr=c("CK","L"))

}
