\name{OutlierD}
\alias{OutlierD}
\title{Outlier dectection using quantile regression on the M-A scatterplots of high-throughput data}
\description{
This detects outliers using quantile regression on the M-A scatterplots of high-throughput data.
}
\usage{
     OutlierD(x1, x2, k=1.5, method="nonlin")
}

\arguments{

   \item{x1}{one n-by-1 vector for data (n= number of peptides, proteins, or genes} 
   \item{x2}{the other n-by-1 vector for data (n= number of peptides, proteins, or genes} 
   \item{k}{parameter in Q1-k*IQR and Q3+k*IQR, IQR=Q3-Q1, k=1.5 (default)}
   \item{method}{one of constant, linear, nonlinear, and nonparametric quantile regression}

}

\value{
    \item{x}{data and results for outliers}
}

\author{
   HyungJun Cho
}


\examples{

data(lcms)
x <- log2(lcms) #log2-tranformation, do normalization if necessary

fit1 <- OutlierD(x1=x[,1], x2=x[,2], method="constant")
fit2 <- OutlierD(x1=x[,1], x2=x[,2], method="linear")
fit3 <- OutlierD(x1=x[,1], x2=x[,2], method="nonlin")
fit4 <- OutlierD(x1=x[,1], x2=x[,2], method="nonpar")

fit3$x[1:10,]

plot(fit3$x$A, fit3$x$M, pch=".", xlab="A", ylab="M")
i <- sort.list(fit3$x$A)
lines(fit3$x$A[i], fit3$x$Q3[i], lty=2); lines(fit3$x$A[i], fit3$x$Q1[i], lty=2)
lines(fit3$x$A[i], fit3$x$LB[i]); lines(fit3$x$A[i], fit3$x$UB[i])
title("Nonlinear")

}

\keyword{models}


