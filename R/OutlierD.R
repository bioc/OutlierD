##########################################################################
#
#        
#    Outlier dectection using quantile regression on M-A scatterplots of high-throughput data
#
#             by
#
#      HyungJun Cho
#      Deparment of Statistics 
#      and Department of Biostatistics
#      Korea University
#
#      August 2007
#
##########################################################################

.First.lib <- function(lib, pkg) {  
     invisible()
     if(.Platform$OS.type=="windows" && require(Biobase) && interactive() 
        && .Platform$GUI=="Rgui") { addVigs2WinMenu("OutlierD") }
}

###################################################################################
#
# Main function  for users         
#
###################################################################################
OutlierD <- function(x1, x2, k=1.5, method="nonlin")
{

     ##########
     #Preparation
     x <- cbind(x1, x2) 
     n.col <- ncol(x) 
     n.row <- nrow(x) 
     if(n.row < 100) stop('Too few rows')
     if(is.data.frame(x)==FALSE) x <- data.frame(x)
     colnames(x) <- 1:ncol(x)           
     rownames(x) <- 1:nrow(x) 

     require(quantreg)
     cat("Please wait...")

    fit <- switch(method,
               constant = quant.const(x),
               linear      = quant.linear(x),
               nonlin     = quant.nonlin(x),
               nonpar   = quant.nonpar(x)
            )

     #Outputs
     x <- cbind(x, fit$A, fit$M, fit$Q1, fit$Q3, 
     fit$Q1-k*(fit$Q3-fit$Q1), fit$Q3+k*(fit$Q3-fit$Q1))
     colnames(x) <- c("X1","X2","A","M","Q1","Q3","LB","UB")

     i <-  which((x$M < x$LB)|(x$M > x$UB))
     n.outliers <- length(i)
     Outlier <- rep(FALSE, n.row)
     if(length(i) >0) Outlier[i] <- TRUE
     x <- cbind(Outlier, x)

     cat(" Done. \n")
     return(list(x=x, k=k, n.outliers=n.outliers, method=method)) 

}


####################################################################
quant.const <- function(x){

        M <- (x[,1]-x[,2])
        A <- (x[,1]+x[,2])/2
        Q2 <- quantile(M, probs=0.50)
        M <- M-Q2

        Q1 <- quantile(M, probs=0.25)
        Q3 <- quantile(M, probs=0.75)

        Q1 <- rep(Q1, length(A))
        Q3 <- rep(Q3, length(A))

        return(list(M=M, A=A, Q1=Q1, Q3=Q3))
}

####################################################################
quant.linear <- function(x){

        M <- (x[,1]-x[,2])
        A <- (x[,1]+x[,2])/2
        Q2 <- rq(M~A, tau=0.50)
        Q2 <- Q2$coef[1]+A*Q2$coef[2]
        M <- M-Q2

        Q3 <- rq(abs(M)~A, tau=0.50)
        Q3 <- Q3$coef[1]+A*Q3$coef[2]
        Q1 <- -Q3

        return(list(M=M, A=A, Q1=Q1, Q3=Q3))
}

####################################################################
quant.nonlin <- function(x){
 
        M <- (x[,1]-x[,2])
        A <- (x[,1]+x[,2])/2

        tmp <- try(nlrq(M ~ SSlogis(A, Asym, mid, scal), tau=0.50), silent = TRUE)
        if(class(tmp) !="try-error") {
           Q2 <- nlrq(M ~ SSlogis(A, Asym, mid, scal), tau=0.50)
           Q2 <- predict(Q2, newdata=list(x=A))
           M <- M-Q2
        }

  
        tmp <-  try(nlrq(abs(M) ~ SSlogis(A, Asym, mid, scal), tau=0.50), silent = TRUE)
        if(class(tmp)=="try-error") {
               fit <- quant.linear(x)
               return(list(M=fit$M, A=fit$A, Q1=fit$Q1, Q3=fit$Q3))
        }

        Q3 <- nlrq(abs(M) ~ SSlogis(A, Asym, mid, scal), tau=0.50)
        Q3 <- predict(Q3, newdata=list(x=A))
        Q1 <- -Q3

        if(length(which(abs(M) < abs(Q1))) < length(M)/4) {
               fit <- quant.linear(x)
               return(list(M=fit$M, A=fit$A, Q1=fit$Q1, Q3=fit$Q3))
        }

        return(list(M=M, A=A, Q1=Q1, Q3=Q3))
}


####################################################################
quant.nonpar <- function(x){

        lbda <- 1 

        #First Quantile
        M <- (x[,2]-x[,1])
        A <- (x[,2]+x[,1])/2
        Q1 <- rqss(M~qss(A, lambda=lbda), tau=0.75)
        i <- match(A, Q1$qss[[1]]$xyz[-1,1])
        y  <- Q1$coef[-1]
        Q1 <- Q1$coef[1] + y[i]
        Q1 <- -Q1

        #Third Quantile
        M <- (x[,1]-x[,2])
        A <- (x[,1]+x[,2])/2
        Q3 <- rqss(M~qss(A, lambda=lbda), tau=0.75) 
        i <- match(A, Q3$qss[[1]]$xyz[-1,1])
        y  <- Q3$coef[-1]
        Q3 <- Q3$coef[1] + y[i]

        return(list(M=M, A=A, Q1=Q1, Q3=Q3))
}

#END################################################################