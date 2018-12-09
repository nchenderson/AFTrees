
#ConcIndex <- function(times, pred.times, status) {
#   indec <- which(status==1)
#  
#   IndexMat <- combn(indec, 2)
#   T1 <- times[IndexMat[1,]]
#   T2 <- times[IndexMat[2,]]
#    
#   That1 <- pred.times[IndexMat[1,]]
#   That2 <- pred.times[IndexMat[2,]]
#    
#   ans <- mean((T1 > T2)*(That1 > That2) + (T2 <= T1)*(That2 <= That1))
#   return(ans)
#}

Cstat <- function(mydata, risk.score, tau) 
{
    cens <- kmcens(mydata[, 1], mydata[, 2], tau)  ### look up this function
    GXi <- cens$surv[match(mydata[, 1], cens$distinct, nomatch = 1)]
    Wi <- (1/GXi/GXi) * mydata[, 2] * as.numeric(mydata[, 1] < tau)*mydata[,2]

    rs = risk.score
   
    cstat = concordance(mydata[, 1], mydata[, 2], rs, Wi)
    Z = list()
    Z$Dhat = cstat
    Z$rs = rs
    Z$cens.surv = cens$surv
    Z$cens.psii = cens$psii
    Z$distinct = cens$distinct
    Z$wt = Wi
    return(Z)
}

#> conc
#function (X, D, W, R) 
#{
#    out <- .Fortran("conc", n = as.integer(length(X)), time = as.integer(X * 
#        1000), status = as.integer(D), rs = as.integer(R * 1e+05), 
#        weight = as.double(W), hfwt = as.double(W/2), WK1 = as.double(0), 
#        CSTAT = as.double(0), PACKAGE = "survC1")
#    return(out$CSTAT)
#}


concordance <- function(time, status, rs, weight) {
     ### rs- risk score
     
      n <- length(time)
      USEP=0
      WK1=0
      ind <- status==1
      for(i in 1:n) {
          for(j in 1:n) {
             if(time[i] < time[j]) {
                   USEP=USEP+weight[i]
                   if(rs[i] > rs[j]) { 
                      WK1 = WK1 + weight[i] 
                   } 
                   if(rs[i] == rs[j]) { 
                      WK1 = WK1 + weight[i]/2 
                   }
             }
         }
      }
      CSTAT=WK1/USEP
      return(CSTAT) 
}

