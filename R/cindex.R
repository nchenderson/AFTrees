
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

Cstat <- function(times, status, risk.score, time.stop=NULL)
{
    if(is.null(time.stop)) {
        tau <- max(times)
    }
    else {
        tau <- time.stop
    }
    keep <- times < tau
    cens <- KMEstCens(times, status)  ### look up this function
    GXi <- cens$surv[match(times, cens$distinct, nomatch = 1)]
    Wi <- rep(0, length(GXi))
    Wi[keep] <- (1/(GXi[keep]*GXi[keep]))*status[keep]

    cstat = concordance(times, status, risk.score, Wi)
    ans <- list()
    ans$Chat <- cstat
    ans$cens.surv <- cens$surv
    return(ans)
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
KMEstCens <- function(time, status)  {
   distinct <- unique(sort(time))
   t <- length(distinct)
   n <- length(time)
   surv <- rep(0, t)
   yi <- sum(as.numeric(time >= distinct[1]))
   di <- sum(as.numeric(time == distinct[1] & status == 0))
   surv[1] <- 1 * (1 - di/yi)
   for (i in 2:t) {
      yi <- sum(as.numeric(time >= distinct[i]))
      di <- sum(as.numeric(time == distinct[i] & status == 0))
      surv[i] <- surv[i - 1] * (1 - di/yi)
   }
   return(list(surv = surv, distinct = distinct))
}


concordance <- function(time, status, rs, weight) {
     ### rs- risk score
     ### Need to write a C++ function to do this!

      n <- length(time)
      USEP=0
      WK1=0
      ind <- status==1
      for(i in 1:n) {
          for(j in 1:n) {
             if(time[i] < time[j] & i!= j) {
                  USEP <- USEP + weight[i]
                  if(rs[i] > rs[j]) {
                      WK1 <- WK1 + weight[i]
                  }
             }
             if((time[i] == time[j] || rs[i] == rs[j]) && i != j) {
                  USEP <- USEP + weight[i]/2
                  WK1 <- WK1 + weight[i]/4
             }
         }
      }
      CSTAT <- WK1/USEP
      return(CSTAT)
}

