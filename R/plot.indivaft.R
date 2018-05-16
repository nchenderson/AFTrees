plot.indivaft <- function(x, ...) {

     lower.quantile <- apply(x$Theta, 2, function(x) quantile(x, prob=.025))
     upper.quantile <- apply(x$Theta, 2, function(x) quantile(x, prob=.975))
     pmedian <- apply(x$Theta, 2, median)
     ind <- order(pmedian)
     lower.quantile <- lower.quantile[ind]
     upper.quantile <- upper.quantile[ind]
     pmedian <- pmedian[ind]
     n <- length(pmedian)
     yax <- seq(1, n, length.out=n)

     xlow <- min(lower.quantile)
     xhigh <- max(upper.quantile)
     plot(0,0, type="n", ylim=c(1,2142), xlim=c(xlow, xhigh), las=1, cex.axis=.8, xlab="treatment effect",
          ylab="Patient index")
     for(k in 1:n) {
       points(pmedian[k],yax[k], col="red", pch=16, cex=.6)
       lines(c(lower.quantile[k], upper.quantile[k]), c(yax[k], yax[k]),lwd=.8)
     }
     abline(v=0, lwd=3, col="grey")
}
