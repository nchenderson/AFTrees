RiskScore <- function(object, type=NULL) {
     rs <- exp(-object$m.train.mean)
     
     rs.test <- NULL
     if(!is.null(object$m.test)) {
        rs.test <- exp(-object$m.test.mean)
     }
     ans <- list(train=rs, test=rs.test)
     return(ans)
}
