SurvivalProb <- function(object, time.points=NULL, test.only=FALSE, train.only=FALSE) {
     ### First check that all elements of time.point > 0
     if(sum(time.points <= 0)) {
         stop("All time points must be positive")
     }
     if(test.only & is.null(object$m.test)) {
         stop("Must have fitted for test set when using test.only=TRUE")
     }
     if(is.null(time.points)) {
         log.times <- log(object$times)
         time.points <- seq(min(log.times), max(log.times), length.out=11)
     }
     if(is.null(object$m.test)) {
          train.only = TRUE
     }
     ### Computation of mean survival curves.

     nsubjects <- ncol(object$m.train)
     ### fitted.values is nsamp x nsubjects
     nsamples <- nrow(object$locations)
     ngrid <- length(time.points)
     log.time.points <- log(time.points)

     if(test.only) {
         ntest <- ncol(object$m.test)
         if(length(time.points) > 1) {
             SS.test <- matrix(0.0, nrow=ntest, ncol=ngrid)

         ### TauMat should be nsamples x nclusters
             for(i in 1:ntest) {
                 for(k in 1:ngrid) {
                     Amat <- (log.time.points[k] - object$locations - object$m.test[,i])/object$sigma
                     SS.test[i,k] <- sum(pnorm(Amat, lower.tail=FALSE)*object$mix.prop)/nsamples
                 }
             }
             SS.test.mean <- colMeans(SS.test)
         }
         else {
             SS.test <- rep(0.0, ntest)
             for(i in 1:ntest) {
                  Amat <- (log.time.points - object$locations - object$m.test[,i])/object$sigma
                  SS.test[i] <- sum(pnorm(Amat, lower.tail=FALSE)*object$mix.prop)/nsamples
             }
         }
         SS.train <- SS.train.mean <- NULL
     } else if(train.only) {
         if(length(time.points) > 1) {
           SS.train <- matrix(0.0, nrow=nsubjects, ncol=ngrid)

           ### TauMat should be nsamples x nclusters
            for(i in 1:nsubjects) {
              for(k in 1:ngrid) {
                 Amat <- (log.time.points[k] - object$locations - object$m.train[,i])/object$sigma
                 SS.train[i,k] <- sum(pnorm(Amat, lower.tail=FALSE)*object$mix.prop)/nsamples
              }
            }
            SS.train.mean <- colMeans(SS.train)
         }
         else {
            SS.train <- rep(0.0, nsubjects)
            for(i in 1:nsubjects) {
               Amat <- (log.time.points - object$locations - object$m.train[,i])/object$sigma
               SS.train[i] <- sum(pnorm(Amat, lower.tail=FALSE)*object$mix.prop)/nsamples
            }
            SS.train.mean <- mean(SS.train)
         }
         SS.test <- SS.test.mean <- NULL
     } else if(!train.only & !test.only) {
        ntest <- ncol(object$m.test)

        if(length(time.points) > 1) {
             SS.test <- matrix(0.0, nrow=ntest, ncol=ngrid)
             SS.train <- matrix(0.0, nrow=ntrain, ncol=ngrid)

           ### TauMat should be nsamples x nclusters
            for(i in 1:ntest) {
              for(k in 1:ngrid) {
                  Amat <- (log.time.points[k] - object$locations - object$m.test[,i])/object$sigma
                  SS.test[i,k] <- sum(pnorm(Amat, lower.tail=FALSE)*object$mix.prop)/nsamples
              }
            }
            for(i in 1:nsubjects) {
              for(k in 1:ngrid) {
                  Amat <- (log.time.points[k] - object$locations - object$m.train[,i])/object$sigma
                  SS.train[i,k] <- sum(pnorm(Amat, lower.tail=FALSE)*object$mix.prop)/nsamples
              }
            }
            SS.test.mean <- colMeans(SS.test)
            SS.train.mean <- colMeans(SS.train)
         }
         else {
            SS.test <- rep(0.0, ntest)
            for(i in 1:ntest) {
                Amat <- (log.time.points - object$locations - object$m.test[,i])/object$sigma
                SS.test[i] <- sum(pnorm(Amat, lower.tail=FALSE)*object$mix.prop)/nsamples
            }
            SS.train <- rep(0.0, nsubjects)
            for(i in 1:nsubjects) {
                Amat <- (log.time.points - object$locations - object$m.train[,i])/object$sigma
                SS.train[i] <- sum(pnorm(Amat, lower.tail=FALSE)*object$mix.prop)/nsamples
            }
            SS.test.mean <- mean(SS.test)
            SS.train.mean <- mean(SS.train)
         }
      }
      ans <- list()
      class(ans) <- "aftsurvcurves"
      ans$Surv.train <- SS.train
      ans$Surv.test <- SS.test
      ans$Surv.train.mean <- SS.train.mean
      ans$Surv.test.mean <- SS.test.mean
      ans$time.points <- time.points
      return(ans)
}



