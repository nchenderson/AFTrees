cv.AFTree <- function(X, y, status, nfold=10, weights=NULL, ndpost=1000, nskip=100, keepevery=1,
                      ntree=200, sigdf=3, sigquant=.90, k=2.0) {

  nobs <- nrow(X)
  if(is.null(weights)) {
    weights <- rep(1, nobs)
  }
  ss <- sample(1:nobs)
  fold_num <- cut(ss, breaks=nfold, labels=F)

  CV_score <- rep(0, nfold)
  for(j in 1:nfold) {
    ind <- fold_num == j

    ytest <- y[ind]
    status_test <- status[ind]
    Xtest <- X[ind,]
    inv_wtest <- weights[ind]

    Xtrain <- X[!ind,]
    ytrain <- y[!ind]
    status_train <- status[!ind]

    cat("fold", j, "\n")
    aftFit <- AFTrees(Xtrain, ytrain,status=status_train,x.test=Xtest, ndpost=ndpost,nskip=nskip,keepevery=keepevery,
                      ntree=ntree, sigdf=sigdf, sigquant=sigquant, k=k)

    fitted.vals <- aftFit$yhat.test.mean
    CV_score[j] <- mean(status_test*inv_wtest*abs(fitted.vals - log(ytest)))
  }
  return(CV_score)
}
