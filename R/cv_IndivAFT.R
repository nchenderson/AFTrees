

cv.IndivAFT <- function(X, y, status, Trt, nfold=10, weights=NULL, ndpost=1000, nskip=100, keepevery=1,
                        ntree=200, sigdf=3, sigquant=.90, k=2.0) {

  nobs <- nrow(X)
  if(is.null(weights)) {
    weights <- rep(1, nobs)
  }

  X_aug <- cbind(X, Trt)  ## augmented design matrix which includes treatment
  X_base <- cbind(X, rep(0,nobs)) ## "test" matrix used to generate baseline predictions

  ind1 <- Trt==1
  ind0 <- Trt==0
  nobs1 <- sum(ind1)
  nobs0 <- sum(ind0)

  ss <- sample(1:nobs)
  fold_num <- cut(ss, breaks=nfold, labels=F)

  fold_num1 <- fold_num  ## fold labels for those with Trt==1
  fold_num0 <- fold_num  ## fold labels for those with Trt==0
  fold_num1[Trt==0] <- 0
  fold_num0[Trt==1] <- 0

  CV_score <- rep(0, nfold)
  for(j in 1:nfold) {
     ind1 <- fold_num1 == j
     ind0 <- fold_num0 == j
     indd1 <- which(fold_num1==j)
     indd0 <- which(fold_num0==j)

     ind <- ind1 | ind0  ### union of the two sets

     ytest <- y[indd1]
     Xtest <- X[indd1,]
     Xtest_base <- X_base[ind,]

     Xtrain <- X[!ind,]
     Xtrain_aug <- X_aug[!ind,]
     ytrain <- y[!ind]
     status_train <- status[!ind]

     cat("fold", j, "\n")
     aftFit <- AFTrees(Xtrain_aug, ytrain,status=status_train,x.test=X_base, ndpost=ndpost,nskip=nskip,keepevery=keepevery,
                    ntree=ntree, sigdf=sigdf, sigquant=sigquant, k=k)
     m0 <- aftFit$yhat.test.mean

     R <- abs(outer(m0[ind1],m0[ind0], FUN="-"))

     ### Only fit a test set at the covariate values in ind1
     tau <- length(indd1)
     Xtest0 <- cbind(X[indd1,], rep(1, tau))
     Xtest1 <- cbind(X[indd1,], rep(0, tau))
     Xtmp <- rbind(Xtest0, Xtest1)
     aftFit <- AFTrees(Xtrain_aug, ytrain,status=status_train,x.test=Xtmp, ndpost=ndpost,nskip=nskip,keepevery=keepevery,
                       ntree=ntree, sigdf=sigdf, sigquant=sigquant, k=k)
     obs_diff <- fitted_vals <- status_test <- inv_wtest <- rep(0, length(indd1))
     for(h in 1:tau) {
       tmp <- which.min(R[h,])
       if(length(tmp)==0) {
         print(R[h,])
       }
       a0 <- indd0[tmp]
       a1 <- indd1[h]

       status_test[h] <- status[a1]*status[a0]
       inv_wtest[h] <- weights[a1]*weights[a0]
       obs_diff[h] <- log(y[a1]) - log(y[a0])
       fitted_vals[h] <- aftFit$yhat.test.mean[h] - aftFit$yhat.test.mean[h + tau]
     }
     CV_score[j] <- mean(status_test*inv_wtest*abs(fitted_vals - obs_diff))
  }
  return(CV_score)
}

