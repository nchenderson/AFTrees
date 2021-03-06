IndivAFTSeparate = function(
   x.train, y.train, status, Trt, x.test=matrix(0.0,0,0),
   sigest=NA, sigdf=3, sigquant=.5,
   k=2.0,power=2.0, base=.95, nonparametric=TRUE,
   ntree=200,
   ndpost=1000, nskip=100,
   printevery=100, keepevery=1, keeptrainfits=TRUE,
   usequants=TRUE, numcut=100, printcutoffs=0,
   verbose=TRUE, scale="log")
{

   npatients <- length(y.train)
   x.train1 <- x.train[Trt==1,]
   x.train0 <- x.train[Trt==0,]
   x.test1 <- x.train[Trt==0,]
   x.test0 <- x.train[Trt==1,]
   
   y.train0 <- y.train[Trt==0]
   status0 <- status[Trt==0]
   y.train1 <- y.train[Trt==1]
   status1 <- status[Trt==1]
   tmp0 <- AFTrees(x.train0, y.train0, status0, x.test0, sigest, sigdf, sigquant, 
                  k, power, base, nonparametric,ntree,ndpost, nskip, printevery, 
                  keepevery, keeptrainfits, usequants, numcut, printcutoffs, verbose)
   tmp1 <- AFTrees(x.train1, y.train1, status1, x.test1, sigest, sigdf, sigquant, 
                  k, power, base, nonparametric,ntree,ndpost, nskip, printevery, 
                  keepevery, keeptrainfits, usequants, numcut, printcutoffs, verbose)
                  
   ControlEsts <- TrtEsts <- matrix(NA, nrow=ndpost, ncol=npatients)
   ControlEsts[,Trt==0] <- tmp0$yhat.train
   ControlEsts[,Trt==1] <- tmp0$yhat.test
   TrtEsts[,Trt==1] <- tmp1$yhat.train
   TrtEsts[,Trt==0] <- tmp1$yhat.test
   ThetaMat <- TrtEsts - ControlEsts
   
   ans <- list()
   ans$Theta <- ThetaMat

   return(ans)
}
