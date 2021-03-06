IndivAFTNew = function(
   x.train, y.train, status, Trt, x.test=NULL,
   sigest=NA, sigdf=3, sigquant=.5,
   k=2.0,power=2.0, base=.95, nonparametric=TRUE,
   ntree=200,
   ndpost=1000, nskip=100,
   printevery=100, keepevery=1, keeptrainfits=TRUE,
   usequants=TRUE, numcut=100, printcutoffs=0,
   verbose=TRUE, scale="log")
{

   if(is.vector(x.train) | is.factor(x.train)) x.train = data.frame(x=x.train)

   if(!is.matrix(x.train)) {
        stop("x.train must be a matrix")
   }
   if(is.null(x.test)) {
      x.train1 <- matrix(x.train[Trt==1,], nrow=sum(Trt==1), ncol=ncol(x.train))
      x.train0 <- matrix(x.train[Trt==0,], nrow=sum(Trt==0), ncol=ncol(x.train))
      print("hello")
      x.test1 <- x.train0
      x.test0 <- x.train1
   } else {
      x.train1 <- x.train[Trt==1,]
      x.train0 <- x.train[Trt==0,]
      x.test1 <- x.test[Trt==1,]
      x.test0 <- x.test[Trt==0,]
   }
   if(ncol(x.train)==1) {
      x.train1 <- matrix(x.train1, nrow=sum(Trt==1), ncol=1)
      x.train0 <- matrix(x.train0, nrow=sum(Trt==0), ncol=1)
      x.test1 <- matrix(x.test1, nrow=sum(Trt==0), ncol=1)
      x.test0 <- matrix(x.test0, nrow=sum(Trt==1), ncol=1)
   }
   print(dim(x.train0))
   #check input arguments:

   ncskip = floor(nskip/keepevery)
   ncpost = floor(ndpost/keepevery)
   nctot = ncskip + ncpost
   totnd = keepevery*nctot

   ################################################################################
   ### Fit things in Control Arm
   y.train0 <- y.train[Trt==0]
   status0 <- status[Trt==0]
   null_lm0 <- survreg(Surv(y.train0, status0) ~ 1, dist="lognormal")
   null_sig0 <- null_lm0$scale
   null_intercept0 <- null_lm0$coefficients


   y_centered_log0 <- log(y.train0) - null_intercept0
   imr0 <- dnorm(log(y.train0), mean=null_intercept0, sd=null_sig0)/pnorm(log(y.train0), mean=null_intercept0, sd=null_sig0, lower.tail=FALSE)
   #y_unobs <- y.train*exp(.02*(1-status)*null_sig)   ### rough, estimate of unobserved survival times
   y_unobs0 <- exp( status0*log(y.train0) + (1 - status0)*null_sig0 + (1 - status0)*imr0)

   #y.train.log = log(y.train)  ### transform to log-survival times
   rgy = range(log(y_unobs0))
   zeta_tmp0 <- rgy[2] - rgy[1]
   zeta0 <- 4*null_sig0
   #aa = (rgy[1] + rgy[2])/(2*(rgy[2] - rgy[1]))




   # Need to include the survival package
   if (is.na(sigest)) {
       tmplm <- survreg(Surv(y.train0, status0) ~ x.train0, dist="lognormal")
       #tmplm <- survreg(Surv(y.train, status) ~ x.train, dist="lognormal")
       sigest0 <- tmplm$scale
   } else {
       sigest0 <- sigest #put input sigma estimate on transformed scale
   }

   nclust = 50
   kappa0 <- FindKappa(q=sigquant, sigsq.hat=sqrt(sigest0), nu=sigdf)

   print(nrow(x.train0))
   print(nctot)
   npind = ifelse(nonparametric, 1, 0)
   cres0 = .C('mbart',as.integer(nrow(x.train0)), as.integer(ncol(x.train0)), as.integer(nrow(x.test0)),
                   as.double(x.train0), as.double(y_centered_log0),
                   as.double(x.test0), as.integer(status0),
                   as.double(sigest0),   as.integer(sigdf), as.double(sigquant),
                   as.double(k),
		   as.double(power), as.double(base),
		   as.integer(ntree),      as.integer(totnd),
                   as.integer(printevery), as.integer(keepevery),  as.integer(keeptrainfits),
                   as.integer(numcut), as.integer(usequants), as.integer(printcutoffs),
		   as.integer(verbose),
                   sdraw=double(nctot),
                   trdraw=double(nrow(x.train0)*nctot),
                   tedraw=double(nrow(x.test0)*nctot),
                   vcdraw=integer(ncol(x.train0)*nctot), as.integer(nclust),
                   mixdraw=double(nclust*nctot),
                   locdraw=double(nclust*nctot),
                   Mdraw=double(nctot),npind=as.integer(npind),as.double(kappa0),
                   as.double(zeta0))
    print("clear1")
    ################################################################################
    ### Fit things in Control Arm
   y.train1 <- y.train[Trt==1]
   status1 <- status[Trt==1]
   null_lm1 <- survreg(Surv(y.train1, status1) ~ 1, dist="lognormal")
   null_sig1 <- null_lm1$scale
   null_intercept1 <- null_lm1$coefficients


   y_centered_log1 <- log(y.train1) - null_intercept1
   imr1 <- dnorm(log(y.train1), mean=null_intercept1, sd=null_sig1)/pnorm(log(y.train1), mean=null_intercept1, sd=null_sig1, lower.tail=FALSE)
   #y_unobs <- y.train*exp(.02*(1-status)*null_sig)   ### rough, estimate of unobserved survival times
   y_unobs1 <- exp( status1*log(y.train1) + (1 - status1)*null_sig1 + (1 - status1)*imr1)

   #y.train.log = log(y.train)  ### transform to log-survival times
   rgy = range(log(y_unobs1))
   zeta_tmp1 <- rgy[2] - rgy[1]
   zeta1 <- 4*null_sig1
   #aa = (rgy[1] + rgy[2])/(2*(rgy[2] - rgy[1]))


   # Need to include the survival package
    if (is.na(sigest)) {
        #tmplm <- survreg(Surv(y.train, status) ~ x.train, dist="lognormal")
        tmplm <- survreg(Surv(y.train1, status1) ~ x.train1, dist="lognormal")
        sigest1 <- tmplm$scale
    } else {
        sigest1 <- sigest #put input sigma estimate on transformed scale
    }

    nclust = 50
    kappa1 <- FindKappa(q=sigquant, sigsq.hat=sqrt(sigest1), nu=sigdf)

    print("hello2")
    cres1 = .C('mbart',as.integer(nrow(x.train1)), as.integer(ncol(x.train1)), as.integer(nrow(x.test1)),
                   as.double(x.train1), as.double(y_centered_log1),
                   as.double(x.test1), as.integer(status1),
                   as.double(sigest1),   as.integer(sigdf), as.double(sigquant),
                   as.double(k),
		   as.double(power), as.double(base),
		   as.integer(ntree),      as.integer(totnd),
                   as.integer(printevery), as.integer(keepevery),  as.integer(keeptrainfits),
                   as.integer(numcut), as.integer(usequants), as.integer(printcutoffs),
		   as.integer(verbose),
                   sdraw=double(nctot),
                   trdraw=double(nrow(x.train1)*nctot),
                   tedraw=double(nrow(x.test1)*nctot),
                   vcdraw=integer(ncol(x.train1)*nctot), as.integer(nclust),
                   mixdraw=double(nclust*nctot),
                   locdraw=double(nclust*nctot),
                   Mdraw=double(nctot),npind=as.integer(npind),as.double(kappa1),
                   as.double(zeta1))

     print("hello3")
     ##################################################################################
   # now read in the results...
       ## look at this: sigma is multiplied by (rgy[2] - rgy[1])
       sigma = cres0$sdraw
      # sigma = cres$sdraw
       first.sigma = sigma[1:ncskip] # we often want the sigma draws
       sigma = sigma[ncskip+(1:ncpost)]

       # put sigest on the original y scale for output purposes
       sigest = sigest0

   yhat.train0 = yhat.test0 = yhat.train.mean0 = yhat.test.mean0 = NULL
   yhat.train1 = yhat.test1 = yhat.train.mean1 = yhat.test.mean1 = NULL
   varcount = NULL

   if (keeptrainfits) {
      yhat.train0 <- matrix(cres0$trdraw, nrow=nctot, byrow=TRUE)[(ncskip+1):nctot,]
      yhat.train0 <- yhat.train0 + null_intercept0
      yhat.test0 = matrix(cres0$tedraw,nrow=nctot,byrow=TRUE)[(ncskip+1):nctot,]
      yhat.test0 = yhat.test0 + null_intercept0
      #yhat.train.mean0 <- colMeans(yhat.train0)

      yhat.train1 <- matrix(cres1$trdraw, nrow=nctot, byrow=TRUE)[(ncskip+1):nctot,]
      yhat.train1 <- yhat.train1 + null_intercept1
      yhat.test1 = matrix(cres1$tedraw,nrow=nctot,byrow=TRUE)[(ncskip+1):nctot,]
      yhat.test1 = yhat.test1 + null_intercept1
      #yhat.train.mean1 <- colMeans(yhat.train1)
   }
 #  if(nrow(x.test)) {
 #     print("hello")
 #  }
   nsamps <- nrow(yhat.train0)
   npatients <- nrow(x.train)

   ControlEsts <- TrtEsts <- matrix(NA, nrow=ndpost, ncol=npatients)
   ControlEsts[,Trt==0] <- yhat.train0
   ControlEsts[,Trt==1] <- yhat.test0
   TrtEsts[,Trt==1] <- yhat.train1
   TrtEsts[,Trt==0] <- yhat.test1



   mix.prop = matrix(cres0$mixdraw,nrow=nctot,byrow=TRUE)[(ncskip+1):nctot,]
   locations = matrix(cres0$locdraw,nrow=nctot,byrow=TRUE)[(ncskip+1):nctot,]
   mass = cres0$Mdraw

   varcount = matrix(cres0$vcdraw,nrow=nctot,byrow=TRUE)[(ncskip+1):nctot,]


   Theta.test <- NULL

   if(scale=="log") {
       ThetaMat <- TrtEsts - ControlEsts
   }
   else if(scale=="time") {
      mgfs <- rowSums(exp(locations)*mix.prop)
      mgfs <- 1
      for(k in 1:npatients) {
        if(Trt[k]==1) {
           #ThetaMat[,k] <- mgfs*(exp(yhat.train[,k] - yhat.test[,k]))
           ThetaMat[,k] <- 0
        }
        else {
           #ThetaMat[,k] <- mgfs*(exp(yhat.test[,k] - yhat.train[,k]))
            ThetaMat[,k] <- 1
        }
      }
   }

   retval = list(
      call=match.call(),
      first.sigma=first.sigma,
      sigma=sigma,
      sigest=sigest,
      varcount=varcount,
      y = y.train,
      mix.prop=mix.prop,
      locations=locations,
      mass = mass,
      Theta = ThetaMat,
      Theta.test = Theta.test
   )
   class(retval) = 'indivaft'
   return(invisible(retval))
}

