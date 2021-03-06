AFTrees = function(
   x.train, y.train, status, x.test=matrix(0.0,0,0),
   sigest=NA, sigdf=3, sigquant=.90,
   k=2.0,
   power=2.0, base=.95,
   nonparametric=TRUE,
   ntree=200, max.error.clust=200,
   ndpost=1000,  nskip=100,
   printevery=100, keepevery=1, keeptrainfits=TRUE,
   keepburnin=FALSE, usequants=TRUE, numcut=100, 
   printcutoffs=0, verbose=TRUE) {

   if(is.vector(x.train) | is.factor(x.train)) x.train = data.frame(x=x.train)
   if(is.vector(x.test) | is.factor(x.test)) x.test = data.frame(x=x.test)

   if(is.data.frame(x.train)) {
      if(nrow(x.test)) {
         if(!is.data.frame(x.test)) {
             stop('x.train is a data frame so x.test must be also')
         }
	     xtemp = rbind(x.train,x.test)
	     xmtemp = makeind(xtemp)
	     x.train = xmtemp[1:nrow(x.train),,drop=FALSE]
	     x.test = xmtemp[nrow(x.train) + 1:nrow(x.test),,drop=FALSE]
      } else {
         x.train = makeind(x.train)
      }
   }

   ## check input arguments:
   if((!is.matrix(x.train)) || (typeof(x.train)!="double")) stop("argument x.train must be a double matrix")
   if((!is.matrix(x.test)) || (typeof(x.test)!="double")) stop("argument x.test must be a double matrix")
   if((!is.vector(y.train)) || (typeof(y.train)!="double")) stop("argument y.train must be a double vector")
   if(nrow(x.train) != length(y.train)) stop("number of rows in x.train must equal length of y.train")
   if((nrow(x.test) >0) && (ncol(x.test)!=ncol(x.train))) stop("input x.test must have the same number of columns as x.train")
   if((!is.na(sigest)) && (typeof(sigest)!="double")) stop("input sigest must be double")
   if((!is.na(sigest)) && (sigest<0.0)) stop("input sigest must be positive")
   if((mode(sigdf)!="numeric") || (sigdf<0)) stop("input sigdf must be a positive number")
   if((mode(printevery)!="numeric") || (printevery<0)) stop("input printevery must be a positive number")
   if((mode(keepevery)!="numeric") || (keepevery<0)) stop("input keepevery must be a positive number")
   if((mode(sigquant)!="numeric") || (sigquant<0)) stop("input sigquant must be a positive number")
   if((mode(ntree)!="numeric") || (ntree<0)) stop("input ntree must be a positive number")
   if((mode(ndpost)!="numeric") || (ndpost<0)) stop("input ndpost must be a positive number")
   if((mode(nskip)!="numeric") || (nskip<0)) stop("input nskip must be a positive number")
   if((mode(k)!="numeric") || (k<0)) stop("input k must be a positive number")
   if(mode(numcut)!="numeric") stop("input numcut must be a numeric vector")
   if(length(numcut)==1) numcut = rep(numcut,ncol(x.train))
   if(length(numcut) != ncol(x.train)) stop("length of numcut must equal number of columns of x.train")
   numcut = as.integer(numcut)
   if(min(numcut)<1) stop("numcut must be >= 1")
   if(typeof(usequants)  != "logical") stop("input usequants must a logical variable")
   if(typeof(keeptrainfits)  != "logical") stop("input keeptrainfits must a logical variable")
   if(typeof(verbose)  != "logical") stop("input verbose must a logical variable")
   if(mode(printcutoffs)  != "numeric") stop("input printcutoffs must be numeric")
   printcutoffs = as.integer(printcutoffs)
   if(printcutoffs <0) stop("input printcutoffs must be >=0")
   if(power <= 0) stop("power must be positive")
   if(base <= 0) stop("base must be positive")
   ############################################################################## 

   null_lm <- survreg(Surv(y.train, status) ~ 1, dist="lognormal")
   null_sig <- null_lm$scale
   null_intercept <- null_lm$coefficients

   y_centered_log <- log(y.train) - null_intercept
   y_unobs <- y.train*exp(.02*(1-status)*null_sig)   ### rough, estimate of unobserved survival times

   #y.train.log = log(y.train)  ### transform to log-survival times
   rgy = range(log(y_unobs))
   zeta <- rgy[2] - rgy[1]
   #aa = (rgy[1] + rgy[2])/(2*(rgy[2] - rgy[1]))
   #bb = 1/(rgy[2] - rgy[1])
   #y = bb*log(y.train) - aa   ## transformed log-responses
   #y = -.5 + (y.train.log-rgy[1])/(rgy[2]-rgy[1])
   #y = y.train

   # if sigest=NA, fit a lm to training data to get the value of sigest...
   # sigest is on the scale of the transformed y, so we do the lm after the scaling above...

   # Need to include the survival package
   if (is.na(sigest)) {
       tmplm <- survreg(Surv(y.train, status) ~ x.train, dist="lognormal")
       sigest <- tmplm$scale
   } else {
       sigest = sigest #put input sigma estimate on transformed scale
   }

   ncskip = floor(nskip/keepevery)
   ncpost = floor(ndpost/keepevery)
   nctot = ncskip + ncpost
   totnd = keepevery*nctot

   ## put sigquant argument here
   kappa <- FindKappa(q=sigquant, sigsq.hat=sigest*sigest, nu=sigdf)

   #nclust = 200
   nclust = max.error.clust
   npind = ifelse(nonparametric, 1, 0)
   cres = .C('mbart',as.integer(nrow(x.train)), as.integer(ncol(x.train)), as.integer(nrow(x.test)),
                   as.double(x.train), as.double(y_centered_log),
                   as.double(x.test), as.integer(status),
                   as.double(sigest),   as.integer(sigdf), as.double(sigquant),
                   as.double(k),as.double(power), as.double(base), as.integer(ntree), as.integer(totnd),
                   as.integer(printevery), as.integer(keepevery),  as.integer(keeptrainfits),
                   as.integer(numcut), as.integer(usequants), as.integer(printcutoffs),
		           as.integer(verbose), sdraw=double(nctot),
                   trdraw=double(nrow(x.train)*nctot),
                   tedraw=double(nrow(x.test)*nctot),
                   vcdraw=integer(ncol(x.train)*nctot), as.integer(nclust),
                   mixdraw=double(nclust*nctot),
                   locdraw=double(nclust*nctot),
                   Mdraw=double(nctot), npind=as.integer(npind),as.double(kappa),
                   as.double(zeta))

   # now read in the results...
   ## look at this: sigma is multiplied by (rgy[2] - rgy[1])
   sigma = cres$sdraw

   first.sigma = NULL
   first.m.train = NULL
   first.m.test = NULL
   
   if(keepburnin) {
       first.sigma = sigma[1:ncskip] # we often want the sigma draws
   }
   sigma = sigma[ncskip+(1:ncpost)]

   yhat.train = yhat.test = yhat.train.mean = yhat.test.mean = NULL
   varcount = NULL

   if(keeptrainfits) {
       tmp <- matrix(cres$trdraw, nrow=nctot,byrow=TRUE) ## this should be nctot x n
       if(keepburnin) {
          first.m.train = tmp[1:ncskip,] + null_intercept
       }
       yhat.train = tmp[(ncskip+1):nctot,] + null_intercept
       yhat.train.mean <- colMeans(yhat.train)
   }
   if (nrow(x.test)) {
       tmp <- matrix(cres$tedraw, nrow=nctot, byrow=TRUE) ## 
       if(keepburnin) {
          first.m.test <- tmp[1:ncskip,] + null_intercept
       }
       yhat.test = tmp[(ncskip+1):nctot,] + null_intercept
       yhat.test.mean <- colMeans(yhat.test)
   }
   mix.prop = matrix(cres$mixdraw,nrow=nctot,byrow=TRUE)[(ncskip+1):nctot,]
   locations = matrix(cres$locdraw, nrow=nctot, byrow=TRUE)[(ncskip+1):nctot,]
   mass = cres$Mdraw

   varcount = matrix(cres$vcdraw,nrow=nctot,byrow=TRUE)[(ncskip+1):nctot,]

   retval = list(
      call=match.call(),
      first.sigma=first.sigma,
      first.m.train=first.m.train,
      first.m.test=first.m.test,
      sigma=sigma,
      sigest=sigest,
      m.train=yhat.train,
      m.train.mean=yhat.train.mean,
      m.test=yhat.test,
      m.test.mean=yhat.test.mean,
      varcount=varcount,
      times = y.train,
      status = status,
      mix.prop=mix.prop,
      locations=locations,
      mass = mass
   )
   class(retval) = 'aftree'
   return(invisible(retval))
}
