IndivAFT = function(
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
   Trt.alt <- 1 - Trt ## counterfactual treatment
   if(is.null(x.test)) {
      x.test <- cbind(x.train, Trt.alt)
      x.train <- cbind(x.train, Trt)
      n.test <- 0
   } else {
      n.test <- nrow(x.test)
      t1 <- cbind(x.train, Trt.alt)
      t2 <- rbind(cbind(x.test, rep(1, n.test)), cbind(x.test, rep(0,n.test)))
      #print(dim(x.train))
      #print(length(Trt.alt))
      x.train <- cbind(x.train, Trt)
      x.test <- rbind(t1, t2)
   }

   #check input arguments:
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

   null_lm <- survreg(Surv(y.train, status) ~ 1, dist="lognormal")
   null_sig <- null_lm$scale
   null_intercept <- null_lm$coefficients

   y_centered_log <- log(y.train) - null_intercept
   imr <- dnorm(log(y.train), mean=null_intercept, sd=null_sig)/pnorm(log(y.train), mean=null_intercept, sd=null_sig, lower.tail=FALSE)
   #y_unobs <- y.train*exp(.02*(1-status)*null_sig)   ### rough, estimate of unobserved survival times
   y_unobs <- exp( status*log(y.train) + (1 - status)*null_sig + (1 - status)*imr)

   #y.train.log = log(y.train)  ### transform to log-survival times
   rgy = range(log(y_unobs))
   zeta_tmp <- rgy[2] - rgy[1]
   zeta <- 4*null_sig
   #print(c(zeta,zeta_tmp))
   #aa = (rgy[1] + rgy[2])/(2*(rgy[2] - rgy[1]))


   # Need to include the survival package
   if (is.na(sigest)) {
       #templm = lm(y~x.train)
       #sigest = summary(templm)$sigma
       tmplm <- survreg(Surv(y.train, status) ~ x.train, dist="lognormal")
       sigest <- tmplm$scale
   } else {
       sigest = sigest #put input sigma estimate on transformed scale
   }

   ncskip = floor(nskip/keepevery)
   ncpost = floor(ndpost/keepevery)
   nctot = ncskip + ncpost
   totnd = keepevery*nctot

   nclust = 200
   kappa <- FindKappa(q=sigquant, sigsq.hat=sqrt(sigest), nu=sigdf)

   npind = ifelse(nonparametric, 1, 0)
   cres = .C('mbart',as.integer(nrow(x.train)), as.integer(ncol(x.train)), as.integer(nrow(x.test)),
                   as.double(x.train), as.double(y_centered_log),
                   as.double(x.test), as.integer(status),
                   as.double(sigest),   as.integer(sigdf), as.double(sigquant),
                   as.double(k),
		   as.double(power), as.double(base),
		   as.integer(ntree),      as.integer(totnd),
                   as.integer(printevery), as.integer(keepevery),  as.integer(keeptrainfits),
                   as.integer(numcut), as.integer(usequants), as.integer(printcutoffs),
		   as.integer(verbose),
                   sdraw=double(nctot),
                   trdraw=double(nrow(x.train)*nctot),
                   tedraw=double(nrow(x.test)*nctot),
                   vcdraw=integer(ncol(x.train)*nctot), as.integer(nclust),
                   mixdraw=double(nclust*nctot),
                   locdraw=double(nclust*nctot),
                   Mdraw=double(nctot),npind=as.integer(npind),as.double(kappa),
                   as.double(zeta))
   # now read in the results...
       ## look at this: sigma is multiplied by (rgy[2] - rgy[1])
       sigma = cres$sdraw
      # sigma = cres$sdraw
       first.sigma = sigma[1:ncskip] # we often want the sigma draws
       sigma = sigma[ncskip+(1:ncpost)]

       # put sigest on the original y scale for output purposes
       sigest = sigest

   m.train = m.test = m.train.mean = m.test.mean = NULL
   varcount = NULL

   if (keeptrainfits) {
      m.train = matrix(cres$trdraw,nrow=nctot,byrow=T)[(ncskip+1):nctot,]
      m.train = m.train + null_intercept
      m.train.mean <- colMeans(m.train)
   }
   m.test = matrix(cres$tedraw,nrow=nctot,byrow=T)[(ncskip+1):nctot,]
   m.test = m.test + null_intercept

   mix.prop = matrix(cres$mixdraw,nrow=nctot,byrow=T)[(ncskip+1):nctot,]
   locations = matrix(cres$locdraw,nrow=nctot,byrow=T)[(ncskip+1):nctot,]
   mass = cres$Mdraw

   varcount = matrix(cres$vcdraw,nrow=nctot,byrow=T)[(ncskip+1):nctot,]

   npatients <- nrow(x.train)
   nsamps <- nrow(m.train)
   ThetaMat <- matrix(0, nrow=nsamps, ncol=npatients)
   Theta.test <- NULL

   if(scale=="log") {
       for(k in 1:npatients) {
           if(Trt[k]==1) {
                ThetaMat[,k] <- m.train[,k] - m.test[,k]
           }
           else {
               ThetaMat[,k] <- m.test[,k] - m.train[,k]
           }
       }
       if(n.test > 0) {
           Theta.test <- matrix(0, nrow=nsamps, ncol=n.test)
           for(h in 1:n.test) {
              Theta.test[,h] <- m.test[,npatients + h] - m.test[,npatients + h + n.test]
           }
       }
   }
   else if(scale=="time") {
      mgfs <- rowSums(exp(locations)*mix.prop)
      mgfs <- 1
      for(k in 1:npatients) {
        if(Trt[k]==1) {
           ThetaMat[,k] <- mgfs*(exp(m.train[,k] - m.test[,k]))
        }
        else {
           ThetaMat[,k] <- mgfs*(exp(m.test[,k] - m.train[,k]))
        }
      }
   }

   retval = list(
      call=match.call(),
      first.sigma=first.sigma,
      sigma=sigma,
      sigest=sigest,
      m.train=m.train,
      m.train.mean=m.train.mean,
      m.test=m.test,
      m.test.mean=m.test.mean,
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
