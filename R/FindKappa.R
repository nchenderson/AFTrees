dinvchisq <- function(x, nu) {
    ans <- (nu/2)*log(nu/2) - lgamma(nu/2) - (nu/2 + 1)*log(x) - nu/(2*x)
    return(exp(ans))
}

pchiconv <- function(p, nu) {
    ### CDF for the random variable defined as
    ###  nu/X + Z  where X is chi_sq with nu degrees of freedom and Z ~ N(1,1) 
    ff <- function(v,x) {
         ans <- pnorm(x - v, mean=1,sd=1)*dinvchisq(v, nu=nu)
         return(ans)
    }
    tau <- length(p)
    a <- rep(0, tau)
    for(k in 1:tau) {
        a[k] <- integrate(ff, x=p[k], lower=0, upper=Inf)$value
    }
    return(a)
}

qchiconv <- function(q, nu) {
    gg <- function(x) {
         pchiconv(x, nu=nu) - q
    }
    a <- uniroot(gg, interval=c(0,200))
    return(a$root)
}

FindKappa <- function(q, sigsq.hat, nu) {
    ### don't compute q directly if equals .9 or .99
    if(q < .905 & q > .895) {
         Q <- 6.1423
         print('hello')
    }
    else if(q < .995 & q > .985) {
         Q <- 27.127
    }
    else if(q < .8 & q > .7) {
         print('hell75')
         Q <- 3.489662
    }
    else if(q < .55 & q > .45) {
         print('hell5')
         Q <- 2.287994
    }
    else if(q < .3 & q > .2) {
         print('hell25')
         Q <- 1.735869
    }
    else {
         Q <- qchiconv(q, nu=nu)
    }
    Kap <- sigsq.hat/Q
    return(Kap)
}

## sigquant = .1 case
#qchiconv(.1, 3)
#[1] 1.447128
#### How to compute using cubature package and 2-dimensional integration
#dinvchisq <- function(x, nu) {
#    l1 <- (nu/2)*log(nu) - (nu/2)*log(2) - lgamma(nu/2) 
#    l2 <- (nu/2 + 1)*log(x) + nu/(2*x)
#    ans <- exp(l1 - l2)
#    return(ans)
#} 

#pchiconv <- function(p, nu) {
  ### CDF for the random variable defined as
  ###  nu/X + Z  where X is chi_sq with nu degrees of freedom and Z ~ N(1,1) 
#  ff <- function(v,x) {
#      ss <- 1/sqrt(2*(v[2]+1))
#      ans <- pnorm(x - v[1], mean=1,sd=ss)*dinvchisq(v[1], nu=nu)*dgamma(v[2], shape=2, rate=.1)
#      return(ans)
#  }
#  tau <- length(p)
#  a <- rep(0, tau)
#  for(k in 1:tau) {
    #a[k] <- integrate(ff, x=p[k], lower=0, upper=Inf)$value
#    a[k] <- adaptIntegrate(ff, lowerLimit=c(0,0), upperLimit=c(40,600), tol=1e-8, x=p[k])$integral
#  }
#  return(a)
#}
#
#qchiconv <- function(q, nu) {
#  gg <- function(x) {
#    pchiconv(x, nu=nu) - q
#  }
#  a <- uniroot(gg, interval=c(0,100))
#  return(a$root)
#}

#qchiconv(.99, 3)

### 6.142334  - for 0.9
### 27.12703  - for 0.99
### 3.489662   - for 0.75
### 2.287994  - for 0.5
 


