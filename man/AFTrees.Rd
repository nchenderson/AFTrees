\name{AFTrees}
\title{Accelerated Failure Regression Trees}
\alias{AFTrees}
\alias{plot.aftree}
\description{
Fits a Nonparametric Accelerated Failure Time Model
}
\usage{
AFTrees(
   x.train, y.train, status, x.test=matrix(0.0,0,0),
   sigest=NA, sigdf=3, sigquant=.90,
   k=2.0,
   power=2.0, base=.95,
   nonparametric=TRUE,
   ntree=200, max.error.clust=200,
   ndpost=1000, nskip=100,
   printevery=100, keepevery=1, keeptrainfits=TRUE,
   keepburnin=FALSE, usequants=TRUE, numcut=100,
   printcutoffs=0, verbose=TRUE)
\method{plot}{aftree}(
   x,
   plquants=c(.05,.95), cols =c('blue','black'),
   \dots)
}
\arguments{
   \item{x.train}{
   Explanatory variables for training (in sample) data.\cr
   May be a matrix or a data frame,
   with (as usual) rows corresponding to observations and columns to variables.\cr
   If a variable is a factor in a data frame, it is replaced with dummies.
   Note that q dummies are created if q>2 and
   one dummy is created if q=2, where q is the number of levels of the factor.
   }
   \item{y.train}{
    A vector follow-up times.
   }
   \item{status}{
    A vector of event indicators; status = 1 indicates that the event was observed
    while status = 0 indicates the observation was right-censored.
   }
   \item{x.test}{
   Explanatory variables for test (out of sample) data.\cr
   Should have same structure as x.train.\cr
   \code{AFTrees} will generate draws of \eqn{m(x)} for each \eqn{x} which is a row of x.test.
   }
   \item{sigest}{
   Rough, preliminary estimate for the residual standard deviation. If sigest=NA then this will be calculated internally.
   }
   \item{sigdf}{
   Degrees of freedom for residual variance prior.
   }
   \item{sigquant}{
   The quantile of the prior that the rough estimate (see sigest) is placed at.
   The closer the quantile is to 1,
   the more aggresive the fit will be as you are putting more prior weight
   on error standard deviations (\eqn{\sigma}{sigma}) less than the rough estimate.
   Not used if y is binary.
   }
   \item{k}{
   For numeric y,
   k is the number of prior standard deviations \eqn{E(log T|x) = m(x)} is away from +/-.5.
   The response (y.train) is internally scaled to range from -.5 to .5.
   }
   \item{power}{
   Power parameter for tree prior.
   }
   \item{base}{
   Base parameter for tree prior.
   }
   \item{nonparametric}{
   An indicator of whether to use a nonparametric vs. semi-parametric model. A semi-parametric
   model assumes a parametric model for the residual distribution.
   }
   \item{ntree}{
   The number of trees used for the regression function \eqn{m(x)}.
   }
   \item{max.error.clust}{
     The maximum number of clusters (or mixture components) that the error distribution may have.
   }
   \item{ndpost}{
   The number of posterior draws after burn in, ndpost/keepevery will actually be returned.
   }
   \item{nskip}{
   Number of MCMC iterations to be treated as burn in.
   }
   \item{printevery}{
   As the MCMC runs, a message is printed every printevery draws.
   }
   \item{keepevery}{
   Every keepevery draw is kept to be returned to the user.\cr
   A \dQuote{draw} will consist of values of the residual standard deviation (\eqn{\sigma}{sigma})
   and \eqn{m^*(x)}{m*(x)}
   at \eqn{x} = rows from the train(optionally) and test data, where \eqn{m^*}{m*} denotes
   the current draw of \eqn{m}.
   }
   \item{keepburnin}{
   If true, we keep every \code{keepevery} draw from the burn in phase.
   }
   \item{keeptrainfits}{
   If true, the draws of \eqn{m(x)} for \eqn{x} = rows of x.train are returned.
   }
   \item{usequants}{
   Decision rules in the tree are of the form
   \eqn{x \le c}{x <= c}  vs. \eqn{x > c}
   for each variable corresponding to a column of x.train.
   usequants determines how the set of possible c is determined.
   If usequants is true, then the c are
   a subset of the values (xs[i]+xs[i+1])/2 where xs is
   unique sorted values obtained from
   the corresponding column of x.train.
   If usequants is false, the cutoffs are equally spaced across the
   range of values taken on by the corresponding
   column of x.train.
   }
   \item{numcut}{
   The number of possible values of c (see usequants).
   If a single number if given, this is used for all variables.
   Otherwise a vector with length equal to ncol(x.train) is required,
   where the \eqn{i^{th}}{i^th} element gives the number of c used for
   the \eqn{i^{th}}{i^th} variable in x.train.
   If usequants is false, numcut equally spaced cutoffs
   are used covering the range of values in the corresponding
   column of x.train.  If usequants is true, then  min(numcut, the number of unique values in the
   corresponding columns of x.train - 1) c values are used.
   }
   \item{printcutoffs}{
   The number of cutoff rules c to printed to screen before the MCMC is run.
   Give a single integer, the same value will be used for all variables.
   If 0, nothing is printed.
   }
   \item{verbose}{
   Logical, if FALSE supress printing.
   }
   \item{x}{
   Value returned by \code{AFTree} which contains the information to be plotted.
   }
   \item{plquants}{
   In the plots, beliefs about \eqn{m(x)} are indicated by plotting the
   posterior median and a lower and upper quantile.
   plquants is a double vector of length two giving the lower and upper quantiles.
   }
   \item{cols}{
   Vector of two colors. First color is used to plot the median of \eqn{m(x)} and
   the second color is used to plot the lower and upper quantiles.
   }
   \item{\dots}{
   Additional arguments passed on to plot.
   }
}
\details{
   The aftrees function fits an Accelerated Failure Time (AFT) model where
   the regression function is modeled with Bayesian Additive Regression
   Trees (BART) while the residual distribution is modeled nonparametrically.
   Computation is done via MCMC. At each MCMC interation, we produce a draw from the joint posterior
   \eqn{(m,f) | (x,y)}{(f,sigma) \| (x,y)}.
}
\value{
   \code{AFTree} returns a list assigned class \sQuote{aftree}.
   The list has the following components:

   \item{m.train}{
   A matrix with (ndpost/keepevery) rows and nrow(x.train) columns.
   Each row corresponds to a draw \eqn{m^*}{m*} from the posterior of the regression function \eqn{m}
   and each column corresponds to a row of x.train. The function \eqn{m} is the regression function. This represents
   expected log-survival time in an AFT model.
   The \eqn{(i,j)} value is \eqn{m^*(x)}{m*(x)} for the \eqn{i^{th}}{i\^th} kept draw of \eqn{m}
   and the \eqn{j^{th}}{j\^th} row of x.train.\cr
   Burn-in is dropped.
   }

   \item{m.test}{Same as m.train except that the x's are the rows of the test data.}
   \item{m.train.mean}{train data fits = mean of m.train columns.}
   \item{m.test.mean}{test data fits = mean of m.test columns.}
   \item{sigma}{post burn in draws of sigma, length = ndpost/keepevery.}
   \item{first.sigma}{burn-in draws of sigma.}
   \item{mix.prop}{posterior draws of the vector of mixture proportions in the DP-mixture residual distribution.}
   \item{locations}{posterior draws of the vector of mixture locations in the DP-mixture residual distribution.}
   \item{varcount}{a matrix with (ndpost/keepevery) rows and nrow(x.train) columns.
   Each row is for a draw. For each variable (corresponding to the columns),
   the total count of the number of times
   that variable is used in a tree decision rule (over all trees) is given.}

   \item{sigest}{
   The rough error standard deviation (\eqn{\sigma}{sigma}) used in the prior.
   }

   \item{y}{
   The input dependent vector of values for the dependent variable.\cr
   This is used in plot.bart.
   }
}
\author{
Nicholas Henderson
}
\seealso{
\code{\link{IndivAFT}}
}
\examples{

set.seed(53718)

### AFT example, with Normal-mixture error distribution
n <- 500
age <- runif(n)
age <- sort(age)
muvec <- muvec_trt <- rep(0, n)
for(k in 1:n) {
  muvec[k] <- 1/2 + pnorm(5*(age[k] - .5))
  muvec_trt[k] <- .35 + muvec[k]
}
cens_times <- runif(2*n, min=1, max=3)

W1 <- W2 <- rep(0, n)
zz <- sample(1:2, size=n, replace=TRUE)
n1 <- sum(zz==1)
n2 <- sum(zz==2)
W1[zz==1] <- rnorm(n1, mean=-1,sd=.5)
W1[zz==2] <- rnorm(n2, mean=1,sd=.5)
zz <- sample(1:2, size=n, replace=TRUE)
n1 <- sum(zz==1)
n2 <- sum(zz==2)
W2[zz==1] <- rnorm(n1, mean=-1,sd=.5)
W2[zz==2] <- rnorm(n2, mean=1,sd=.5)


#### generate log-survival times
a0 <- muvec + W1
a1 <- muvec_trt + W2
yy <- c(a0, a1)
status <- as.numeric(yy <= cens_times)


yy[status==0] <- cens_times[status==0]
Trt <- c(rep(0,n), rep(1,n))
XX <- cbind(rep(age,2), Trt)

system.time(aftFit <- AFTrees(XX, exp(yy),status=status,ndpost=500,nskip=500,printevery=100))

### Compute fitted expected log-survival for each treatment group
muhat <- colMeans(aftFit$m.train)

### Compare with true values

plot(age, muvec, type="n", ylim=c(.4,2))
lines(age, muvec)
lines(age, muvec_trt)
lines(age, muhat[Trt==0], col="red")
lines(age, muhat[Trt==1], col="red")


### plot true error density and posterior error densities
tgrid <- seq(-3, 3, length.out=200)
plot(tgrid, .5*dnorm(tgrid, mean=-1,sd=.5) + .5*dnorm(tgrid,mean=1,sd=.5), type="n",
    ylab="Density")
lines(tgrid, .5*dnorm(tgrid, mean=-1,sd=.5) + .5*dnorm(tgrid,mean=1,sd=.5), col="red",lwd=3)
for(k in 110:120) {
  FF <- dnorm(outer(tgrid, aftFit$locations[k,], FUN="-"), sd=aftFit$sigma[k])
  ff <- aftFit$mix.prop[k,]
  fhat <- FF \%*\% ff
  lines(tgrid, fhat)
}


}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
