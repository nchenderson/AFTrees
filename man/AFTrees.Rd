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
   ntree=200,
   ndpost=1000, nskip=100,
   printevery=100, keepevery=1, keeptrainfits=TRUE,
   usequants=TRUE, numcut=100, printcutoffs=0,
   verbose=TRUE)
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
   \code{makeind} is used to generate the dummies.
   \code{bart} will generate draws of \eqn{f(x)} for each \eqn{x}
   which is a row of x.train.
   }

   \item{y.train}{
    A vector follow-up times
   }

   \item{status}{
    A vector of censoring indicators; status = 1 indicates that the event was observed
    while status = 0 indicates the observation was right-censored.
   }
   \item{x.test}{
   Explanatory variables for test (out of sample) data.\cr
   Should have same structure as x.train.\cr
   \code{bart} will generate draws of \eqn{f(x)} for each \eqn{x} which is a row of x.test.
   }

   \item{sigest}{
   Rough estimate for the error standard deviation.
   If sigest=NA then this will be calculated.
   }

   \item{sigdf}{
   Degrees of freedom for error variance prior.
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
   k is the number of prior standard deviations \eqn{E(Y|x) = f(x)} is away from +/-.5.
   The response (y.train) is internally scaled to range from -.5 to .5.
   For binary y,
   k is the number of prior standard deviations \eqn{f(x)} is away from +/-3.
   In both cases, the bigger k is, the more conservative the fitting will be.
   }

   \item{power}{
   Power parameter for tree prior.
   }

   \item{base}{
   Base parameter for tree prior.
   }

   \item{nonparametric}{
   An indicator of whether to use a nonparametric vs. semi-parametric fit
   }
   \item{ntree}{
   The number of trees in the sum.
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
   A \dQuote{draw} will consist of values of the error standard deviation (\eqn{\sigma}{sigma})
   and \eqn{f^*(x)}{f*(x)}
   at \eqn{x} = rows from the train(optionally) and test data, where \eqn{f^*}{f*} denotes
   the current draw of \eqn{f}.
   }

   \item{keeptrainfits}{
   If true the draws of \eqn{f(x)} for \eqn{x} = rows of x.train are returned.
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
   Vector of two colors.  First color is used to plot the median of \eqn{f(x)} and
   the second color is used to plot the lower and upper quantiles.
   }

   \item{\dots}{
   Additional arguments passed on to plot.
   }
}
\details{
   BART is an Bayesian MCMC method.
   At each MCMC interation, we produce a draw from the joint posterior
   \eqn{(f,\sigma) | (x,y)}{(f,sigma) \| (x,y)} in the numeric \eqn{y} case
   and just \eqn{f} in the binary \eqn{y} case.

   Thus, unlike a lot of other modelling methods in R, we do not produce a single model object
   from which fits and summaries may be extracted.  The output consists of values
   \eqn{f^*(x)}{f*(x)} (and \eqn{\sigma^*}{sigma*} in the numeric case) where * denotes a particular draw.
   The \eqn{x} is either a row from the training data (x.train) or the test data (x.test).
}
\value{
   The \code{plot} method sets mfrow to c(1,2) and makes two plots.\cr
   The first plot is the sequence of kept draws of \eqn{\sigma}{sigma}
   including the burn-in draws.  Initially these draws will decline as BART finds fit
   and then level off when the MCMC has burnt in.\cr
   The second plot has \eqn{y} on the horizontal axis and posterior intervals for
   the corresponding \eqn{f(x)} on the vertical axis.


   \code{AFTree} returns a list assigned class \sQuote{aftree}.
   In the numeric \eqn{y} case, the list has components:

   \item{yhat.train}{
   A matrix with (ndpost/keepevery) rows and nrow(x.train) columns.
   Each row corresponds to a draw \eqn{f^*}{f*} from the posterior of \eqn{f}
   and each column corresponds to a row of x.train.
   The \eqn{(i,j)} value is \eqn{f^*(x)}{f*(x)} for the \eqn{i^{th}}{i\^th} kept draw of \eqn{f}
   and the \eqn{j^{th}}{j\^th} row of x.train.\cr
   Burn-in is dropped.
   }

   \item{yhat.test}{Same as yhat.train but now the x's are the rows of the test data.}
   \item{yhat.train.mean}{train data fits = mean of yhat.train columns.}
   \item{yhat.test.mean}{test data fits = mean of yhat.test columns.}
   \item{sigma}{post burn in draws of sigma, length = ndpost/keepevery.}
   \item{first.sigma}{burn-in draws of sigma.}
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

   In the binary \eqn{y} case, the returned list has the components
   yhat.train, yhat.test, and varcount as above.  In addition the list
   has a binaryOffset component giving the value used.

   Note that in the binary \eqn{y}, case yhat.train and yhat.test are
   \eqn{f(x)} + binaryOffset.  If you want draws of the probability
   \eqn{P(Y=1 | x)} you need to apply the normal cdf (\code{pnorm})
   to these values.
}
\author{
Nicholas Henderson
}
\seealso{
\code{\link{pdbart}}
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
plot(tgrid, .5*dnorm(tgrid, mean=-1,sd=.5) + .5*dnorm(tgrid,mean=1,sd=.5), type="n")
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
