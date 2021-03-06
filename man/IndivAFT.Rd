\name{IndivAFT}
\alias{IndivAFT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Individualized Treatment Effects
}
\description{
Estimation of Individualized Treatment Effects with Nonparametric or
Semiparametric Accelerated Failure Time Models
}
\usage{
IndivAFT(
   x.train, y.train, status, Trt, x.test=NULL,
   sigest=NA, sigdf=3, sigquant=.5,
   k=2.0,power=2.0, base=.95, nonparametric=TRUE,
   ntree=200,
   ndpost=1000, nskip=100,
   printevery=100, keepevery=1, keeptrainfits=TRUE,
   usequants=TRUE, numcut=100, printcutoffs=0,
   verbose=TRUE, scale="log")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.train}{
   Explanatory variables for training (in sample) data.\cr
   May be a matrix or a data frame,
   with (as usual) rows corresponding to observations and columns to variables.\cr
   If a variable is a factor in a data frame, it is replaced with dummies.
   Note that q dummies are created if q>2 and
   one dummy is created if q=2, where q is the number of levels of the factor.
   \code{makeind} is used to generate the dummies.
   }


   \item{y.train}{
   A vector follow-up times
   }

   \item{status}{
   A vector of censoring indicators; status = 1 indicates that the event was observed
   while status = 0 indicates that the observation was right-censored.
   }
   \item{Trt}{
   Treatment assignment
   }
   \item{x.test}{
   Explanatory variables for test data.
   }
   \item{sigest}{
   The prior for the error variance  (\eqn{\sigma^2}{sigma\^2}) is inverted chi-squared
   (the standard conditionally conjugate prior).
   The prior is specified by choosing the degrees of freedom, a rough estimate of the
   corresponding standard deviation and a quantile to put this rough estimate at.
   If sigest=NA then the rough estimate will be the usual least squares estimator.
   Otherwise the supplied value will be used.
   Not used if y is binary.
   }

   \item{sigdf}{
   Degrees of freedom for error variance prior.
   Not used if y is binary.
   }

   \item{sigquant}{
   The quantile of the prior that the rough estimate (see sigest) is placed at.
   The closer the quantile is to 1,
   the more aggresive the fit will be as you are putting more prior weight
   on error standard deviations (\eqn{\sigma}{sigma}) less than the rough estimate.
   Not used if y is binary.
   }

   \item{k}{
   Shrinkage parameter
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
   \item{scale}{
   if equal to log, ITEs are returned on the log-difference scale. Otherwise,
   ITEs are returned as ratios of expected failure times.
   }
}
\value{

   \code{IndivAFT} returns a list of class ()
   with (at least) the following components:

}
\author{
Nicholas Henderson
}
\seealso{
\code{\link{AFTrees}}
}
\examples{

n <- 500
age <- runif(n)
age <- sort(age)
muvec <- muvec_trt <- rep(0, n)
for(k in 1:n) {
  muvec[k] <- 1/2 + pnorm(5*(age[k] - .5))
  muvec_trt[k] <- .35 + muvec[k]
}

## generate responses without censoring
a0 = muvec + rnorm(n)
a1 = muvec_trt + rnorm(n)
yy = c(a0, a1)

Trt <- c(rep(0,n), rep(1,n))
XX <- matrix(rep(age,2), ncol=1)

aftFit <- IndivAFT(XX,exp(yy),status=rep(1,2*n),Trt=Trt,ndpost=100,nskip=100,printevery=100)

### get posterior means of treatment effect:
post.means <- colMeans(aftFit$Theta)

#plot(aftFit)

}
\keyword{ nonparametric }
\keyword{ trees }
