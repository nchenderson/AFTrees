\name{SurvivalProb}
\alias{SurvivalProb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Survival Probabilities
}
\description{
Computes P(T > t|x) for a covariate vector x and desired values of t
}
\usage{
SurvivalProb(object, time.points=NULL, test.only=FALSE, train.only=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class indivAFT
}
\item{time.points}{
the point (or points) at which the survival probability is computed
}
}
\value{
A list with the following two components
  \item{Surv}{a vector or matrix of survival probabilities}
   \item{time.point}{the time point entered}
}
\author{
Nicholas Henderson
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

\dontrun{
## generate responses without censoring
a0 = muvec + rnorm(n)
a1 = muvec_trt + rnorm(n)
yy = c(a0, a1)

Trt <- c(rep(0,n), rep(1,n))
XX <- matrix(rep(age,2), ncol=1)

aftFit <- IndivAFT(XX,exp(yy),status=rep(1,2*n),Trt=Trt,ndpost=100,nskip=100,printevery=100)

ss <- SurvivalProb(aftFit, time.point=2)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
