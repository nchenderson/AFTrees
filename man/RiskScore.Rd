\name{RiskScore}
\alias{RiskScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Individual-specific risk scores
}
\description{
From a fitted aftree object, computes covariate-specific risk scores.
}
\usage{
RiskScore(object, type=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
 A fitted aftree object
}
}
\details{
The risk score for individual i is defined as the exponential of the expected
log-failure time. Specifically, the risk score is defined as exp[ E( log T | x[i] )  ].
}
\value{
A list with two components
  \item{train}{Risk scores for those in the training set.}
  \item{test}{Risk scores for those in the test set (if the fitted object
  contains fitted values for a test set).}
}
\author{
Nicholas Henderson
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{AFTrees}}.
}
\examples{
data(solvd_sim)

## Compute baseline risk scores as a function of ejection fraction
control.group <- solvd_sim$trt=="Placebo"
x.train <- matrix(solvd_sim$lvef[control.group], ncol=1)
y.train <- solvd_sim$ttodthorchfhosp[control.group]
status <- solvd_sim$deathorchfhosp[control.group]

aft.obj <- AFTrees(x.train, y.train, status)
rs <- RiskScore(aft.obj)

plot(x.train[,1], rs$train)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
