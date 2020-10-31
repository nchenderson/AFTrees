\name{Cstat}
\alias{Cstat}
\title{
C-statistic
}
\description{
Uno's C-statistic that is adjusted for censoring
}
\usage{
Cstat(times, status, risk.score, time.stop = NULL)
}
\arguments{
  \item{times}{
    A vector follow-up times
}
\item{status}{
 A vector of censoring indicators; status = 1 indicates that the event was observed
    while status = 0 indicates the observation was right-censored.
}
  \item{risk.score}{
A vector of risk scores. This should have the same length as both \code{times} and \code{status}
}
  \item{time.stop}{
A truncation time point which provides an upper limit for the time range (0, \code{time.stop}) under consideration.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with the following two components
  \item{Chat}{The c-statistic}
  \item{comp2 }{The estimated values of the censoring survival function}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
library(TH.data)
library(survival)
data("GBSG2")
ind <- sample(1:nrow(GBSG2), size=500)
GB.train <- GBSG2[ind,]
GB.test <- GBSG2[-ind,]

#########################################################
## Compute c-statistics on the test set using only training data
x.train <- model.matrix(~ horTh + age + menostat + tsize + tgrade + pnodes
                          + progrec + estrec + 0, data=GB.train)
x.test <- model.matrix(~ horTh + age + menostat + tsize + tgrade + pnodes
                       + progrec + estrec + 0, data=GB.test)

y.train <- as.double(GB.train$time)
y.test <- as.double(GB.test$time)

status.train <- GB.train$cens
status.test <- GB.test$cens

cox.train <- coxph(Surv(time, cens) ~ horTh + age + menostat + tsize + tgrade
                    + pnodes + progrec + estrec, x=TRUE, y=TRUE, method="breslow",
                    data=GB.train)

cox.pred <- predict(cox.train, newdata=GB.test)
aft.obj <- AFTrees(x.train, y.train, status.train, x.test=x.test, ndpost=500)
rs.bart <- RiskScore(aft.obj)$test

## Compare c-statistics
tau <- quantile(y.test, probs=0.95)
Cstat(y.test, status.test, rs.bart, time.stop=tau)$Chat
Cstat(y.test, status.test, cox.pred, time.stop=tau)$Chat

###################################################################3
## Compute in-sample c-statistics for full data sets.
x.full <- model.matrix(~ horTh + age + menostat + tsize + tgrade + pnodes
                        + progrec + estrec + 0, data=GBSG2)
y.full <- as.double(GBSG2$time)
status.full <- GBSG2$cens

aft.obj <- AFTrees(x.full, y.full, status.full, ndpost=500)

rs.bart <- RiskScore(aft.obj)$train

cox.full <- coxph(Surv(time, cens) ~ horTh + age + menostat + tsize + tgrade
                  + pnodes + progrec + estrec, x=TRUE, y=TRUE, method="breslow",
                  data=GBSG2)

cox.pred <- predict(cox.full, newdata=GBSG2)

## Compare in-sample c-statistics
tau <- quantile(y.full, probs=0.95)
Cstat(y.full, status.full, rs.bart, time.stop=tau)$Chat
Cstat(y.full, status.full, cox.pred, time.stop=tau)$Chat
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
