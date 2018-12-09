\name{Cstat}
\alias{Cstat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
C-statistic
}
\description{
Uno's C-statistic that is adjusted for censoring
}
\usage{
Cstat(mydata, risk.score, tau)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{tau}{
A truncation time
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with the following two components
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
a=2
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
