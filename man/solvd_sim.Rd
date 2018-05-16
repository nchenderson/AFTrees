\name{solvd_sim}
\alias{solvd_sim}
\docType{data}
\title{
Simulated Data Set
}
\description{
Simulated Data Set to have the same variables as in the SOLVD Treatment and Prevention Trials Data.
Note that this is just a simulated data set and does not contain any data from the SOLVD trials.
}
\usage{data("solvd_sim")}
\format{
  A data frame with 6249 observations on the following 22 variables.
  \describe{
    \item{\code{ttodthorchfhosp}}{a numeric vector}
    \item{\code{trt}}{a factor with levels \code{Placebo} \code{treatment}}
    \item{\code{study}}{a factor with levels \code{SOLVD Prevention} \code{SOLVD Treatment}}
    \item{\code{deathorchfhosp}}{a numeric vector}
    \item{\code{age}}{a numeric vector}
    \item{\code{smoke}}{a factor with levels \code{Current} \code{Former} \code{Never}}
    \item{\code{diabet}}{a numeric vector}
    \item{\code{lvef}}{a numeric vector}
    \item{\code{himi}}{a numeric vector}
    \item{\code{nyha}}{a numeric vector}
    \item{\code{depedema}}{a numeric vector}
    \item{\code{crackles}}{a numeric vector}
    \item{\code{weight}}{a numeric vector}
    \item{\code{anydiurbl}}{a numeric vector}
    \item{\code{avsys}}{a numeric vector}
    \item{\code{avdia}}{a numeric vector}
    \item{\code{beat}}{a numeric vector}
    \item{\code{creatinine}}{a numeric vector}
    \item{\code{sodium}}{a numeric vector}
    \item{\code{copd}}{a numeric vector}
    \item{\code{histk}}{a numeric vector}
    \item{\code{gender}}{a factor with levels \code{female} \code{male}}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(solvd_sim)
head(solvd_sim)
}
\keyword{datasets}
