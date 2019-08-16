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
    \item{\code{ttodthorchfhosp}}{Time to death or hospitalization.}
    \item{\code{trt}}{Treatment assignment. Factor with levels: \code{Placebo} and \code{treatment}}
    \item{\code{study}}{a factor representing one of the two SOLVD studies.
    This has levels \code{SOLVD Prevention} \code{SOLVD Treatment}.}
    \item{\code{deathorchfhosp}}{Event indicator. A binary numeric vector.}
    \item{\code{age}}{Age. A numeric vector.}
    \item{\code{smoke}}{Smoking history. A factor with levels: \code{Current}, \code{Former}, and \code{Never}.}
    \item{\code{diabet}}{Diabetic indicator. A binary numeric vector.}
    \item{\code{lvef}}{Ejection fraction. A numeric vector.}
    \item{\code{himi}}{History of Myocardial Infarction. A numeric vector.}
    \item{\code{nyha}}{New York heart association score. A numeric vector.}
    \item{\code{depedema}}{Dependent edema. A binary numeric vector.}
    \item{\code{crackles}}{Lung crackles. A binary numeric vector.}
    \item{\code{weight}}{Body weight. A numeric vector.}
    \item{\code{anydiurbl}}{Any use of a diuretic. A binary numeric vector.}
    \item{\code{avsys}}{Systolic blood pressure. A numeric vector.}
    \item{\code{avdia}}{Diastolic blood pressure. A numeric vector.}
    \item{\code{beat}}{Beats per minute. A numeric vector.}
    \item{\code{creatinine}}{Creatinine levels. A numeric vector.}
    \item{\code{sodium}}{Sodium. A numeric vector.}
    \item{\code{copd}}{Chronic Obstructive Pulmonary Disease. A binary numeric vector.}
    \item{\code{histk}}{History of stroke. A binary numeric vector.}
    \item{\code{gender}}{A factor with levels: \code{female} and \code{male}.}
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
summary(solvd_sim$ttodthorchfhosp)
}
\keyword{datasets}
