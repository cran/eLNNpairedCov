\name{esDiff}
\alias{esDiff}
\docType{data}
\title{
	An ExpressionSet Object Storing a Simulated Data
}
\description{
	An ExpressionSet object storing a simulated data of log2 difference of expression levels with 1000 probes, 20 subjects, and 2 covariates.
}
\usage{data("esDiff")}
\details{
	This dataset was generated from the mixture of 3-component
	Bayesian hierarchical models. For true parameters, please
	refer to the manual for the R function \code{\link{genSimDat}}.
}

\examples{
data(esDiff)

print(esDiff)

}
\keyword{datasets}
