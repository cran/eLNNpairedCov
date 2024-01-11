\name{genSimDat}
\alias{genSimDat}
\title{
Generate Simulated Data
}
\description{
Generate a simulated dataset from a mixture of Bayesian hierarchical models
with two covariates: age and sex.
}
\usage{
genSimDat(G, n, psi, t_pi, m.age = 50, sd.age = 5, p.female = 0.5)
}
\arguments{
  \item{G}{
	  integer. Number of probes.
}
  \item{n}{
	  integer. Number of samples.
}
  \item{psi}{
	  numeric. A vector of model hyper-parameters with elements
\eqn{\alpha_1}, \eqn{\beta_1}, \eqn{k_1}, 
\eqn{\eta_{1,intercept}}, \eqn{\eta_{1,age}}, \eqn{\eta_{1,sex}},
\eqn{\alpha_2}, \eqn{\beta_2}, \eqn{k_2}, 
\eqn{\eta_{2,intercept}}, \eqn{\eta_{2,age}}, \eqn{\eta_{2,sex}},
\eqn{\alpha_3}, \eqn{\beta_3}, \eqn{k_3}, 
\eqn{\eta_{3,intercept}}, \eqn{\eta_{3,age}}, \eqn{\eta_{3,sex}}.
}
  \item{t_pi}{
	  numeric. A vector of mixture proportions: \eqn{\pi_1} (proportion
	  for probes over-expressed in cases); 
	  \eqn{\pi_2} (proportion for probes under-expressed in cases).
}
  \item{m.age}{
	  numeric. mean age.
}
  \item{sd.age}{
	  numeric. standard deviation of age.
}
  \item{p.female}{
	  numeric. proportion of females.
}
}

\value{
	An ExpressionSet object.
}
\references{
Zhang Y, Liu W, Qiu W.
A model-based clustering via mixture of hierarchical models with covariate adjustment for detecting differentially expressed genes from paired design.
\emph{BMC Bioinformatics} 24, 423 (2023)
}
\author{
	Yixin Zhang \email{zhyl133@gmail.com},
	Wei Liu \email{liuwei@mathstat.yorku.ca},
	Weiliang Qiu \email{weiliang.qiu@sanofi.com}
}
\note{
	Age will be mean-centered and scaled so that it will
	have mean zero and variance one.
}

\examples{

set.seed(1234567)

true.psi = c(2, 1, 0.8,
             0.1, -0.01, -0.1,
             2, 1, 0.8,
             -0.1, -0.01, -0.1,
             2, 1, 0.8,
             -0.01, -0.1)
names(true.psi)=c("alpha1", "beta1", "k1",
                  "eta1.intercept", "eta1.Age", "eta1.Sex",
                  "alpha2", "beta2", "k2",
                  "eta2.intercept", "eta2.Age", "eta2.Sex",
                  "alpha3", "beta3", "k3",
                  "eta3.Age", "eta3.Sex")
true.pi=c(0.1, 0.1)
names(true.pi)=c("pi.OE", "pi.UE")
par.true=c(true.pi, true.psi)

esDiff = genSimDat(G = 1000, 
	       n = 20, 
	       psi = true.psi, 
	       t_pi = true.pi, 
	       m.age = 0, # scaled age 
	       sd.age = 1, # scaled age
	       p.female = 0.5)

print(esDiff)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }% use one of  RShowDoc("KEYWORDS")
