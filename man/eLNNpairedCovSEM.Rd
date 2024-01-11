\name{eLNNpairedCovSEM}
\alias{eLNNpairedCovSEM}
\title{
	Model-Based Clustering for Paired Data Adjusting for Covariates Using Simulated Annealing Modified EM
}
\description{
Model-based clustering based on extended log-normal normal model for paired data adjusting for covariates.
}
\usage{
eLNNpairedCovSEM(
  EsetDiff,
  fmla = ~Age + Sex,
  probeID.var = "probeid",
  gene.var = "gene",
  chr.var = "chr",
  scaleFlag = TRUE,
  Maxiter =10,
  maxIT = 10,
  b=c(2,2,2),
  converge_threshold = 1e-3,
  optimMethod = "L-BFGS-B",
  bound.alpha = c(0.001, 6),
  bound.beta = c(0.001, 6),
  bound.k = c(0.001, 0.9999),
  bound.eta = c(-10, 10),
  mc.cores = 1,
  temp0 = 2, 
  r_cool=0.9,
  verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{EsetDiff}{
	  An ExpressionSet object storing the log2 difference
	  between post-treatment and pre-treatment.
}
  \item{fmla}{
A formula without outcome variable.
}
  \item{probeID.var}{
	  character. Indicates the probe id.
}
  \item{gene.var}{
	  character. Indicates the gene symbol.
}
  \item{chr.var}{
	  character. Indicates the chromosome.
}
  \item{scaleFlag}{
  logical. Indicating if rows (probes) need to be scaled (but not centered).
  }
  \item{Maxiter}{
	  integer. 
	 The max allowed number of iterations for EM algorithm. Default value is maxRT = 100. 
}
  \item{maxIT}{
	  integer. The max allowed number of iterations in R built-in function optim. Default value is maxIT = 100.
	 \code{maxIT}. 
}
  \item{b}{
	  numeric. 
	 A vector of concentration parameters used in Dirichlet distribution. Default value is b = c(2,2,2). 
}
  \item{converge_threshold}{
	  numeric.
	  One of the two termination criteria of iteration. 
The smaller this value is set, the harder the optimization procedure in eLNNpaired will be considered to be converged. Default value is converge_threshold \code{= 1e-6}.
}
  \item{optimMethod}{
character. Indicates the method for optimization. {\code{\link{optim}}}.
  }
  \item{bound.alpha}{
  numeric. A vector of 2 positive numbers to specify lower and upper bound of estimate of \eqn{\alpha_c}, \code{c="OE", "UE", or "NE"}.
  }
    \item{bound.beta}{
  numeric. A vector of 2 positive numbers to specify lower and upper bound of estimate of \eqn{\beta_c}, \code{c="OE", "UE", or "NE"}.
  }
    \item{bound.k}{
  numeric. A vector of 2 positive numbers to specify lower and upper bound of estimate of \eqn{k_c}, \code{c="OE", "UE", or "NE"}.
  }
    \item{bound.eta}{
  numeric. A vector of \code{p+1} positive numbers to specify lower and upper bound of estimate of \eqn{\eta_c}, \code{c="OE", "UE", or "NE"},
where \code{p} is the number of covariates.
  }

  \item{mc.cores}{
  integer. A positive integer specifying number of computer cores to be used by parallel computing.
  }
  \item{temp0}{
  numeric. Initial temperature in simulated-annealing modified EM.
  }
  \item{r_cool}{
  numeric. Cooling rate in simulated-annealing modified EM, which is inside interval \eqn{(0, 1)}.
  }
  \item{verbose}{
	  logic.
	  An indicator variable telling if print out intermediate results: FALSE for not printing out, TRUE for printing out. Default value is verbose = False.
}
}
\details{
A gene will be assigned to cluster \dQuote{NE} if its posterior probability for non-differentially expressed gene cluster is the largest. 
A gene will be assigned to cluster \dQuote{OE} if its posterior probability for over-expressed gene cluster is the largest. A gene will be assigned to cluster \dQuote{UE} if its responsibility for under-expressed gene cluster is the largest.
}
\value{
	A list of 9 elementes: 
  \item{par.ini}{initial estimate of parameter}
  \item{par.final}{A vector of the estimated model parameters in original scale.}
  \item{memGenes}{probe cluster membership based on eLNNpairedCovSEM algorithm.}
  \item{memGenes2}{probe cluster membership based on eLNNpairedCovSEM algorithm. 2-categories: "DE" indicates differentially expressed; "NE" indicates non-differentially expressed.}
  \item{memGenes.limma}{probe cluster membership based on limma.}
  \item{res.ini}{results of limma analysis}
  \item{update_info}{object returned by \code{\link{optim}} function}
  \item{wmat}{matrix of responsibilities}
  \item{iter.EM}{number of EM iterations.}
  \item{tempFinal}{final temperature in simulated-annealing modification EM}
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

\examples{
data(esDiff)

res.SEM = eLNNpairedCovSEM(EsetDiff = esDiff, 
		    fmla = ~Age + Sex, 
		    probeID.var = "probeid", 
		    gene.var = "gene", 
		    chr.var = "chr",
		    scaleFlag = FALSE,
		    mc.cores = 1,
        verbose = TRUE)

# true probe cluster membership
memGenes.true = fData(esDiff)$memGenes.true
print(table(memGenes.true))

# probe cluster membership
memGenes.limma = res.SEM$memGenes.limma
print(table(memGenes.limma))

# final probe cluster membership
memGenes.SEM = res.SEM$memGenes
print(table(memGenes.SEM))

# cross tables
print(table(memGenes.true, memGenes.limma))
print(table(memGenes.true, memGenes.SEM))

# accuracies
print(mean(memGenes.true == memGenes.limma))
print(mean(memGenes.true == memGenes.SEM))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }% use one of  RShowDoc("KEYWORDS")
