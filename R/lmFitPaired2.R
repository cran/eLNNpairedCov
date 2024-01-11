# copied and modified from iCheck package

lmFitPaired2=function (esDiff, formula = ~1, pos.var.interest = 0, 
                       pvalAdjMethod = "fdr", 
          alpha = 0.05, probeID.var = "ProbeID", gene.var = "Symbol", 
          chr.var = "Chromosome", verbose = TRUE) 
{
  diffDat <- exprs(esDiff)
  pDat <- pData(esDiff)
  fDat <- fData(esDiff)
  probeIDs = fDat[, c(probeID.var)]
  geneSymbols = fDat[, c(gene.var)]
  chr = fDat[, c(chr.var)]
  
  pDat2 <- pDat
  rownames(pDat2) <- 1:nrow(pDat2)
  if (is.null(formula)) {
    designMat <- matrix(rep(1, ncol(diffDat)), ncol = 1)
    rownames(designMat) <- 1:nrow(designMat)
    colnames(designMat) <- "(Intercept)"
  }
  else {
    flist = as.list(formula)
    if (length(flist) > 1) {
      if (flist[[2]] != 1) {
        designMat <- model.matrix(formula, data = pDat2)
      }
      else {
        designMat <- matrix(rep(1, ncol(diffDat)), ncol = 1)
        rownames(designMat) <- 1:nrow(designMat)
        colnames(designMat) <- "(Intercept)"
      }
    }
    else {
      
      stop("Examples of formula are ~1 or ~cov1+cov2\nformula not correct.\n")
    }
  }
  rn <- as.numeric(rownames(designMat))
  diffDat <- diffDat[, rn, drop = FALSE]
  rownames(designMat) <- colnames(diffDat)
  if (verbose) {
    cat("Running lmFit for paired data ...\n")
  }
  fit = lmFit(diffDat, designMat)
  if (verbose) {
    cat("Running eBayes...\n")
  }
  ebFit = eBayes(fit)
  if (verbose) {
    cat("Preparing output...\n")
  }
  pval<- ebFit$p.value[,1]
  stats = ebFit$t[,1]
  p.adj = p.adjust(pval, method=pvalAdjMethod)
 
  frame.unsorted=data.frame(probeIDs = probeIDs,
    gene=geneSymbols, chr=chr, stats=stats, pval=pval,
    p.adj = p.adj, pos=1:length(pval)) 

  frame = frame.unsorted[order(frame.unsorted$pval),]

  memGenes=rep(2, nrow(esDiff))
  memGenes[which(stats>0 & p.adj<alpha)] = 1
  memGenes[which(stats<0 & p.adj<alpha)] = 3

  memGenes2=rep(1, nrow(esDiff))
  memGenes2[which(memGenes==2)]=0

  res = list(memGenes=memGenes, memGenes2=memGenes2,
    frame.unsorted=frame.unsorted, frame=frame,
    ebFit=ebFit)

  invisible(res)
}

