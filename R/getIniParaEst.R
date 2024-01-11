# created on Sept. 25, 2021
#  (1) obtain initial parameter estimates

getIniParaEst = function(
  EsetDiff,
  fmla = ~Age + Sex,
  probeID.var = "probeid",
  gene.var = "gene",
  chr.var = "chr",
  bound.alpha = c(0.001, 6),
  bound.beta = c(0.001, 6),
  bound.k = c(0.001, 0.9999),
  scaleFlag = TRUE
){ 
  # get expression levels
  # rows of 'dat' are probes; columns are samples
  dat=exprs(EsetDiff)
  if(scaleFlag)
  {
    dat.s=t(apply(dat, 1, function(x) { return(x/sd(x, na.rm=TRUE))}))
    rownames(dat.s)=rownames(dat)
    colnames(dat.s)=colnames(dat)
    exprs(EsetDiff)=dat.s
    dat=exprs(EsetDiff)
  }
  
  # number of probes
  G = nrow(EsetDiff)
  
  # number of samples
  n = ncol(EsetDiff)
  
  # get phenotype data
  pDat=pData(EsetDiff)
  
  # get design matrix
  # the first column of 'tX' is a vector of ones for intercept
  tX=model.matrix(fmla, data=pDat) #n*(p+1)
  
  # QWL: obtain covariate names
  covNames=colnames(tX)
  
  # QWL: calculate p = number of covariates
  p = ncol(tX)-1
  
  ###Get initial estimates from limma
  result_limma = lmFitPaired2(
    esDiff = EsetDiff,
    formula = fmla, 
    pos.var.interest = 0, # focus on intercept
    pvalAdjMethod = "fdr", # control fdr < 0.05
    alpha = 0.05, 
    probeID.var=probeID.var, 
    gene.var = gene.var, 
    chr.var = chr.var, 
    verbose = FALSE)
  
  ####calculate the initial values of para's for each cluster/get initial cluster with limma
  frame_unsorted = result_limma$frame.unsorted
  
  # first check adjust p-values
  over_expressed_sub_script = which(frame_unsorted$stats > 0 & frame_unsorted$p.adj < 0.05)
  under_expressed_sub_script = which(frame_unsorted$stats < 0 & frame_unsorted$p.adj < 0.05)

  cutoff.nSig = 5
  if(length(over_expressed_sub_script) <= cutoff.nSig || length(under_expressed_sub_script) <= cutoff.nSig)
  {
    ## use un-adjusted p-values
    #over_expressed_sub_script = which(frame_unsorted$stats > 0 & frame_unsorted$pval < 0.05)
    #under_expressed_sub_script = which(frame_unsorted$stats < 0 & frame_unsorted$pval < 0.05)


    #if(length(over_expressed_sub_script) <= cutoff.nSig || length(under_expressed_sub_script) <= cutoff.nSig)
    #{
      # use top 100 genes
      frame = result_limma$frame
      frame100 = frame[1:100,]
      over_expressed_sub_script=frame100$pos[which(frame100$stats>0) ]
      under_expressed_sub_script=frame100$pos[which(frame100$stats<0) ]
    #}
  }
  
  myset = 1:G
  
  non_diff_sub_script = myset[-c(over_expressed_sub_script, under_expressed_sub_script)]
  
  # initial memGenes record probe cluster membership
  memGenes.ini = rep("OE",G)
  memGenes.ini[under_expressed_sub_script]="UE"
  memGenes.ini[non_diff_sub_script]="NE"
  
  # pi.overExpressed, pi.UnderExpressed, pi.non-differentially expressed
  pi_prior = c(
    length(over_expressed_sub_script)/G,
    length(under_expressed_sub_script)/G,
    0)
  pi_prior[3] = 1 - pi_prior[1] - pi_prior[2]
  # add names for pi_prior
  names(pi_prior) = c("pi.OE", "pi.UE", "pi.NE")
  
  over_expressed_EsetDiff = EsetDiff[over_expressed_sub_script,]
  under_expressed_EsetDiff = EsetDiff[under_expressed_sub_script,]
  non_diff_EsetDiff = EsetDiff[non_diff_sub_script,]
  
  ####################
  ##cluster 1###
  ###parameters in cluster 1: alpha1, beta1, k1, eta1
  data_matrix_of_EsetDiff = exprs(over_expressed_EsetDiff) # G1 x n matrix
  
  # QWL: for each subject, get median expression levels as an estimate average of hat{mu}_g
  muhat = apply(data_matrix_of_EsetDiff, 2, median, na.rm= TRUE) # n x 1 vector
  
  # QWL: hat{tau}_g^{(0)}, g=1,..., G1
  temp_tau = 1 / (apply(data_matrix_of_EsetDiff, 1, mad, na.rm=TRUE)^2) # G1 x 1 vector
  
  # QWL: add description
  # omega = k1*average(tau_g^{-1})
  # so k1 = omega*average(tau_g)
  omega = mad(muhat)^2
  k1 = omega * median(temp_tau, na.rm=TRUE)
  
  tt.median = median(temp_tau, na.rm=TRUE)
  tt.mad = mad(temp_tau, na.rm=TRUE)
  alpha1 = (tt.median/tt.mad)^2
  beta1 = tt.median/(tt.mad^2)
  
  if(alpha1 > bound.alpha[2] || alpha1 < bound.alpha[1])
  {
    alpha1 = runif(1, min = bound.alpha[1], max = bound.alpha[2])
  }
  if(beta1 > bound.beta[2] || beta1 < bound.beta[1])
  {
    beta1 = runif(1, min = bound.beta[1], max = bound.beta[2])
  }
  
  # QWL: estimate eta_1
  # to handle missing values, we can use lm()
  pDat$logmu1=log(abs(muhat)+0.001)
  fmla1 = as.formula(paste("logmu1", paste(as.character(fmla), collapse=""), collapse = ""))
  tt1=lm(formula=fmla1, data=pDat)
  eta1=tt1$coefficients
  
  ## End of cluster1##
  
  
  ####################
  ##cluster 2 ##
  ###parameters in cluster 2: alpha2, beta2, k2, eta2
  data_matrix_of_EsetDiff = exprs(under_expressed_EsetDiff)
  
  # QWL: for each subject, get median expression levels as an estimate average of hat{mu}_g
  muhat = apply(data_matrix_of_EsetDiff, 2, median, na.rm= TRUE) # n x 1 vector
  
  # estimate of tau_{g}
  temp_tau = 1 / (apply(data_matrix_of_EsetDiff, 1, mad, na.rm=TRUE)^2)
  
  # QWL: add description
  # omega = k2*average(tau_g^{-1})
  # so k2 = omega*average(tau_g)
  omega = mad(muhat)^2
  k2 = omega * median(temp_tau, na.rm=TRUE)
  
  tt.median = median(temp_tau, na.rm=TRUE)
  tt.mad = mad(temp_tau, na.rm=TRUE)
  alpha2 = (tt.median/tt.mad)^2
  beta2 = tt.median/(tt.mad^2)
  
  if(alpha2 > bound.alpha[2] || alpha2 < bound.alpha[1])
  {
    alpha2 = runif(1, min = bound.alpha[1], max = bound.alpha[2])
  }
  if(beta2 > bound.beta[2] || beta2 < bound.beta[1])
  {
    beta2 = runif(1, min = bound.beta[1], max = bound.beta[2])
  }
 
  # get initial eat from fomula
  
  # QWL: estimate eta_2
  # to handle missing values, we can use lm()
  pDat$lognmu2=log(abs(-muhat)+0.001)
  fmla2 = as.formula(paste("lognmu2", paste(as.character(fmla), collapse=""), collapse = ""))
  tt2=lm(formula=fmla2, data=pDat)
  eta2=tt2$coefficients
  
  ##End of cluster 2##
  
  ####################
  ####################
  ## cluster 3####
  ######parameters in cluster 3: alpha3, beta3, k3, eta3
  data_matrix_of_EsetDiff = exprs(non_diff_EsetDiff)
  
  # estimate of tau_{g}
  temp_tau = 1 / (apply(data_matrix_of_EsetDiff, 1, mad, na.rm=TRUE)^2)
  
  tt.median = median(temp_tau, na.rm=TRUE)
  tt.mad = mad(temp_tau, na.rm=TRUE)
  alpha3 = (tt.median/tt.mad)^2
  beta3 = tt.median/(tt.mad^2)

  if(alpha3 > bound.alpha[2] || alpha3 < bound.alpha[1])
  {
    alpha3 = runif(1, min = bound.alpha[1], max = bound.alpha[2])
  }
  if(beta3 > bound.beta[2] || beta3 < bound.beta[1])
  {
    beta3 = runif(1, min = bound.beta[1], max = bound.beta[2])
  }
 
  # QWL: we can use limma results to roughly estimate theta_g 
  ebFit=result_limma$ebFit
  coef=ebFit$coefficients
  coef3=coef[non_diff_sub_script,,drop=FALSE]
  theta_g = coef3[, -1, drop=FALSE]

  eta3 = apply(theta_g, 2, median, na.rm = TRUE)
  omega=sum(diag(var(theta_g, use = "na.or.complete")), na.rm=TRUE)
  k3 = omega * median(temp_tau, na.rm=TRUE)
  
  ######end of cluster 3 ###
  
  if(k1 > bound.k[2] || k1 < bound.k[1])
  {
    k1 = runif(1, min = bound.k[1], max = bound.k[2])
  }
  if(k2 > bound.k[2] || k2 < bound.k[1])
  {
    k2 = runif(1, min = bound.k[1], max = bound.k[2])
  }
  if(k3 > bound.k[2] || k3 < bound.k[1])
  {
    k3 = runif(1, min = bound.k[1], max = bound.k[2])
  }
  
  #psi = c(alpha1, beta1, k1, eta1, alpha2, beta2, k2, eta2, alpha3, beta3, k3, eta3)
  # get more stable initial parameter estimates
  psi = c((alpha1+alpha2)/2, (beta1+beta2)/2, (k1+k2)/2, (eta1+eta2)/2, 
          (alpha1+alpha2)/2, (beta1+beta2)/2, (k1+k2)/2, (eta1+eta2)/2, 
          alpha3, beta3, k3, eta3)
  
  # QWL: need to have adaptive length for 'eta1', 'eta2' and 'eta3'
  names(psi)=c("alpha1", "beta1", "k1", 
               paste("eta1", c(0, 1:p), sep="."), 
               "alpha2", "beta2", "k2", 
               paste("eta2", c(0, 1:p), sep="."), 
               "alpha3", "beta3", "k3", 
               paste("eta3", 1:p, sep=".")) 
  
  
  t_pi = pi_prior
  
  memGenes.limma = rep("NE", G)
  memGenes.limma[which(result_limma$memGenes == 1)] = "OE"
  memGenes.limma[which(result_limma$memGenes == 3)] = "UE"
  

  res = list(psi = psi, 
             t_pi = t_pi, 
             memGenes.ini = memGenes.ini, 
             memGenes.limma = memGenes.limma,
             tX = tX,
             res.limma = result_limma)
  invisible(res)
}
