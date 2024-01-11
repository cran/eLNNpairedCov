###modified by Weiliang, 2021/09/27
#  simplify functions

###created by Weiliang, 2019/10/19
# (1) use EM algorithm to udpate hyperparameters

###created by Weiliang, 2019/10/15
#  (1) different clusters have different hyperparameters
#
###created by Weiliang, 2019/10/14
#  (1) add scaling option 'scaleFlag'. By default, will do probe-wise
#       scaling, but not centering
#  (2) set alpha1=alpha3, and alpha2=alpha3; beta1=beta3, beta2=beta3
#  (3) set k1=k2=k3=(k1+k2)/2
#  (4) add a convergence criterion: negative Q func (new) > negative Q func (old)

###created by Weiliang, 2019/10/13
#  (1) only update mixture proportion

###modifed by Weiliang, 2019/10/12
#  (1) do not use EM algorithm to estimate model parameters,
#      except for cluster mixture proportions
#  (2) for clusters 1 (OE) and 2 (UE), we use the same hyperparameters
#
###modifed by Weiliang, 2019/04/15
###Created by Yixin, 2018/05/14


eLNNpairedCov = function(
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
  verbose=FALSE
){ 
  # get expression levels
  # rows of 'dat' are probes; columns are samples
  dat=exprs(EsetDiff) # G x n matrix
  
  if(scaleFlag) # recommended
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

  # get initial parameter estimates
  res.ini = getIniParaEst(
    EsetDiff = EsetDiff,
    fmla = fmla,
    probeID.var = probeID.var,
    gene.var = gene.var,
    chr.var = chr.var,
    bound.alpha = bound.alpha,
    bound.beta = bound.beta,
    bound.k = bound.k,
    scaleFlag = FALSE
  )

  # transpose of design matrix
  tX = res.ini$tX # n x (p+1) matrix
  p = ncol(tX) - 1 # number of covariates

  # calculate d_g^T*d_g (G x 1 vector)
  dg2Vec = apply(dat, 1, function(dg)
  {
      dg2 = sum(dg^2, na.rm = TRUE)
      return(dg2)
  }) 

  # initial estimates of cluster mixture proportions
  t_pi = res.ini$t_pi
  nm.pi = c("OE", "UE", "NE")
  names(t_pi) = nm.pi
  
  # initial parameter estimates
  psi = res.ini$psi
  nm.psi=c("alpha1", "beta1", "k1", 
               paste("eta1", c(0, 1:p), sep="."), 
               "alpha2", "beta2", "k2", 
               paste("eta2", c(0, 1:p), sep="."), 
               "alpha3", "beta3", "k3", 
               paste("eta3", 1:p, sep=".")) 
  names(psi) = nm.psi
  
  #####
  ##EM algorithm#
  #####
  # use EM algorithm only to estimate mixture proportions
  par.ini=c(t_pi, psi)
  par.old = par.ini

  if(verbose)
  {
    cat("initial parameter estimates>>>\n")
    print(par.ini)
    cat("\n")
  }

  iter = 0
  const.b1 = lgamma(b[1]+b[2]+b[3]) - lgamma(b[1]) - lgamma(b[2]) - lgamma(b[3])
  
  n.psi = length(psi)

  lower.psi = rep(-bound.eta[1], n.psi)
  lower.psi[1] = bound.alpha[1] # for alpha1
  lower.psi[2] = bound.beta[1] # for beta1
  lower.psi[3] = bound.k[1] # for k1
  lower.psi[(3+p+1)+1] = bound.alpha[1] # for alpha2
  lower.psi[(3+p+1)+2] = bound.beta[1] # for beta2
  lower.psi[(3+p+1)+3] = bound.k[1] # for k2
  lower.psi[(2*p+8)+1] = bound.alpha[1] # for alpha3
  lower.psi[(2*p+8)+2] = bound.beta[1] # for beta3
  lower.psi[(2*p+8)+3] = bound.k[1] # for k3
  
  upper.psi = rep(bound.eta[2], n.psi)
  upper.psi[1] = bound.alpha[2] # for alpha1
  upper.psi[2] = bound.beta[2] # for beta1
  upper.psi[3] = bound.k[2] # for k1
  upper.psi[(3+p+1)+1] = bound.alpha[2] # for alpha2
  upper.psi[(3+p+1)+2] = bound.beta[2] # for beta2
  upper.psi[(3+p+1)+3] = bound.k[2] # for k2
  upper.psi[(2*p+8)+1] = bound.alpha[2] # for alpha3
  upper.psi[(2*p+8)+2] = bound.beta[2] # for beta3
  upper.psi[(2*p+8)+3] = bound.k[2] # for k3
  

  updatePiOnly = FALSE

  while (iter<Maxiter)
  {
    iter = iter + 1
    
    # QWL: E-step
    tilde_z = Get_tilde_z(
      t_pi = t_pi, # cluster mixture proportions (pi.OE, pi.UE, pi.NE)
      psi = psi, # vector of model parameters
      dat = dat, # G x n expression data matrix
      tX = tX, # n x (p+1) design matrix
      dg2Vec = dg2Vec, # G x 1 vector of d_g^T*d_g
      mc.cores = mc.cores
    )

    if(!updatePiOnly)
    {
      # M-step
      # first update psi
      const.b = const.b1 + sum((b-1) * log(t_pi), na.rm = TRUE)
      
      update_info = optim(par = psi, 
                          fn = Qfun,
                          gr = dQfun,
                          t_pi = t_pi, 
                          tilde_z = tilde_z, 
                          dat = dat, 
                          tX = tX,
                          dg2Vec = dg2Vec,
                          mc.cores = mc.cores,
                          const.b = const.b,
                          method = optimMethod,
                          lower = lower.psi,
                          upper = upper.psi,
                          control = list(maxit = maxIT, fnscale = -1))
  
     # QWL: par.new includes both piVec and psi
      psi.new=update_info$par
      names(psi.new) = nm.psi

  
    } else {
      
      psi.new = psi
    }

    # then update piVec
    t_pi.new = Get_t_pi(
      G = G,
      tilde_z = tilde_z,
      b = b)
    names(t_pi.new) = nm.pi
    
   
    par.new=c(t_pi.new , psi.new)
    names(par.new) = c(nm.pi, nm.psi)
    
    #sumSqDiff=sum((par.old - par.new)^2, na.rm=TRUE)
    # proportional change
    pChange = max(abs((par.old - par.new)/par.old), na.rm = TRUE)
    if(verbose)
    {
      # cat(" iter=", iter)
      print(c("EM iteration:", iter))
      cat('par >>\n',par.new,'\n')
      
      #cat("sum(|diff|^2)=", sumSqDiff, "\n")
      cat("pChange=", pChange, "\n")
    }
    
    
    #if(sumSqDiff < converge_threshold)
    if(pChange < converge_threshold) 
    {
      break
    }
    
    # QWL: update parameters


    flag1 = sum((psi.new[1] %in% bound.alpha) | (psi.new[2] %in% bound.beta) | (psi.new[3] %in% bound.k) | any(psi.new[(3+1):(3+p+1)] %in% bound.eta), na.rm = TRUE)
    flag2 = sum((psi.new[(3+p+1)+1] %in% bound.alpha) | (psi.new[(3+p+1)+2] %in% bound.beta) | (psi.new[(3+p+1)+3] %in% bound.k) | any(psi.new[((p+7)+1):((p+7)+(p+1))] %in% bound.eta), na.rm = TRUE)
    flag3 = sum((psi.new[(2*p+8)+1] %in% bound.alpha) | (psi.new[(2*p+8)+2] %in% bound.beta) | (psi.new[(2*p+8)+3] %in% bound.k) | any(psi.new[((2*p+11)+1):((2*p+11)+p)] %in% bound.eta), na.rm = TRUE)
    flag = flag1 + flag2 + flag3

    if(flag)
    {
      updatePiOnly = TRUE

      par.new=c(t_pi, psi)
      names(par.new) = c(nm.pi, nm.psi)
    } else {
      t_pi=t_pi.new
      psi = psi.new
    }

    par.old = par.new
   
  }

  memGenes = getMemFunc(tilde_z = tilde_z)
  
  # QWL: create 2-category membership: "NE" - non-DE; "DE" - DE
  memGenes2=rep("NE", length(memGenes))
  pos.sig=which(memGenes != "NE")
  if(length(pos.sig))
  {
    memGenes2[pos.sig]="DE"
  }
  
  pi.new = rep(NA, 3)
  pi.new[1] = mean(memGenes == "OE", na.rm = TRUE)
  pi.new[2] = mean(memGenes == "UE", na.rm = TRUE)
  pi.new[3] = 1 - pi.new[1] - pi.new[2]
  par.new[1:3] = pi.new

  result=list(
    par.ini=par.ini,
    par.final = par.new, # in original scale
    memGenes = memGenes,
    memGenes2 = memGenes2,
    memGenes.limma = res.ini$memGenes.limma,
    res.ini = res.ini,
    update_info = update_info,
    wmat = tilde_z,
    iter=iter)
  
  # QWL: add summary
  if(verbose)
  {
    cat("\n***** Summary of the analysis>>>>\n")
    cat("\nNumber of subjects>>", ncol(EsetDiff), "\n")
    cat("\nNumber of probes>>", nrow(EsetDiff), "\n")
    cat("\nNumber of covariates>>", p, "\n")
    cat("\nformula for adjusting covariates>>", as.character(fmla), "\n")
    cat("\nNumber of EM iteration>>", iter, "\n")
    cat("\nFrequencies of the 3 clusters>>\n")
    print(table(result$memGenes, useNA="ifany"))
    cat("\nParameter estimates>>\n")
    print(result$par.final)
    cat("\n************************\n")
  } 
  
  invisible(result)

}

