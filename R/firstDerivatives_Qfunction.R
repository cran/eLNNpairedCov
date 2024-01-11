###created by Weiliang, 2021/09/25
# obtain first derivatives of Q function
#


##################################################
# first partial derivatives of logf: d logf / d psi
##################################################

dlogf123 <- function(
  psi, # vector of model parameters
  dat, # G x n expression data matrix
  tX, # n x (p+1) design matrix
  dg2Vec, # G x 1 vector of d_g^T*d_g
  mc.cores = 4
)
{
  G = nrow(dat) # number of genes
  p = ncol(tX)-1 # number of covariates
  n = nrow(tX) # number of subjects
  
  nlog2pi.2 = n * log(2*pi)/2
  nd2 = n/2
  X = t(tX) # (p+1) x n matrix
  
  ######
  dlogf1.deta1.n = paste("dlogf1.deta1", c(0, 1:p), sep=".")
  dlogf2.deta2.n = paste("dlogf2.deta2", c(0, 1:p), sep=".")
  dlogf3.deta3.n = paste("dlogf3.deta3", c(1:p), sep=".")
  
  ####
  # cluster 1
  ####
  ## dg~N(mu_g, tau_g^{-1} I_n), mu_g~N(exp(X^T*eta1), k1 tau_g^{-1} I_n)
  # tau_g ~ Gamma(alpha_1, beta_1)
  ###parameters in cluster 1
  alpha1 = psi[1]
  beta1 = psi[2]
  k1 = psi[3]
  eta1 = psi[(3+1):(3+p+1)] # (p+1) x 1 vector
  
  # xi1 = exp(X^T*eta1) - n x 1 vector
  # just in case of missing values in tX
  xi1 = apply(tX, 1, function(xg) {
    res = exp(sum(xg*eta1, na.rm = TRUE))
    return(res)
  })

  # d_g^T * xi_i: G x 1 vector
  dgTxi1 = dat%*% xi1
  
  # xi1^T * xi1 - 1x1 scalar
  xi1.2 = sum(xi1^2, na.rm = TRUE)

  # M1: Gx1 vector
  k1.1 = k1+1
  k1.2 = 2*k1.1
  M1 = -dg2Vec/k1.2 + dgTxi1/k1.1 - xi1.2/k1.2
  beta1mM1 = beta1 - M1 # G x 1

  dlogf1.dalpha1 = log(beta1) + digamma(nd2 + alpha1) - log(beta1mM1) - digamma(alpha1)
  dlogf1.dbeta1 = alpha1/beta1 - (nd2+alpha1)/(beta1mM1)

  # 
  k1sq=(k1+1)^2
  k1sq.2=2*k1sq
  
  dM1.dk1 = dg2Vec/k1sq.2 - dgTxi1/k1sq + xi1.2/k1sq.2 # G x 1
  dlogf1.dk1 = - nd2/k1.1+(nd2+alpha1)*dM1.dk1/beta1mM1 # G x 1
  
  #
  xi1Mat = matrix(rep(xi1, G), ncol = n, byrow = TRUE) # G x n matrix
  diag.dg.xi1 = dat*xi1Mat # G x n matrix
  mat1 = diag.dg.xi1 %*% tX # G x (p+1)
  vec1 = X %*% diag(xi1) %*% xi1 # (p+1) x 1
  vec1Mat = matrix(rep(vec1, G), ncol=p+1, byrow=TRUE) # G x (p+1)
  dM1.deta1 = (mat1-vec1Mat)/k1.1 # G x (p+1)
  
  beta1mM1mat = matrix(rep(beta1mM1, p+1), ncol = p +1, byrow = FALSE) # G x (p+1)
  dlogf1.deta1 = (nd2+alpha1) * dM1.deta1 / (beta1mM1mat) # G x (p+1)

  ####
  # cluster 2
  ###
  ### dg~N(mu_g, tau_g^{-1} I_n), mu_g~N(-exp(X^T*eta2), k2 tau_g^{-1} I_n)
  # tau_g ~ Gamma(alpha_2, beta_2)
  ##parameters in cluster 2
  alpha2=psi[(3+p+1)+1]
  beta2=psi[(3+p+1)+2]
  k2=psi[(3+p+1)+3]
  eta2=psi[((p+7)+1):((p+7)+(p+1))] # (p+1) x 1 vector
  
  # xi2 = exp(X^T*eta2) - n x 1 vector
  xi2 = apply(tX, 1, function(xg) {
    res = exp(sum(xg*eta2, na.rm = TRUE))
    return(res)
  })
  
  # d_g^T * xi_i: G x 1 vector
  dgTxi2 = dat%*% xi2
  
  # xi2^T * xi2 - 1 x 1 scalar
  xi2.2 = sum(xi2^2, na.rm = TRUE)
  
  # M2: Gx1 vector
  k2.1 = k2+1
  k2.2 = 2*k2.1
  M2 = -dg2Vec/k2.2 - dgTxi2/k2.1 - xi2.2/k2.2
  beta2mM2 = beta2 - M2 # G x 1
  
  dlogf2.dalpha2 = log(beta2) + digamma(nd2 + alpha2) - log(beta2mM2) - digamma(alpha2)
  dlogf2.dbeta2 = alpha2/beta2 - (nd2+alpha2)/(beta2mM2)
  
  # 
  k2sq=(k2+1)^2
  k2sq.2=2*k2sq
  
  dM2.dk2 = dg2Vec/k2sq.2 + dgTxi2/k2sq + xi2.2/k2sq.2 # G x 1
  dlogf2.dk2 = - nd2/k2.1+(nd2+alpha2)*dM2.dk2/beta2mM2 # G x 1
  
  #
  xi2Mat = matrix(rep(xi2, G), ncol = n, byrow = TRUE) # G x n matrix
  diag.dg.xi2 = dat*xi2Mat # G x n matrix
  mat2 = diag.dg.xi2 %*% tX # G x (p+1)
  vec2 = X %*% diag(xi2) %*% xi2 # (p+1) x 1
  vec2Mat = matrix(rep(vec2, G), ncol=p+1, byrow=TRUE) # G x (p+1)
  dM2.deta2 = (-mat2-vec2Mat)/k2.1 # G x (p+1)
  
  beta2mM2mat = matrix(rep(beta2mM2, p+1), ncol = p +1, byrow = FALSE) # G x (p+1)
  dlogf2.deta2 = (nd2+alpha2) * dM2.deta2 / (beta2mM2mat) # G x (p+1)
  
    

  ####
  # cluster 3
  ####
  ### d_g~N(U^T*theta_g, tau_g^{-1} I_n), theta_g~N(eta3, k_3 tau_g^{-1} I_p)
  #   tau_g ~ Gamma(alpha3, beta3)
  alpha3 = psi[(2*p+8)+1] 
  beta3 = psi[(2*p+8)+2]
  k3 = psi[(2*p+8)+3]
  eta3=psi[((2*p+11)+1):((2*p+11)+p)] # p x 1 vector
  
  ###density for genes in cluster 3 (non-diff-expressed) 
  # remove the 1st column of tX
  tU = tX[, -1,drop=FALSE] # n x p
  U = t(tU) # p x n
  Ipk3UtU = diag(p)+k3*U%*%tU # p x p
  iIpk3UtU = MASS::ginv(Ipk3UtU) # p x p
  det.Ipk3UtU = det(Ipk3UtU)
  eta3.2 = sum(eta3^2, na.rm = TRUE) # 1x1 scalar
  
  # U*d_g
  Udg = dat%*%tU #Gxp matrix
  eta3Mat = matrix(rep(eta3, G), ncol = p, byrow = TRUE) # Gxp matrix
  #  k3*U*d_g + eta3
  Udgpeta3 = k3*Udg + eta3Mat #Gxp matrix
  
  # Cg3: Gx1 vector
  ttLst = mclapply(X = 1:G, FUN = function(g) {
    Udgpeta3.g = c(Udgpeta3[g,]) # p x 1
    Cg3 = c(t(Udgpeta3.g) %*% iIpk3UtU %*% Udgpeta3.g)
    
    Ud.g = c(Udg[g,]) # p x 1
    Cg4 = c(t(Ud.g) %*% iIpk3UtU %*% Udgpeta3.g)
    
    Cg5 = c(t(Udgpeta3.g) %*% iIpk3UtU %*% U %*% tU %*% iIpk3UtU %*% Udgpeta3.g)
    
    res = c(Cg3, Cg4, Cg5)
    names(res) = c("Cg3", "Cg4", "Cg5")
    return(res)
  }, mc.cores = mc.cores)
  
  Cg345Mat = t(sapply(ttLst, function(x){x})) # G x 3
  Cg3 = c(Cg345Mat[,1]) # (k3 U dg + eta3)^T (Ip + k3 U U^T)^{-1} (k3 U dg + eta3)
  Cg4 = c(Cg345Mat[,2]) # (U dg)^T (Ip + k3 U U^T)^{-1} (k3 U dg + eta3)
  Cg5 = c(Cg345Mat[,3]) # (k3 U dg + eta3)^T (Ip + k3 U U^T)^{-1} U U^T (Ip + k3 U U^T)^{-1} (k3 U dg + eta3)
  
  # M3: Gx1 vector
  k3.2 = 2*k3
  M3 = dg2Vec/2 + eta3.2/k3.2 - Cg3/k3.2 # G x 1
  beta3pM3 = beta3 + M3 # G x 1  
  
  dlogf3.dalpha3 = log(beta3) + digamma(nd2 + alpha3) - log(beta3pM3) - digamma(alpha3)
  dlogf3.dbeta3 = alpha3/beta3 - (nd2+alpha3)/(beta3pM3)
  
  # dM3.deta3
  dM3.deta3 = eta3Mat/k3 - Udgpeta3 %*% iIpk3UtU / k3 # G x p
  beta3pM3Mat = matrix(rep(beta3pM3, p), ncol = p, byrow = FALSE) # G x p
  dlogf3.deta3 = (nd2 + alpha3) * dM3.deta3 / beta3pM3Mat # G x p
  
  # dM3.dk3: G x 1
  k3sq = k3^2
  k3sq.2 = 2*k3sq
  
  part1 = -eta3.2/k3sq.2 + Cg3/k3sq.2 # G x 1
  part2 = - Cg4 / k3 # G x 1
  part3 = Cg5 / k3sq.2 # G x 1
  
  dM3.dk3 = part1 + part2 + part3 # G x 1
    
  dlogf3.dk3 = (nd2 + alpha3) * dM3.dk3 / beta3pM3 # G x 1

  ##################
  
  dlogfMat = cbind(dlogf1.dalpha1, dlogf1.dbeta1, dlogf1.dk1, dlogf1.deta1,
                   dlogf2.dalpha2, dlogf2.dbeta2, dlogf2.dk2, dlogf2.deta2,
                   dlogf3.dalpha3, dlogf3.dbeta3, dlogf3.dk3, dlogf3.deta3
                   )
  colnames(dlogfMat) = c(
    "dlogf1.dalpha1",
    "dlogf1.dbeta1",
    "dlogf1.dk1",
    dlogf1.deta1.n,
    "dlogf2.dalpha2",
    "dlogf2.dbeta2",
    "dlogf2.dk2",
    dlogf2.deta2.n,
    "dlogf3.dalpha3",
    "dlogf3.dbeta3",
    "dlogf3.dk3",
    dlogf3.deta3.n
  )
  
  return (dlogfMat)
}


##############################
# first partial derivatives of Q: d Q / d psi.r
############################
dQfun <- function(
  # QWL: psi should be the first parameter, which is required by 'optim'
  psi, # parameters to be optimized
  t_pi, # cluster mixture proportions (pi.OE, pi.UE, pi.NE) at t-th iteration of EM algorithm
  tilde_z, # E(z | d, pi^{(t)}, Psi^{(t)}) t-th iteration of EM algorithm
  dat, # G x n expression data matrix
  tX, # n x (p+1) design matrix
  mc.cores = 4,
  dg2Vec, # G x 1 vector of d_g^T*d_g
  const.b # const.b = lgamma(b[1]+b[2]+b[3]) - lgamma(b[1]) - lgamma(b[2]) - lgamma(b[3]) + sum((b-1) * log(t_pi), na.rm = TRUE)
)
{
  
  p = ncol(tX) - 1

  # d logf / d Psi
  dlogf=dlogf123(psi, dat, tX, dg2Vec, mc.cores)
  
  #######
  alpha1 = psi[1]
  beta1 = psi[2]
  k1 = psi[3] # k1 = Phi(delta1)
  #delta1 = psi.r[3] # d k1 / d delta1 = phi(delta1)
  eta1 = psi[(3+1):(3+p+1)] # (p+1) x 1 vector
  
  alpha2=psi[(3+p+1)+1]
  beta2=psi[(3+p+1)+2]
  k2=psi[(3+p+1)+3] # k2 = Phi(delta2)
  #delta2 = psi.r[(3+p+1)+3] # d k2 / d delta2 = phi(delta2)
  eta2=psi[((p+7)+1):((p+7)+(p+1))] # (p+1) x 1 vector
  
  alpha3 = psi[(2*p+8)+1] 
  beta3 = psi[(2*p+8)+2]
  k3 = psi[(2*p+8)+3] # k3 = Phi(delta3)
  #delta3 = psi.r[(2*p+8)+3] # d k3 / d delta3 = phi(delta3)
  eta3=psi[((2*p+11)+1):((2*p+11)+p)] # p x 1 vector
  
  ######
  tildez1 = tilde_z[,1]
  dQ.dalpha1 = sum(tildez1 * dlogf[,1], na.rm = TRUE)
  dQ.dbeta1 = sum(tildez1 * dlogf[,2], na.rm = TRUE)
  dQ.dk1 = sum(tildez1 * dlogf[,3], na.rm = TRUE)
  dQ.deta1 = apply(dlogf[, c((3+1):(3+p+1)), drop = FALSE], 2, function(dlogf.i) {
    return( sum(tildez1 * dlogf.i, na.rm = TRUE) )
  })
  
  ######
  tildez2 = tilde_z[,2]
  dQ.dalpha2 = sum(tildez2 * dlogf[,(3+p+1)+1], na.rm = TRUE)
  dQ.dbeta2 = sum(tildez2 * dlogf[,(3+p+1)+2], na.rm = TRUE)
  dQ.dk2 = sum(tildez2 * dlogf[,(3+p+1)+3], na.rm = TRUE)
  dQ.deta2 = apply(dlogf[, ((p+7)+1):((p+7)+(p+1)), drop = FALSE], 2, function(dlogf.i) {
    return( sum(tildez2 * dlogf.i, na.rm = TRUE) )
  })
  
  ######
  tildez3 = tilde_z[,3]
  dQ.dalpha3 = sum(tildez3 * dlogf[,(2*p+8)+1], na.rm = TRUE)
  dQ.dbeta3 = sum(tildez3 * dlogf[,(2*p+8)+2], na.rm = TRUE)
  dQ.dk3 = sum(tildez3 * dlogf[,(2*p+8)+3], na.rm = TRUE)
  dQ.deta3 = apply(dlogf[, ((2*p+11)+1):((2*p+11)+p), drop = FALSE], 2, function(dlogf.i) {
    return( sum(tildez3 * dlogf.i, na.rm = TRUE) )
  })
  
  dQ.deta1.n = paste("dQ.deta1", c(0, 1:p), sep=".")
  dQ.deta2.n = paste("dQ.deta2", c(0, 1:p), sep=".")
  dQ.deta3.n = paste("dQ.deta3", c(1:p), sep=".")
  #####
  res = c(dQ.dalpha1, dQ.dbeta1, dQ.dk1, dQ.deta1,
          dQ.dalpha2, dQ.dbeta2, dQ.dk2, dQ.deta2,
          dQ.dalpha3, dQ.dbeta3, dQ.dk3, dQ.deta3)
  names(res) = c("dQ.dalpha1", "dQ.dbeta1", "dQ.dk1", dQ.deta1.n,
                 "dQ.dalpha2", "dQ.dbeta2", "dQ.dk2", dQ.deta2.n,
                 "dQ.dalpha3", "dQ.dbeta3", "dQ.dk3", dQ.deta3.n)
  return ( res )
  
}
