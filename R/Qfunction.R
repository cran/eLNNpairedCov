###modified by Yixin, 2023/06/26
#add simulated annealing modified EM version to calculate \tilde z
#

###Modified by Weiliang, 2021/10/18
# (1) use matrix operation to speed up
#
###Modified by Weiliang, 2021/09/27
#  simplify functions
#
###Modified by Weiliang, 2019/10/13
# logf123: use alpha1, beta1, k1, instead of alpha, beta, and k
#
###Created by Yixin, 2018/03/28
# likelihood function and Q function
##sum_{i=1}^{nGenes} log10(pi.1*f1(xi)+pi.2*f2(xi)+pi.3*f3(xi))
##dat 1000*300 matrix
##design matrix X 3*300 matrix
  ##here eta3 is 2*1 matrix(which does not contain the intercept)
### To avoid matrix operations to handle missing values, use sum() to replace the inner product of the matrix

# 1
###
# Calculate the log density for each cluster for a given gene g
###
#' @param psi parameters in each cluster
#     psi=c(alpha1,beta1,k_1, beta10, beta11, beta12, 
#       alpha2,beta2,k_2, beta20, beta21, beta22, 
#       alpha3,beta3,k3, beta31, beta32)
# QWL: dat is a G x n matrix
# QWL: tX = X^T is a n x (p+1) design matrix
logf123 <- function(
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

  # xi1^T * xi1 - 1 x 1 scalar
  xi1.2 = sum(xi1^2, na.rm = TRUE)

  # M1: Gx1 vector
  k1.1 = k1+1
  k1.2 = 2*k1.1
  M1 = -dg2Vec/k1.2 + dgTxi1/k1.1 - xi1.2/k1.2

  logf1 = alpha1 * log(beta1) + lgamma(nd2+alpha1) - n*log(k1.1)/2 - (nd2 + alpha1)*log(beta1-M1) 
  logf1 = logf1 - nlog2pi.2 - lgamma(alpha1) 


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

  # d_g^T * xi_2: G x 1 vector
  dgTxi2 = dat%*% xi2
  
  # xi2^T * xi2 - 1 x 1 scalar
  xi2.2 = sum(xi2^2, na.rm = TRUE)

  # M2: Gx1 vector
  k2.1 = k2+1
  k2.2 = 2*k2.1
  M2 = -dg2Vec/k2.2 - dgTxi2/k2.1 - xi2.2/k2.2

  logf2 = alpha2 * log(beta2) + lgamma(nd2+alpha2) - nd2*log(k2.1) -(nd2 + alpha2)*log(beta2 - M2) 
  logf2 = logf2 - nlog2pi.2 - lgamma(alpha2) 
  
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
  tU = tX[,-1,drop=FALSE] # n x p
  U = t(tU) # p x n
  Ipk3UtU = diag(p)+k3*U%*%tU
  iIpk3UtU = MASS::ginv(Ipk3UtU) # p x p 
  det.Ipk3UtU = det(Ipk3UtU)
  eta3.2 = sum(eta3^2, na.rm = TRUE) # 1x1 scalar

  # U*d_g: Gxp matrix
  Udg = dat%*%tU
  eta3Mat = matrix(rep(eta3, G), ncol = p, byrow = TRUE) # Gxp matrix
  #  k3*U*d_g + eta3: Gxp matrix
  Udgpeta3 = k3*Udg + eta3Mat

  # Cg3: Gx1 vector
  Cg3 = unlist(mclapply(X = 1:G, FUN = function(g) {
    Udgpeta3.g = c(Udgpeta3[g,])
    res = t(Udgpeta3.g) %*% iIpk3UtU %*% Udgpeta3.g
    return(res)
  }, mc.cores = mc.cores))

  # M3: Gx1 vector
  k3.2 = 2*k3
  M3 = dg2Vec/2 + eta3.2/k3.2 - Cg3/k3.2 # G x 1

  logf3 = alpha3 * log(beta3) + lgamma(nd2+alpha3) - nlog2pi.2 - lgamma(alpha3) 
  logf3 = logf3 - log(det.Ipk3UtU)/2
  logf3 = logf3 - (nd2 + alpha3)*log(beta3 + M3)
  
  logfMat = cbind(logf1, logf2, logf3)
  colnames(logfMat) = c("logf1", "logf2", "logf3")
  
  return (logfMat)
}

#2. calculate \tilde z, the expected missing value given observed value
Get_tilde_z <- function(
  t_pi, # cluster mixture proportions (pi.OE, pi.UE, pi.NE)
  psi, # vector of model parameters
  dat, # G x n expression data matrix
  tX, # n x (p+1) design matrix
  dg2Vec, # G x 1 vector of d_g^T*d_g
  mc.cores = 4
)
{
  logf = logf123(psi, dat, tX, dg2Vec, mc.cores)
  max_logf = apply(logf, 1, max, na.rm = TRUE)
  
  pi1=t_pi[1]
  pi2=t_pi[2]
  pi3=1-pi1-pi2
  
  t1 = pi1 * exp(logf[,1] - max_logf)
  t2 = pi2 * exp(logf[,2] - max_logf)
  t3 = pi3 * exp(logf[,3] - max_logf)
  total = t1 + t2 + t3
  
  result = cbind(t1, t2, t3)/total
  colnames(result) = c("OE", "UE", "NE")
  
  return(result)
}



##added by Yixin: get\tilda z by Simulated annealing modification EM
Get_tilde_z_SEM <- function(
    t_pi, # cluster mixture proportions (pi.OE, pi.UE, pi.NE)
    psi, # vector of model parameters
    dat, # G x n expression data matrix
    tX, # n x (p+1) design matrix
    dg2Vec, # G x 1 vector of d_g^T*d_g
    mc.cores = 4,
    t_temp ## "temperature" introduced by SEM and iterated at each step
    #  r_coolrate = 0.9
)
{
  logf = logf123(psi, dat, tX, dg2Vec, mc.cores)
  max_logf = apply(logf, 1, max, na.rm = TRUE)
  
  pi1=t_pi[1]
  pi2=t_pi[2]
  pi3=1-pi1-pi2
  
  mypower = 1/t_temp
  
  t1 = pi1 * exp(logf[,1] - max_logf)
  t1 = t1^mypower
  
  t2 = pi2 * exp(logf[,2] - max_logf)
  t2 = t2^mypower
  
  t3 = pi3 * exp(logf[,3] - max_logf)
  t3 = t3^mypower
  
  total = t1 + t2 + t3
  
  result = cbind(t1, t2, t3)/total
  
  colnames(result) = c("OE", "UE", "NE")
  
  return(result)
}


#3. calculate pi
Get_t_pi <- function(
  G,
  tilde_z,
  b
)
{
  denominator = G + sum(b, na.rm=TRUE) - 3
  
  t1 = (sum(tilde_z[,1], na.rm=TRUE) + b[1] - 1) / denominator
  t2 = (sum(tilde_z[,2], na.rm=TRUE) + b[2] - 1) / denominator
  t3 = (sum(tilde_z[,3], na.rm=TRUE) + b[3] - 1) / denominator
  
  res=c(t1, t2, t3)
  names(res)=c("pi.OE", "pi.UE", "pi.NE")
  return(res)
}

#4. calculate the Q function
Qfun <- function(
  # QWL: psi should be the first parameter, which is required by 'optim'
  psi, # parameters to be optimized
  t_pi, # cluster mixture proportions (pi.OE, pi.UE, pi.NE) at t-th iteration of EM algorithm
  tilde_z, # E(z | d, pi^{(t)}, Psi^{(t)}) t-th iteration of EM algorithm
  dat, # G x n expression data matrix
  tX, # n x (p+1) design matrix
  dg2Vec, # G x 1 vector of d_g^T*d_g
  mc.cores = 4,
  const.b # const.b = lgamma(b[1]+b[2]+b[3]) - lgamma(b[1]) - lgamma(b[2]) - lgamma(b[3]) + sum((b-1) * log(t_pi), na.rm = TRUE)
  )
{

  p = ncol(tX) - 1

  logf=logf123(psi, dat, tX, dg2Vec, mc.cores)
  
  result = 0
  result = result + sum(tilde_z[,1] * logf[,1], na.rm = TRUE)
  result = result + sum(tilde_z[,2] * logf[,2], na.rm = TRUE)
  result = result + sum(tilde_z[,3] * logf[,3], na.rm = TRUE)
  result = result + sum(tilde_z[,1] * log(t_pi[1]), na.rm = TRUE)
  result = result + sum(tilde_z[,2] * log(t_pi[2]), na.rm = TRUE)
  result = result + sum(tilde_z[,3] * log(t_pi[3]), na.rm = TRUE)
  
  result = result + const.b
  return ( result)
}


