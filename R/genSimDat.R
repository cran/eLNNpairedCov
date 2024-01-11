####modified on Sept. 25, 2021
#  (1) use 'OE', 'UE', and 'NE' to repalce '1', '2', and '3'
#
####modified on Apr 19, 2019
#  (1) change the function name to 'genSimDat'
####Created by Yixin, 2018/05/08
  #note: get Age scaled


#####Generate hierarchical data#####
##For cluster 1, 2 and 3, they do not have the covariate effects.
# G - number of probes
# n - number of pairs
#  psi=c(alpha1,beta1,k_1, beta10, beta11, beta12, 
#       alpha2,beta2,k_2, beta20, beta21, beta22, 
#       alpha3,beta3,k3, beta31, beta32)
# t_pi=c(pi_1, pi_2)
#   pi_1 - proportion of over-expressed probes
#   pi_2 - proportion of under-expressed probes
# Covariate formula = ~ age + sex

genSimDat<-function(
  G,
  n, 
  psi, 
  t_pi,
  m.age=50, 
  sd.age=5, 
  p.female=0.5) 
{
  ##1. covariate ###
  ##Age: 0-100##
  #Age = scale(round(rnorm(n, mean=m.age, sd=sd.age)), center = TRUE, scale = TRUE)
  Age = rnorm(n, mean=m.age, sd=sd.age)
 
  ##Sex: male,0,female,1##
  Sex=sample(c("F","M"),size=n,replace=TRUE, prob = c(p.female, 1-p.female))
  
  subjid=paste("subj", 1:n, sep="")
  pDat=data.frame(subjid=subjid, Age=Age, Sex=Sex)
  rownames(pDat)=subjid
  
  # need an intercept term
  # QWL: n x (p+1) matrix
  # keep the original design matrix
  #designMat=cbind(rep(1, n), Age, Sex)
  designMat = model.matrix(~Age+Sex, data=pDat)
  
  # (p+1) by n matrix
  X=t(designMat)
  #####2. Hierarchical gene data
  dat = matrix(NA, nrow = G, ncol = n)
  
  pi.OE=t_pi[1]
  pi.UE=t_pi[2]
  pi.NE=1-pi.OE-pi.UE
  
  memGenes=sample(x=c("OE","UE","NE"), size=G, prob=c(pi.OE, pi.UE, pi.NE), replace=TRUE)
  memGenes2=rep("DE", G)
  pos=which(memGenes=="NE")
  if(length(pos))
  {	   
    memGenes2[pos]="NE"
  }
  
  pos.OE=which(memGenes=="OE")
  pos.UE=which(memGenes=="UE")
  pos.NE=which(memGenes=="NE")
  
  if(length(pos.OE) == 0 || length(pos.UE) == 0 || length(pos.NE)==0)
  {
    stop("each probe cluster should have at least 2 probes!\nPlease re-assign values of t_pi!\n")
  }
  
  G.OE=length(pos.OE)
  G.UE=length(pos.UE)
  G.NE=length(pos.NE)
  
  ###cluster 1##
  alpha1 = psi[1]
  beta1 = psi[2]
  k1 = psi[3]
  eta1 = psi[4:6]
  tau_g = rgamma(G.OE, shape=alpha1, rate=beta1)
  for(g in 1:G.OE)
  {
    mu_g = MASS::mvrnorm(
      1,
      mu = exp(apply(X*eta1, 2, sum, na.rm = TRUE)),
      Sigma = diag(k1/tau_g[g], nrow = n)
    )
    dat[pos.OE[g],] = MASS::mvrnorm(
      1,
      mu = mu_g,
      Sigma = diag(1/tau_g[g], nrow = n)
    )
    
  }
  
  ###cluster 2##
  alpha2 = psi[7]
  beta2 = psi[8]
  k2 = psi[9]
  eta2 = psi[10:12]
  tau_g = rgamma(G.UE, shape=alpha2, rate=beta2)
  
  for (g in 1:G.UE) {
    mu_g = MASS::mvrnorm(
      1,
      mu = - exp(apply(X*eta2, 2, sum, na.rm = TRUE)),
      Sigma = diag(k2/tau_g[g], nrow = n)
    )
    
    dat[pos.UE[g],] = MASS::mvrnorm(
      1,
      mu = mu_g,
      Sigma = diag(1/tau_g[g], nrow = n)
    )
  } 
  
  
  
  ###cluster 3###
  alpha3 = psi[13]
  beta3 = psi[14]
  k3 = psi[15]
  # 'eta' is a p x1 vector, not (p+1) x 1 vector
  eta3 = psi[16:17]
  
  tau_g = rgamma(G.NE, shape=alpha3, rate=beta3)
  
  for (g in 1:G.NE) {
    theta_g = MASS::mvrnorm(
      1,
      mu = eta3,
      Sigma = diag(k3/tau_g[g], nrow = 2) #p=2 is the no. of covariate
    )
    # mu_g is a (p+1) x 1 vector
    mu_g=c(0, theta_g)
    dat[pos.NE[g],] = MASS::mvrnorm(
      1,
      # QWL: transpose operation also take times
      #mu = t(X)%*%mu_g,
      mu = designMat%*%mu_g,
      Sigma = diag(1/tau_g[g], nrow = n)
    )
  }
  
  ##ExpressionSet 
  
  probeid=1:G
  
  rownames(dat) =probeid
  colnames(dat)= subjid

  gene=paste("g", 1:G, sep="")
  fDat=data.frame(probeid=probeid, gene=gene, chr=rep(1, G),
                  memGenes.true=memGenes, memGenes2.true=memGenes2)
  rownames(fDat)=probeid
  
  es=genExprSet2(
    ex=dat, 
    pDat=pDat, 
    fDat = fDat, 
    annotation = "")
  
  return(es)
  
}

