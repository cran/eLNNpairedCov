# created on Oct. 15, 2019
#  (1) function to get gene membership consistent with
#      cluster mixture proportion estimates

# tilde_z - G x 3 matrix; rows are genes; columns are 3 clusters.
#  column 1: for OE genes; column 2: for UE genes; column 3: for NE genes
# t_pi - 3x1 vector of cluster mixture proportion.
#  t_pi[1]: for OE; t_pi[2]: for UE; t_pi[3]: for NE

getMemFunc = function(tilde_z)
{
  G = nrow(tilde_z)
  mem.n = apply(tilde_z, 1, which.max)
  memGenes = rep("NE", G)
  memGenes[which(mem.n==1)] = "OE"
  memGenes[which(mem.n==2)] = "UE"

  invisible(memGenes)
}

getMemFunc2 = function(tilde_z, t_pi)
{
  # obtain position of OE genes
  G = nrow(tilde_z)
  tt = round(G*t_pi)
  G1 = tt[1] # no. of OE genes 
  G2 = tt[2] # no. of UE genes
  G3 = G - G1 - G2 # no. of NE genes

  memGenes=rep("NE", G)

  mat=cbind(tilde_z, 1:G)

  ##
  # OE genes
  ##
  mat1.s=mat[order(mat[,1], decreasing = TRUE),]
  ttpos1 = which(mat1.s[,1]> mat1.s[,2])
  pos1=mat1.s[ttpos1[1:G1],4]
  memGenes[pos1]="OE"

  ##
  # UE genes
  ##
  mat2.s=mat[order(mat[,2], decreasing = TRUE),]
  ttpos2 = which(mat2.s[,2]> mat2.s[,1])
  pos2=mat2.s[ttpos2[1:G2],4]
  memGenes[pos2]="UE"

  invisible(memGenes)
}

