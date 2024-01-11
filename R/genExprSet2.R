# copied from iCheck Bioconductor package
genExprSet2 = function (ex, pDat, fDat = NULL, annotation = "") 
{
  cn.dat <- colnames(ex)
  rn.pdat <- rownames(pDat)
  aa <- match(cn.dat, rn.pdat)
  if (length(cn.dat) != length(rn.pdat)) {
    warning("Warning: No. of columns of ex=", length(cn.dat), 
        "\n")
    warning("not equalt to that of pDat =", length(rn.pdat), 
        "\n")
    diffxy <- setdiff(cn.dat, rn.pdat)
    if (length(diffxy)) {
      warnings("The sample in ex, but not in pDat are>>\n", diffxy,"\n")
    }
    diffyx <- setdiff(rn.pdat, cn.dat)
    if (length(diffyx)) {
      warnings("The sample in pDat, but not in ex are>>\n", diffyx, "\n")
    }
  }
  if (!any(is.na(aa) == TRUE)) {
    pDat2 <- pDat[aa, , drop = FALSE]
    identical(rownames(pDat2), colnames(ex))
    pDat3 <- as(pDat2, "data.frame")
    aa <- new("AnnotatedDataFrame", data = pDat3)
    exprs <- as(ex, "matrix")
    es.raw <- new("ExpressionSet", exprs = exprs, phenoData = aa, 
                  annotation = annotation)
  }
  else {
    stop("Column names of ex != row names of pDat!\n")
  }
  if (!is.null(fDat)) {
    cn.fdat <- colnames(fDat)
    if (identical(sort(rownames(ex)), sort(rownames(fDat)))) {
      cc <- match(rownames(ex), rownames(fDat))
      Biobase::fData(es.raw) = fDat[cc, , drop = FALSE]
    }
    else {
      stop("Row names of ex != row names of dat.control!\n")
    }
  }
  invisible(es.raw)
}
