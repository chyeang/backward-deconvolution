plate_model_loglikelihood_evaluation2 <- function(pmode, Q, P, bulkdata, scdata, sccdfs, genegroups, bulksubtypeinds, cellsubtypeinds, cellcelltypeinds, datamode, bdw) {

valintervals <- seq(0.05, 0.95, by = 0.1)
nvalintervals <- length(valintervals)
dI <- 0.1
xmin <- valintervals[1] - 0.5 * dI
xmax <- valintervals[nvalintervals] + 0.5 * dI
epsilon <- 1e-5
epsilon <- 0.01

# Obtain the dimensions of the data.

sizevec <- dim(Q)
if (length(sizevec) == 2) {
  a <- sizevec[1]
  b <- sizevec[2]
} else if (length(sizevec) == 3) {
  a <- sizevec[1]
  b <- sizevec[3]
} else {
  a <- 0
  b <- 0
}

ngenes <- length(genegroups)
nk <- b

ngenegroups <- max(genegroups)
nsubtypes <- max(bulksubtypeinds)
ncelltypes <- max(cellcelltypeinds)
dim_sccdfs <- dim(sccdfs)
ncells <- dim_sccdfs[2]
dim_bulkdata <- dim(bulkdata)
nbsamples <- dim_bulkdata[2]

# Use the old method to rescale Q.
# Rescale each Q column to make the L2 norm equal to the median of the bulk data column L2 norms.

if (pmode == 1) {
  # Calculate the L2 norms of bulk data columns.
  normvecs <- rep(0, nbsamples)
  for (i in 1:nbsamples) {
    val <- norm(bulkdata[, i], type = "2")
    normvecs[i] <- val
  }
  medval <- median(normvecs)
  
  # Rescale each Q column to make the L2 norm equal to the median of the bulk data column L2 norms.
  rQ <- matrix(0, ngenes, nk)
  for (i in 1:nk) {
    vec <- Q[, i]
    val <- norm(vec, type = "2")
    if (val > 0) {
      val <- medval / val
      rQ[, i] <- vec * val
    }
  }
  
}


# Derive two conditional probability tables from the inputs.

# MP: conditional probabilities of components given subtypes.
MP <- matrix(0, nk, nsubtypes)
for (i in 1:nsubtypes) {
  ss <- which(bulksubtypeinds == i)
  val0 <- sum(sum(P[, ss]))
  for (j in 1:nk) {
    val1 <- sum(P[j, ss])
    val <- val1 / val0
    MP[j, i] <- val
    }
}

# GP: conditional probabilities of marker gene expressions given gene groups and components.
# There are three possible ways to derive GP.
# For pmode=1, choose the global scale of Q to maximize log likelihood score of the bulk data.

GP <- array(0, dim = c(ngenegroups, nvalintervals, nk))

# pmode=1: Estimate the pdf of a gene group from Q and bulk expression data.
if (pmode == 1) {
  # Rescaled Q is already computed.
  
  # Obtain the cdf values of Q by comparing Q entries with bulk gene expression data.
  Qcdf <- matrix(0, ngenes, nk)
  for (i in 1:ngenes) {
    for (j in 1:nk) {
      val <- rQ[i, j]
      k <- sum(bulkdata[i, ] <= val)
      val2 <- k / nbsamples
      Qcdf[i, j] <- val2
      }
  }
  # Obtain GP from Qcdf.
  # Apply kernel density estimation.  The bandwith is specified.
  for (i in 1:nk) {
    for (j in 1:ngenegroups) {
      ss <- which(genegroups == j)
      vals <- Qcdf[ss, i]

      if (bdw <= 0) {
        dens <- density(t(vals))
        f <- approx(dens$x, dens$y, xout = valintervals, rule = 2)
      } else {
        dens <- density(t(vals), bw = bdw)
        f <- approx(dens$x, dens$y, xout = valintervals, rule = 2)
      }

      if (sum(f$y) == 0) {
        cnts <- rep(0, nvalintervals)
        for (l in 1:length(vals)) {
          val <- vals[l]
          ds <- abs(valintervals - val)
          k <- which(ds <= min(ds))
          cnts[k] <- cnts[k] + 1 / length(k)
        }
      f$y <- cnts
      }

      f$y <- f$y / sum(f$y)
      GP[j, , i] <- f$y
    }
  }

# pmode=2: Estimate the pdf of a gene group from the distribution in single-cell RNAseq data. 
} else if (pmode == 2) {
  # Obtain the distributions of sccdfs in each gene group and cell type.
  # Also include zero entries.
  # Apply kernel density estimation.
  for (i in 1:ncelltypes) {
    ss <- which(cellcelltypeinds == i)
    for (j in 1:ngenegroups) {
      tt <- which(genegroups == j)
      vals <- sccdfs[tt, ss]
      dens <- density(vals)
      f <- approx(dens$x, dens$y, xout = valintervals, rule = 2)	  
      f$y <- f$y / sum(f$y)
      GP[j, , i] <- f$y
    }
  }

# pmode=3:  Estimate the pdf of a gene group from Q and single-cell RNAseq data.
} else if (pmode == 3) {
  # Calculate the cdf values of the signature matrix in terms of the sccdfs values.
  Qcdf <- matrix(0, ngenes, ncelltypes)
  for (i in 1:ngenes) {
    for (j in 1:ncelltypes) {
      val <- Q[i, j]
      if (val <= 0) {
        cdfval <- 0
      } else {
        k <- sum(scdata[i, ] <= val)
        cdfval <- k / ncells
      }
      Qcdf[i, j] <- cdfval
    }
  }  

  #Obtain GP from Qcdf.
  # Apply kernel density estimation.
  for (i in 1:nk) {
    for (j in 1:ngenegroups) {
      ss <- which(genegroups == j)
      vals <- Qcdf[ss, i]
      dens <- density(t(vals))
      f <- approx(dens$x, dens$y, xout = valintervals, rule = 2)
      f$y <- f$y / sum(f$y)
      GP[j, , i] <- f$y
    }
  }

# pmode=4: Directly load GP in Q.
} else if (pmode == 4) {
  GP <- Q
}

# If GP still has zero entries, then place epsilon to them and renormalize GP.
newGP <- GP
for (i in 1:nk) {
  for (j in 1:ngenegroups) {
    vec <- GP[j, , i]
    ss <- which(vec < epsilon)
    if (length(ss) > 0) {
      vec[ss] <- epsilon
      val <- sum(vec)
      vec <- vec / val
      newGP[j, , i] <- vec
    }
  }
}
GP <- newGP
rm(newGP)

# Evaluate the joint loglikelihood score of the entire data, the loglikelihood scores of single cells, and inferred states of single cells.
# Include all entries including zero entries.
LL <- 0
celllls <- matrix(0, nk, ncells)
cellstates <- rep(0, ncells)

for (i in 1:ncells) {

#  if (i%%1000 == 0) {
#   print(i)
#	flush.console() 
#  }

  subtypeind <- cellsubtypeinds[i]
  
  # Calculate the log likelihood scores of each cell type in each cell.
  subLs <- rep(0, nk)
  for (j in 1:ngenes) {
    val <- sccdfs[j, i]
    g <- genegroups[j]
    if ((datamode == 1 & val >= 0) | (datamode == 2 & val > 0)) {
      qind <- 0
      if (val <= xmin) {
        qind <- 1
      } else if (val >= xmax) {
        qind <- nvalintervals
      } else {
        qind <- ceiling((val - xmin) / dI)
      }
      for (k in 1:nk) {
        p <- GP[g, qind, k]
        subLs[k] <- subLs[k] + log(p)
      }
    }
  }

  # Add the log prior probabilities to each component of subLs.
  jsubLs <- subLs
  for (j in 1:nk) {
    val <- MP[j, subtypeind]
    val <- log(val)
    jsubLs[j] <- jsubLs[j] + val
  }

  maxval <- max(jsubLs)
  cands <- which(jsubLs >= (maxval - 10))
  if (length(cands) == 1) {
    ll <- maxval
  } else {
    v1 <- jsubLs[cands]
    v2 <- v1 - maxval
    v3 <- exp(v2)
    val <- sum(v3)
    ll <- maxval + log(val)
  }

  LL <- LL + ll

  cand <- which(jsubLs >= maxval)
  if (length(cand) > 1) {
    b <- match(cand, cands)
    vals <- MP[cands[b], subtypeind]
    z <- which(vals >= max(vals))
    cand <- cand[z]
  }

  cellstates[i] <- cand
  v3 <- t(jsubLs)
  celllls[, i] <- v3
}

return(list(LL = LL, celllls = celllls, cellstates = cellstates))
}