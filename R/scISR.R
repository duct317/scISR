#' @title scISR: Single-cell Imputation using Subspace Regression
#' @description Perform single-cell Imputation using Subspace Regression
#'
#' @param data Input matrix or data frame. Rows represent genes while columns represent samples
#' @param ncores Number of cores that the algorithm should use. Default value is \code{1}.
#' @param ncores Seed for reproducibility. Default value is \code{1}.
#
#' @details
#' scISR performs imputation for single-cell sequencing data. scISR identifies the true dropout values in the scRNA-seq dataset using
#' hyper-geomtric testing approach. Based on the result obtained from hyper-geometric testing, the original dataset is segregated into two
#' subsets including training data and imputable data. Next, training data is used for constructing a generalize linear regression model that
#' is used for imputation on the imputable data.
#'
#' @return
#' \code{scISR} returns an imputed single-cell expression matrix where rows represent genes while columns represent samples.
#'
#' @examples
#'
#' # Load the package
#' library(scISR)
#' # Load Goolam dataset
#' data('Goolam'); raw <- Goolam$data; label <- Goolam$label
#' # Log transform the raw data if necessary
#' if(max(data) > 100) raw <- log2(raw + 1)
#'
#' # Perform the imputation
#' imputed <- scISR(data = raw)
#'
#' # Perform PCA and k-means clustering on raw data
#' set.seed(1)
#' # Filter genes that have only zeros from raw data
#' raw_filer <- raw[rowSums(raw != 0) > 0, ]
#' pca_raw <- irlba::irlba(t(raw_filer), nv = 20)$u
#' cluster_raw <- kmeans(pca_raw, length(unique(label)), nstart = 2000, iter.max = 2000)$cluster
#' print(paste('ARI of clusters using raw data:', round(adjustedRandIndex(cluster_raw, label),3)))
#'
#' # Perform PCA and k-means clustering on imputed data
#' set.seed(1)
#' pca_imputed <- irlba::irlba(t(imputed), nv = 20)$u
#' cluster_imputed <- kmeans(pca_imputed, length(unique(label)), nstart = 2000, iter.max = 2000)$cluster
#' print(paste('ARI of clusters using imputed data:', adjustedRandIndex(cluster_imputed, label)))
#'
#' @import stats utils cluster entropy parallel irlba
#' @export

scISR <- function(data, ncores = 1, seed = 1){
  if (.Platform$OS.type != "unix") ncores <- 1

  data <- t(data)
  # Transform data into log-scale
  if (max(data) > 100) data <- log2(data+1)

  # Perform gene filtering
  non.zero.gene <- colSums(data != 0)
  non.zero.prob.gene <- colSums(data != 0)/nrow(data)
  mean.gene <- colMeans(data)
  data <- data[,non.zero.prob.gene > 0.01 & non.zero.gene > 2 & mean.gene > mean(data)/10]

  or.data <- data

  if(nrow(data) > 1e4)
  {
    tmp <- rowSums2(data  != 0)
    if(nrow(data) > 5e4) idx <- order(tmp, decreasing = T)[1:1e4 ] else idx <- order(tmp, decreasing = T)[1:5e3 ]
    data <- data[idx, ]
  }

  zero.cell <- matrix(ncol = 2, nrow = nrow(data))
  zero.cell[,1] <- rowSums(data == 0)
  zero.cell[,2] <- rowSums(data > 0)

  zero.gene <- matrix(ncol = 2, nrow = ncol(data))
  zero.gene[,1] <- colSums(data == 0)
  zero.gene[,2] <- colSums(data > 0)

  # Perform hyper geometric testing

  fisherRes <- myFisher(data, zero.gene, zero.cell, ncores = ncores)
  trust.matrix.pval <- fisherRes$trust
  distrust.matrix.pval <- fisherRes$distrust

  colnames(trust.matrix.pval) <- colnames(distrust.matrix.pval) <- colnames(data)
  rownames(trust.matrix.pval) <- rownames(distrust.matrix.pval) <- rownames(data)

  trust.gene.prob <- colSums(trust.matrix.pval < 0.01)/nrow(trust.matrix.pval)
  distrust.gene.prob <- colSums(distrust.matrix.pval < 0.01)/nrow(distrust.matrix.pval)

  keep <- trust.gene.prob == 1
  suck <- distrust.gene.prob > 0.5

  data <- or.data

  rm(fisherRes)
  rm(or.data)
  gc()

  goodData <- data[,keep]
  junk <- data[,suck]
  inferredData <- data[,!keep & !suck]
  trust.matrix.pval.inferred <- trust.matrix.pval[,!keep & !suck]
  distrust.matrix.pval.inferred <- distrust.matrix.pval[,!keep & !suck]

  needImpute <- ncol(goodData) > 2000 & ncol(inferredData) > 100
  if (needImpute) {
    if(nrow(data) > 1e4)
    {
      set.seed(seed)
      result <- PINSPlus::PerturbationClustering(data = t(goodData), ncore = ncores)

      clusters <- result$cluster
      names(clusters) <- colnames(goodData)
    } else {
      set.seed(seed)
      result <- PINSPerturbationClustering(data = t(goodData), Kmax = 5,
                                       PCAFunction = function(x){irlba::prcomp_irlba(x, n = 10)$x},
                                       iter = 100, kmIter = 100)
      clusters <- result$groups
      names(clusters) <- colnames(goodData)
    }

    centers <- matrix(ncol = length(unique(clusters)), nrow = nrow(goodData))
    for (clus in sort(unique(clusters))) {
      centers[,clus] <- rowMeans(as.matrix(goodData[,clusters==clus]))
    }

    res <- mclapply(as.list(seq(ncol(inferredData))),
                      FUN = linearRegression, im = inferredData, goodData = goodData,
                      trust.matrix.pval.inferred = trust.matrix.pval.inferred, distrust.matrix.pval.inferred = distrust.matrix.pval.inferred,
                      centers = centers, clusters = clusters,
                      mc.cores = ncores)

    newData <- as.data.frame(res)
    fi.res <- as.matrix(cbind(goodData, newData, junk))
  } else {
    message("No need to impute")
    fi.res <- data
  }

  fi.res <- t(fi.res)
  return(fi.res)
}
