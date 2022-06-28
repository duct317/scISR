#' @title scISR: Single-cell Imputation using Subspace Regression
#' @description Perform single-cell Imputation using Subspace Regression
#'
#' @param data Input matrix or data frame. Rows represent genes while columns represent samples
#' @param ncores Number of cores that the algorithm should use. Default value is \code{1}.
#' @param force_impute Always perform imputation.
#' @param do_fast Use fast imputation implementation.
#' @param preprocessing Perform preprocessing on original data to filter out low quality features.
#' @param batch_impute Perform imputation in batches to reduce memory consumption.
#' @param seed Seed for reproducibility. Default value is \code{1}.
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
#' {
#' # Load the package
#' library(scISR)
#' # Load Goolam dataset
#' data('Goolam');
#' # Use only 500 random genes for example
#' set.seed(1)
#' raw <- Goolam$data[sample(seq_len(nrow(Goolam$data)), 500), ]
#' label <- Goolam$label
#'
#' # Perform the imputation
#' imputed <- scISR(data = raw)
#'
#' if(requireNamespace('mclust'))
#' {
#'   library(mclust)
#'   # Perform PCA and k-means clustering on raw data
#'   set.seed(1)
#'   # Filter genes that have only zeros from raw data
#'   raw_filer <- raw[rowSums(raw != 0) > 0, ]
#'   pca_raw <- irlba::prcomp_irlba(t(raw_filer), n = 50)$x
#'   cluster_raw <- kmeans(pca_raw, length(unique(label)),
#'                         nstart = 2000, iter.max = 2000)$cluster
#'   print(paste('ARI of clusters using raw data:',
#'               round(adjustedRandIndex(cluster_raw, label),3)))
#'
#'   # Perform PCA and k-means clustering on imputed data
#'   set.seed(1)
#'   pca_imputed <- irlba::prcomp_irlba(t(imputed), n = 50)$x
#'   cluster_imputed <- kmeans(pca_imputed, length(unique(label)),
#'                             nstart = 2000, iter.max = 2000)$cluster
#'   print(paste('ARI of clusters using imputed data:',
#'               round(adjustedRandIndex(cluster_imputed, label),3)))
#' }
#'}
#' @import stats utils cluster entropy parallel irlba matrixStats PINSPlus
#' @export

scISR <- function(data, ncores = 1, force_impute = FALSE, do_fast = TRUE, preprocessing = TRUE, batch_impute = FALSE, seed = 1){
  if (.Platform$OS.type != "unix") ncores <- 1

  data <- t(data)
  if (is.null(colnames(data))) {
    colnames(data) <- paste0('col', seq(ncol(data)))
  }
  colnames_save <- colnames(data)
  # Transform data into log-scale
  if (max(data) > 100) data <- log2(data+1)

  # Perform gene filtering
  if(preprocessing)
  {
    non.zero.gene <- colSums(data != 0)
    non.zero.prob.gene <- colSums(data != 0)/nrow(data)
    mean.gene <- colMeans(data)
    data <- data[,non.zero.prob.gene > 0.01 & non.zero.gene > 2 & mean.gene > mean(data)/10]
  }

  or.data <- data
  impute_limit <- c(2000, 100)

  if(nrow(data) > 1e4)
  {
    impute_limit <- c(500, 500)
    tmp <- rowSums2(data  != 0)
    if(nrow(data) > 5e4 & !batch_impute) idx <- order(tmp, decreasing = TRUE)[1:1e4 ] else idx <- order(tmp, decreasing = TRUE)[1:5e3 ]
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
  if(force_impute) suck = rep(FALSE, length(keep)) else suck <- distrust.gene.prob > 0.5

  data <- or.data

  rm(fisherRes)
  rm(or.data)
  gc()

  goodData <- data[,keep]
  junk <- data[,suck]
  inferredData <- data[,!keep & !suck]
  trust.matrix.pval.inferred <- trust.matrix.pval[,!keep & !suck]
  distrust.matrix.pval.inferred <- distrust.matrix.pval[,!keep & !suck]

  needImpute <- (ncol(goodData) > impute_limit[1] & ncol(inferredData) > impute_limit[2]) | force_impute
  if (needImpute) {

    set.seed(seed)
    if(nrow(data) > 1e4 | do_fast)
    {
      result <- PINSPlus::PerturbationClustering(data = t(goodData), ncore = ncores)

      clusters <- result$cluster
      names(clusters) <- colnames(goodData)
    } else {
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

    if(batch_impute & nrow(inferredData) > 5e4)
    {
      res <- matrix(ncol = ncol(inferredData), nrow = nrow(inferredData))

      folds <- round(seq(1, nrow(inferredData), length.out = ceiling(nrow(inferredData)/5e4) + 1))
      idx_shuffle <- sample(seq(1, nrow(inferredData)), nrow(inferredData))
      for (i in 2:length(folds)) {
        idx <- idx_shuffle[folds[i-1]:folds[i]]
        inferredData_tmp <- inferredData[idx, ]
        goodData_tmp <- goodData[idx, ]
        centers_tmp = centers[idx, ]
        tmp <- mclapply(as.list(seq(ncol(inferredData))),
                        FUN = linearRegression, im = inferredData_tmp, goodData = goodData_tmp,
                        trust.matrix.pval.inferred = trust.matrix.pval.inferred, distrust.matrix.pval.inferred = distrust.matrix.pval.inferred,
                        centers = centers_tmp, clusters = clusters,
                        mc.cores = ncores)
        tmp <- as.matrix(as.data.frame(tmp))
        res[idx, ] <- tmp
      }

    } else {
      res <- mclapply(as.list(seq(ncol(inferredData))),
                      FUN = linearRegression, im = inferredData, goodData = goodData,
                      trust.matrix.pval.inferred = trust.matrix.pval.inferred, distrust.matrix.pval.inferred = distrust.matrix.pval.inferred,
                      centers = centers, clusters = clusters,
                      mc.cores = ncores)
    }

    res <- as.data.frame(res)
    colnames(res) <- colnames(inferredData)
    res <- as.matrix(cbind(goodData, res, junk))
  } else {
    res <- data
  }
  res <- res[, colnames_save[which(colnames_save %in% colnames(res))] ]
  res <- t(res)
  return(res)
}

