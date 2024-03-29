% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scISR.R
\name{scISR}
\alias{scISR}
\title{scISR: Single-cell Imputation using Subspace Regression}
\usage{
scISR(
  data,
  ncores = 1,
  force_impute = FALSE,
  do_fast = TRUE,
  preprocessing = TRUE,
  batch_impute = FALSE,
  seed = 1
)
}
\arguments{
\item{data}{Input matrix or data frame. Rows represent genes while columns represent samples}

\item{ncores}{Number of cores that the algorithm should use. Default value is \code{1}.}

\item{force_impute}{Always perform imputation.}

\item{do_fast}{Use fast imputation implementation.}

\item{preprocessing}{Perform preprocessing on original data to filter out low quality features.}

\item{batch_impute}{Perform imputation in batches to reduce memory consumption.}

\item{seed}{Seed for reproducibility. Default value is \code{1}.}
}
\value{
\code{scISR} returns an imputed single-cell expression matrix where rows represent genes while columns represent samples.
}
\description{
Perform single-cell Imputation using Subspace Regression
}
\details{
scISR performs imputation for single-cell sequencing data. scISR identifies the true dropout values in the scRNA-seq dataset using
hyper-geomtric testing approach. Based on the result obtained from hyper-geometric testing, the original dataset is segregated into two
subsets including training data and imputable data. Next, training data is used for constructing a generalize linear regression model that
is used for imputation on the imputable data.
}
\examples{
{
# Load the package
library(scISR)
# Load Goolam dataset
data('Goolam');
# Use only 500 random genes for example
set.seed(1)
raw <- Goolam$data[sample(seq_len(nrow(Goolam$data)), 500), ]
label <- Goolam$label

# Perform the imputation
imputed <- scISR(data = raw)

if(requireNamespace('mclust'))
{
  library(mclust)
  # Perform PCA and k-means clustering on raw data
  set.seed(1)
  # Filter genes that have only zeros from raw data
  raw_filer <- raw[rowSums(raw != 0) > 0, ]
  pca_raw <- irlba::prcomp_irlba(t(raw_filer), n = 50)$x
  cluster_raw <- kmeans(pca_raw, length(unique(label)),
                        nstart = 2000, iter.max = 2000)$cluster
  print(paste('ARI of clusters using raw data:',
              round(adjustedRandIndex(cluster_raw, label),3)))

  # Perform PCA and k-means clustering on imputed data
  set.seed(1)
  pca_imputed <- irlba::prcomp_irlba(t(imputed), n = 50)$x
  cluster_imputed <- kmeans(pca_imputed, length(unique(label)),
                            nstart = 2000, iter.max = 2000)$cluster
  print(paste('ARI of clusters using imputed data:',
              round(adjustedRandIndex(cluster_imputed, label),3)))
}
}
}
