# scISR: Single-cell Imputation using Subspace Regression
scISR performs imputation for single-cell sequencing data. scISR identifies the true dropout values in the scRNA-seq dataset using hyper-geomtric testing approach. Based on the result obtained from hyper-geometric testing, the original dataset is segregated into two including training data and imputable data. Next, training data is used for constructing a generalize linear regression model that is used for imputation on the imputable data.  
# How to install  
- The package can be installed from this repository.  
- Install devtools: `utils::install.packages('devtools')`  
- Install the package using: `devtools::install_github('duct317/scISR')`  
# Example   
## Load the Goolam dataset and perform imputation  
- Load the package: `library(scISR)`  
- Load Goolam dataset: `data('Goolam'); raw <- Goolam$data; label <- Goolam$label`  
- Perform the imputation: `imputed <- scISR(data = raw)`  
## Result assessment
- Perform PCA and k-means clustering on raw data:
```R
library(irlba)
library(mclust)
set.seed(1)
# Filter genes that have only zeros from raw data
raw_filer <- raw[rowSums(raw != 0) > 0, ]
pca_raw <- irlba::prcomp_irlba(t(raw_filer), n = 50)$x
cluster_raw <- kmeans(pca_raw, length(unique(label)),
                      nstart = 2000, iter.max = 2000)$cluster
print(paste('ARI of clusters using raw data:', round(adjustedRandIndex(cluster_raw, label),3)))
```
- Perform PCA and k-means clustering on imputed data:
```R
set.seed(1)
pca_imputed <- irlba::prcomp_irlba(t(imputed), n = 50)$x
cluster_imputed <- kmeans(pca_imputed, length(unique(label)),
                          nstart = 2000, iter.max = 2000)$cluster
print(paste('ARI of clusters using imputed data:', round(adjustedRandIndex(cluster_imputed, label),3)))
```
# Citation:
Duc Tran, Bang Tran, Hung Nguyen, Tin Nguyen (2022). A novel method for single-cell data imputation using subspace regression. <i>Scientific Reports</i>, <b>12</b>, 2697. doi: 10.1038/s41598-022-06500-4 ([link](https://www.nature.com/articles/s41598-022-06500-4))
