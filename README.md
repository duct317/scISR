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
- Log transform the raw data if necessary: `if(max(data) > 100) raw <- log2(raw + 1)`  
- Perform the imputation: `imputed <- scISR(data = raw)`  
## Result assessment
- Perform PCA and k-means clustering on raw data:
```R
library(irlba)
library(mclust)
set.seed(1)
pca_raw <- irlba::irlba(t(raw), nv = 20)$u
cluster_raw <- kmeans(pca_raw, length(unique(label)), nstart = 2000, iter.max = 2000)$cluster
print(paste('ARI of clusters using raw data:', round(adjustedRandIndex(cluster_raw, label),2)))
```
- Perform PCA and k-means clustering on imputed data:
```R
set.seed(1)
pca_imputed <- irlba::irlba(t(imputed), nv = 20)$u
cluster_imputed <- kmeans(pca_imputed, length(unique(label)), nstart = 2000, iter.max = 2000)$cluster
print(paste('ARI of clusters using imputed data:', round(adjustedRandIndex(cluster_imputed, label),2)))
```
