PINSPerturbationClustering <- function (data, Kmax=10, noisePercent="med", iter=200, kmIter=20, PCAFunction = NULL) {
  if (is.null(PCAFunction)){
    pca <- prcomp(data)
  } else {
    pca <- list(
      x = PCAFunction(data)
    )
  }

  # get the partitioning from simply clustering the real data, for consecutive k
  message("Building original connectivity matrices");flush.console()
  origPartition <- getOriginalSimilarity(data=pca$x, clusRange = 2:Kmax)
  origS <- origPartition$origS

  #noise
  noise = getNoise(data, noisePercent)
  message ("Noise set to ", noise)

  # get perturbed similarity
  message("Building perturbed connectivity matrices");flush.console()
  pertS <- getPerturbedSimilarity(data = pca$x, clusRange=2:Kmax, iter = iter, noiseSd = noise, kmIter = kmIter)

  # get discrepancy
  #message("Calculate discrepancy between original and perturbed connectivity matrices")
  Discrepancy <- getPerturbedDiscrepancy (origS = origS, pertS = pertS, clusRange = 2:Kmax)


  ret <- NULL
  ret$k=min(which(Discrepancy$AUC==max(Discrepancy$AUC[2:Kmax])))
  ret$Kmax=10
  ret$groups <- origPartition$groupings[[ret$k]]
  ret$kmRes <- origPartition$kmRes[[ret$k]] #scPINS
  ret$origS <- origS
  ret$pertS <- pertS
  ret$Discrepancy <- Discrepancy

  message("Done. \n");flush.console()

  ret
}

confusionMatrix = function(clusters, classes, label="Missing") {
  rowNA = names(classes)[!names(classes)%in%names(clusters)]
  colNA = names(clusters)[!names(clusters)%in%names(classes)]

  rows=sort(unique(clusters))
  cols=sort(unique(classes))

  if (length(rowNA) >0) {
    rows=c(rows,label)
  }
  if (length(colNA) >0) {
    cols=c(cols,label)
  }

  if (length(rowNA)>0) {
    a <- rep(label,length(rowNA));names(a) <- rowNA
    clusters=c(clusters,a)
  }
  if (length(colNA)>0) {
    b <- rep(label,length(colNA));names(b) <- colNA
    classes=c(classes,b)
  }

  confMat <- data.frame(matrix(0,length(rows),length(cols)))
  rownames(confMat)=rows
  colnames(confMat)=cols

  for (i in rows) {
    cluster <- clusters[clusters==i]

    for (j in cols) {
      confMat[i, j] = sum(classes[names(cluster)]==j)
    }
  }

  confMat
}


distance <- function(x1,x2) {
  sqrt(crossprod(x1 - x2))
}



AddNoiseAll <- function(A, Sd) {
  row=nrow(A)
  col=ncol(A)

  epsilon = matrix(data=rnorm(row*col, mean=0,sd=Sd), nrow=row, ncol=col)

  A <- A+epsilon
  A
}

AddNoiseCol <- function(A, Sds) {
  row=nrow(A)
  col=ncol(A)
  for (i in 1:col) {
    A[,i] = A[,i] + rnorm(row,mean=0,sd=Sds[i])
  }
  A
}


getOriginalSimilarity <- function (data, clusRange) {
  pb <- txtProgressBar(min = min(clusRange)-1, max = max(clusRange)-1, style = 3)

  groupings <- list()

  kmRes <- list() #scPINS

  origS <- list()
  N=nrow(data)
  for (clus in clusRange) {
    setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)

    S=matrix(0,N,N)
    km <- kmeansSS(data, clus, nstart = 2000, iter.max = 2000)
    groupings[[clus]] <- km$cluster
    kmRes[[clus]] <- km #scPINS
    for (j in 1:clus){
      X=rep(0,N);
      X[which(km$cluster==j)]=1
      S=S+X%*%t(X)
    }
    origS[[clus]] <- S
    rownames(origS[[clus]])=rownames(data)
    colnames(origS[[clus]])=rownames(data)
  }

  ret <- NULL
  ret$origS <- origS
  ret$groupings <- groupings
  ret$kmRes <- kmRes #scPINS

  cat("\n")

  ret
}

getPerturbedSimilarity <- function (data, clusRange, iter, noiseSd, kmIter) {
  pertS <- list()
  N=nrow(data)

  for (clus in clusRange) {
    pertS[[clus]] = matrix(0, nrow(data), nrow(data))
    rownames(pertS[[clus]])=rownames(data)
    colnames(pertS[[clus]])=rownames(data)
  }

  pb <- txtProgressBar(min = 0, max = iter, style = 3)

  for (i in 1:iter) {
    #cat(paste("iteration ", i, "/", iter, "\n", sep=''))
    setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)

    tmp = AddNoiseAll(data, noiseSd)
    for (clus in clusRange) {
      S=matrix(0,N,N)
      km <- kmeansSS(tmp, clus, nstart = kmIter)

      for (j in 1:clus){
        X=rep(0,N);
        X[which(km$cluster==j)]=1
        S=S+X%*%t(X)
      }

      pertS[[clus]] = pertS[[clus]] + S/iter
    }
  }

  cat("\n")

  pertS
}


getPerturbedDiscrepancy <- function (origS, pertS, clusRange) {
  ret <- list()

  N<-nrow(origS[[2]])
  diff<-NULL
  for (clus in clusRange) {
    diff[clus]=sum(abs(origS[[clus]]-pertS[[clus]]))
  }

  AUC <- NULL
  entries<-list()
  cdfs<-list()
  for (clus in clusRange) {
    S=abs(origS[[clus]]-pertS[[clus]])
    diag(S) <- 0
    #added -10^(-5) for visual purposes
    A = c(-10^(-5), sort(unique(as.numeric(S))))
    if (max(A)<1) A=c(A,1)
    B<-NULL
    for (i in 1:length(A)) {
      B[i]=sum(S<=A[i])/(N*N)
    }
    entries[[clus]] <- A
    cdfs[[clus]] <- B
    #AUC[clus] = auc(A,B)
    AUC[clus] = areaUnderTheCurve(A, B)
  }


  ret$Diff <- round(diff, digits = 10)
  ret$Entry <- entries
  ret$CDF <- cdfs
  ret$AUC <- round(AUC, digits = 10)
  ret
}

areaUnderTheCurve <- function (x, y) {
  area <- 0;
  for (i in (2:length(x))) {
    area <- area + y[i-1]*(x[i]-x[i-1])
  }
  area
}

kmeansSS <- function(A, k, nstart=20, iter.max=1000) {
  N=nrow(A)
  km <- kmeans(A, centers = k, nstart = nstart, iter.max = iter.max)
}

# getNoise <- function(data, noisePercent="med") {
#   if (noisePercent=="med") {
#     sds <- apply(data, 2, sd)
#     noise = median(sds)
#   } else {
#     noise=sd(as.numeric(as.matrix(data)))*sqrt(noisePercent)
#   }
#   noise
# }

getNoise <- function(data, noisePercent="med") {
  if (noisePercent=="med") {
    sds <- apply(data, 2, sd)
    noise = median(sds)
  } else {
    sds <- apply(data, 2, sd)
    sds = sort(sds)
    ind=round(length(sds)*noisePercent)
    noise=sds[ind]
    #noise=sd(as.numeric(as.matrix(data)))*sqrt(noisePercent)
  }
  noise
}

findMaxHeight <- function (hc, maxK) {
  height <- rev(diff(hc$height))[1:(maxK-1)]
  i=which.max(height)
  i+1
}

pam1 <- function(x,k) list(cluster = pam(x,k, diss=TRUE, cluster.only=TRUE))

getSimilarityFromGrouping <- function(g) {
  N=length(g)
  S=matrix(0,N,N)
  colnames(S)=rownames(S)=names(g)
  for (j in unique(g)){
    X=rep(0,N);
    X[which(g==j)]=1
    S=S+X%*%t(X)
  }
  S
}

clusterUsingHierarchical <- function(orig, pert, Kmax, groupings, method="Rand") {
  ret <- list()
  hcO <- hclust(as.dist(1-orig), method="average")
  hcP <- hclust(as.dist(1-pert), method="average")

  diff<-NULL;agree<-NULL
  for (i in 2:Kmax) {
    gO <- cutree(hcO, i)
    gP <- cutree(hcP, i)
    A=getSimilarityFromGrouping(g = gO) + getSimilarityFromGrouping(g = gP);
    A=A/2; diff[i]=1-(sum(A==0)+sum(A==1))/nrow(orig)^2

    #km <- structure(list(cluster=gO), class="kmeans")
    l = groupings; l[[length(l)+1]]=gO
    agree[i]=clusterAgreement(l, nrow(orig))
  }
  range=which(diff==min(diff[2:Kmax]))

  ret$diff=diff
  ret$agree=agree
  ret$k <- range[which.max(agree[range])]

  ret$gO <- cutree(hcO, ret$k)
  ret$gP <- cutree(hcP, ret$k)
  ret$km <- structure(list(cluster=ret$gO), class="kmeans")

  ret
}

clusterUsingPAM <- function(orig, pert, Kmax, groupings, method="Rand") {
  ret <- list()


  diff<-NULL;agree<-NULL
  for (i in 2:Kmax) {
    gO <- pam1(x = 1-orig, k = i)$cluster
    gP <- pam1(x = 1-pert, k = i)$cluster
    A=getSimilarityFromGrouping(g = gO)+getSimilarityFromGrouping(g = gP); A=A/2; diff[i]=1-(sum(A==0)+sum(A==1))/nrow(orig)^2
    #km <- structure(list(cluster=gP), class="kmeans")
    l = groupings; l[[length(l)+1]]=gP
    agree[i]=clusterAgreement(l, nrow(orig))
  }
  range=which(diff==min(diff[2:Kmax]))

  ret$diff=diff
  ret$agree=agree
  ret$k <- range[which.max(agree[range])]

  ret$gO <- pam1(x = 1-orig, k = ret$k)$cluster
  ret$gP <- pam1(x = 1-pert, k = ret$k)$cluster
  ret$km <- structure(list(cluster=ret$gP), class="kmeans")

  ret
}

clusterAgreement <- function(l, N, method="Rand") {
  A <- matrix(0, N, N)
  for (group in l) {
    A <- A + getSimilarityFromGrouping(group)
  }
  A=A/(length(l))
  ret=(sum(A==0)+sum(A==1))/(N^2)

  ret
}
