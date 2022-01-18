myFisher <- function(data, zero.gene, zero.cell, ncores = 1){
    size = floor(nrow(zero.gene)/ncores)
    parts <- lapply(1:ncores, function(i){
        if (i == ncores){
            c((i - 1)*size + 1, nrow(zero.gene))
        } else {
            c((i - 1)*size + 1, i*size)
        }
    })
    res <- mclapply(parts, mc.cores = ncores, function(part){
        zero.gene <- zero.gene[part[1]:part[2],]
        zg <- zero.gene[sort(rep(1:nrow(zero.gene), nrow(zero.cell))),]
        zc <- zero.cell[rep(1:nrow(zero.cell), nrow(zero.gene)),] - 1 - zg
        zc[zc < 0] <- 0
        m <- zg[, 1] + zc[, 1]
        n <- zg[, 2] + zc[, 2]
        k <- zg[, 1] + zg[, 2]
        x <- zg[, 1]
        list(
            trust = phyper(x, m, n, k),
            distrust = phyper(x - 1, m, n, k, lower.tail = FALSE)
        )
    })
    trust = do.call(what = c, lapply(res, function(r) r$trust))
    distrust = do.call(what = c, lapply(res, function(r) r$distrust))
    trust <- matrix(trust, nrow = nrow(zero.cell), ncol = nrow(zero.gene))
    distrust <- matrix(distrust, nrow = nrow(zero.cell), ncol = nrow(zero.gene))
    trust[data > 0] <- 0
    distrust[data > 0] <- 1
    list(
        trust = trust,
        distrust = distrust
    )
}

linearRegression <- function (i, im, goodData, trust.matrix.pval.inferred, distrust.matrix.pval.inferred, centers, clusters) {
    X <- centers - im[,i]
    distances <- sqrt(diag(t(X) %*% X))
    clus <- which.min(distances)

    goodData <- goodData[,clusters==clus]

    im <- im[,i]
    if(length(im) > 5e4)
    {
        idx.trust <- which(trust.matrix.pval.inferred[,i]  < 0.01)
        idx.distrust <- (i:length(im))[-idx.trust]
    } else if (length(im) > 1e4) {
        idx.trust <- which(im != 0)
        idx.distrust <- which(im == 0)
    } else {
        idx.trust <- which(trust.matrix.pval.inferred[,i]  < 0.01)
        idx.distrust <- which(distrust.matrix.pval.inferred[,i]  < 0.01)
    }

    if ((length(idx.trust)) < 10 | (length(idx.distrust) < 3)) return(im)

    tmp <- abs(cor(goodData[idx.trust,], im[idx.trust]))
    idx.im <- order(tmp, decreasing = TRUE)[1:min(10, ncol(goodData))]

    m <-  lm(y ~ x, list(x = goodData[idx.trust,idx.im], y = im[idx.trust]))

    pred <- predict(m, list(x = goodData[idx.distrust,idx.im]))

    im[idx.distrust] <- pred

    if(length(im) > 1e4)
    {
        rm(goodData)
        gc()
    }

    im
}

#' @title Goolam
#'
#' @description Goolam dataset with data and cell types information.The number of genes is reduced to 10,000.
"Goolam"
