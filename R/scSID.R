scSID <- function(data, k=100,h=0.85,n=10){
  
  ## Fano genes for clustering
  Fanodata <- CreateSeuratObject(counts = data)
  Fanodata <- FindVariableFeatures(object = Fanodata, selection.method='vst', nfeatures=dim(data)[1], verbose = F)
  vst <- (Fanodata@assays$RNA@meta.features$vst.variance.standardized)
  den <- density(vst)
  features.vst <- dimnames(data)[[1]][vst > elbow(den$x[which.max(den$y):length(den$x)], den$y[which.max(den$y):length(den$y)])]
  tmp <- data[dimnames(data)[[1]] %in% (features.vst),]
  tmp <- log2(as.matrix(tmp)+1)
  pca <- irlba(t(tmp), nv=min(c(50, dim(tmp)-1))) # More robust no error, contrast to calcul.pca
  pca$pca <-t(pca$d*t(pca$u))
  ##
  knn.res <- Neighbour(pca$pca, pca$pca, k=k,build = "kdtree", cores = 0, checks = 1)
  ##
  dists1<-knn.res$distances
  #首先呢，我需要发现细胞与邻居之间那个间隔大，因为要是大的话，可能是稀有细胞
  
  distance.diff <- (dists1[, -1, drop = FALSE] - dists1[, -ncol(knn.res$distances), drop = FALSE])
  
  dists1_max <- c()
  dists1_rep<-c()
  for ( i in 1:dim(data)[2]) {
    
    dist_rep2<- which(distance.diff[i,]==max(distance.diff[i,]))
    dists1_rep<-c(dists1_rep,dist_rep2)
  }
  # print(dists1_rep)
  rare.cells <- list()
  for(l in 2: k){
    gapMax<-which(dists1_rep==l)
    for (maxgap in gapMax) {
      c<-knn.res$indices[maxgap,1:(l+n)]
      a<-pca$pca[c,]
      d<-dist(a)
      fit1<-hclust(d,method = "single")
      heihts<-max(fit1$height)*h
      bili<-(fit1$height[length(fit1$height)]-fit1$height[length(fit1$height)-2])/fit1$height[length(fit1$height)]
      
      clusterss<-cutree(fit1,h=heihts)
      cell_num<-table(clusterss)
      
      if(cell_num[1]==l&bili>0.2){
        rare_cellss<-knn.res$indices[maxgap, 1:(l)]
        rare.cells[[as.character(maxgap)]] <- rare_cellss
        
      }
      
    }
  }
  return(rare.cells)
}
elbow <- function(x, y){
  n <- length(x)
  firstPoint <- c(x[1], y[1])
  lineVec = c(x[n]-x[1], y[n]-y[1])
  lineVecNorm = lineVec/(sqrt(sum(lineVec^2)))
  vecFromFirst = cbind(x-x[1], y-y[1])
  scalaProd =rowSums(vecFromFirst * cbind(rep(lineVecNorm[1], n), rep(lineVecNorm[2], n)))
  vecFromFirstParallel = outer(scalaProd, lineVecNorm)
  vecToLine = vecFromFirst - vecFromFirstParallel
  distToLine = sqrt(rowSums(vecToLine^2))
  idx = which.max(distToLine)
  return(x[idx])
}
