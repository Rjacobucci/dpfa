sample_Z = function(x_kn,
                    p0,
                    rk,
                    Phi,
                    W,
                    sk,
                    p1,
                    C_k1n,
                    Pi_k,
                    numSample,
                    y_kn){
  numTime <-  ncol(x_kn)/numSample
  K <- nrow(x_kn)
  t1_ix <- seq(1, numSample*numTime, numTime)
  c0 <- C_k1n[,t1_ix]
  C_k1n <-  cbind(C_k1n[, 2:ncol(C_k1n)], 0)
  tn_ix <- seq(numTime, numSample*numTime, numTime)
  C_k1n[,tn_ix] <-  0
  
  lix <- (x_kn == 0 & C_k1n == 0 & y_kn == 0)
  ind <- which(lix, arr.ind=T)
  rix <- ind[,1]
  cix <- ind[,2]
  
  BPL <- function(lambda){
    1-exp(-lambda)
  }
  
  lambda = Phi %*% W;
  Pi <-  BPL(lambda)
  
  
  p_1 = Pi[lix]*(( 1 - p0 )^rk[rix]) * ((1-p1)^sk[rix])
  idx = which(cix %in% tn_ix)
  p_1[idx] = p_1[idx]/(( 1 - p1)^sk[rix[idx]])
  p_0 <-  1 - Pi[lix]
  
  ZZip <- matrix(1, nrow=nrow(x_kn), ncol=ncol(x_kn))
  ZZip[lix] <-   (p_1/( p_1 + p_0 ) ) > runif( length( rix ))
  pZZip <-  matrix(1, nrow=nrow(x_kn), ncol=ncol(x_kn))
  pZZip[lix] <-  (p_1/( p_1 + p_0 ))
  lix <-  (c0 == 0);
  ind <- which(lix, arr.ind=T)
  rix <- ind[,1]
  cix <- ind[,2]
  p_1 <- Pi_k[rix]*( ( 1 - p1 )^sk[rix])
  p_0 <-  1 - Pi[rix];
  z0 <-  matrix(1, nrow=nrow(c0), ncol=ncol(c0))
  z0[lix] <- (p_1/( p_1 + p_0 )) > runif( length( rix ))
  pZ0 <-  matrix(1, nrow=nrow(c0), ncol=ncol(c0))
  pZ0[lix] <- (p_1/( p_1 + p_0 )) 
  pZZip <- cbind(pZ0, pZZip)
  list(ZZip=ZZip, z0=z0, pZZip=pZZip)
}
