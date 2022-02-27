#'
#'
#' Write intro here
#'
#' @param data Name of the dataset
#' @param K The number of classes
#' @param N Sample size
#' @param M Number of words (variables)
#' @param numTime The number of time points
#' @param niter fill in
#' @param burnin fill in
#' @param init_con fill in
#' @param init_con_H fill in
#' @param alpha_psi fill in
#' @param alpha_phi fill in
#' @param p0  fill in
#' @param p1  fill in
#' @param rk_a  fill in
#' @param rk_b  fill in
#' @param sk_a  fill in
#' @param sk_b fill in
#' @param a0  fill in
#' @param b0  fill in
#' @param verbose Print progress bar?
#' @return phi_est
#' @return psi_est
#' @return W_est
#' @return theta_est
#' @return ZZip_est
#' @return ZZip_rowsums
#' @keywords dynamic poisson factor analysis
#' @useDynLib dpfa, .registration=TRUE
#' @import Rcpp
#' @import RcppArmadillo
#' @import psych
#' @import extraDistr
#' @importFrom stats na.omit rgamma runif prcomp rbinom
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export


dpfa <- function(data,
                 K,
                 N,
                 M,
                 numTime,
                 niter=5000,
                 burnin=2000,
                 init_con=1,
                 init_con_H=3,
                 alpha_psi=1,
                 alpha_phi=1,
                 p0 = 0.5,
                 p1 = 0.5,
                 rk_a = 100,
                 rk_b = 1/rk_a,
                 sk_a = 1e5,
                 sk_b=1/sk_a,
                 a0 = 1,
                 b0 = 1,
                 verbose=T){ # inits

  Xmtot = data
  numSample = N
  numTotal = (numTime) * numSample
  timeSelect = 1:numTime
  timeSpan = c(diff(timeSelect),1)

  # hyperparameters
  rk = rep(1,K)
  sk = rgamma(K,sk_a)*sk_b
  bias_0 = rep(1,K)

  res <- list()


  if(verbose==TRUE){
    pb <- txtProgressBar(min = 0, max = niter, style = 3)
  }

  # initialization

      inits <- init_vals(M,
                        K,
                        N,
                        data,
                        numTime,
                        numTotal,
                        alpha_psi,
                        alpha_phi,
                        init_con_psi=init_con,
                        init_con_phi=init_con,
                        init_con_W=init_con,
                        init_con_theta=init_con,
                        init_con_H)#,
                       # phi_tru,
                       # psi_tru,
                        #w_tru,
                        #theta_tru)

  psi_init = inits[[1]]
  phi_init = inits[[2]]
  theta_init = inits[[3]]
  W_init = inits[[4]]
  H_init = inits[[5]]

  Psi = psi_init
  Phi = phi_init
  Theta = theta_init
  W = W_init
  ZZip = H_init
  #ZZip = matrix(1, nrow = K, ncol = numTotal)

  #z0 = matrix(1,K,numSample)
  # ZZip2 = update_Z(z0,ZZip,numTime)
  # W = ZZip2
  #Pi_k <- BPL(Phi %*% (W *ZZip))

  rett <- list()
  #fits <- matrix(NA,niter,2)
  count=0
  for(b in 1:niter){
    out = mult_cpp(Xmtot,Psi,Theta, ZZip)
    x_pk = out[[1]]
    x_kn = out[[2]]

    if(verbose==TRUE){
      setTxtProgressBar(pb, b)
    }

    W_time = sweep(W,2,timeSpan,'/')


    W_3D = matrix_to_array(W_time, K, numSample,numTime)
    ZZip_3D = matrix_to_array(ZZip,K, numSample,numTime)
    C_kn <- calcC_kn(ZZip_3D, bias_0, W_3D, Phi)
    C_kn = matrix(C_kn, nrow=K)

    out2 = mult_cpp(C_kn,Phi,W, ZZip)
    C_kk1 = out2[[1]]
    C_k1n = out2[[2]]

    res=sample_Z(x_kn , p0, rk, Phi,W_time, sk, p1,C_k1n, numSample,ZZip)#Pi_k,
    ZZip = res[[1]]
    #ZZip = matrix(1, nrow = K, ncol = numTotal) # this causes higher correlations

    Psi = matrix(NA,dim(x_pk)[1],dim(x_pk)[2])
    for (i in 1:dim(x_pk)[2]){
      Psi[,i] <-  rdirichlet(1,(alpha_psi+x_pk)[,i])
    }

    # chinese restaurant table distribution
    Lk = crt_cpp(x_kn,rk)

    sumbpi = rowSums(ZZip) * log(1-p0)
    rk = rgamma(K, rk_a + Lk)/( rk_b - sumbpi); # from code

    Theta = calcTheta(rk, ZZip, x_kn, p0)

    Phi = matrix(NA,dim(x_pk)[2],dim(x_pk)[2])
    for (i in 1:dim(x_pk)[2]){
      Phi[,i] <-  rdirichlet(1,(alpha_psi+C_kk1)[,i])
    }

    Lk = crt_cpp(C_k1n,sk)
    sumbpi = rowSums(ZZip) * log(1-p0)
    sk = rgamma(K, sk_a + Lk)/( 1/sk_b - sumbpi); # from code
    W = calcW(sk, ZZip, C_k1n, 0.5)

    # Pi_k = rbeta(K,shape1 = a0 + rowSums(ZZip),shape2=b0 + numTotal*numTime - rowSums(ZZip));
    # Pi_k = matrix(rbeta(K*numTotal,a0,b0),K,numTotal)
    #theta = matrix(rgamma(rk*ZZip+ x_kn),dim(x_kn)) * p0;# from code
    #theta = matrix(rgamma(matrix(rk,K,numTotal) + x_kn,1),dim(x_kn)) # from paper

    # nz = (Xmtot[,1:numTime]!=0)*1
    # lambda = Psi%*%Theta[,1:numTime]
    # RMSE = sum((Xmtot[,1:numTime][nz==1]-lambda[nz==1])^2)/numTotal
    # lambda_avg = lambda/colSums(lambda)
    # negLL = -  sum(Xmtot[nz==1]*log(lambda_avg[nz==1]))/sum(Xmtot[nz==1])

    if(b > burnin){
      count = count + 1
      rett[[count]] <- list(Theta=Theta,Lk=Lk,Psi=Psi,rk=rk, Phi = Phi,  sumbpi = sumbpi, W = W, x_kn = x_kn, ZZip = ZZip,#Xmtot = Xmtot,
                        C_kn = C_kn, C_k1n = C_k1n, C_kk1=C_kk1)
    }
  }


  # --------------------------
  # summarize results
  # --------------------------

  #------------------- analyze results for each condition ------------------------

  Psi.array <- array(NA,dim=c(M,K,(niter-burnin)))

  for(i in 1:(niter-burnin)){
    Psi.array[,,i] <- rett[[i]]$Psi
  }

  psi_est = apply(Psi.array,c(1,2),mean)
  res$psi_est = psi_est

  Phi.array <- array(NA,dim=c(K,K,(niter-burnin)))
  for(i in 1:(niter-burnin)){
    Phi.array[,,i] <- rett[[i]]$Phi
  }

  phi_est = apply(Phi.array,c(1,2),mean)
  res$phi_est = phi_est

  Theta.array <- array(NA,c(c(K,numTotal),(niter-burnin)))
  for(i in 1:(niter-burnin)){
    Theta.array[,,i] <- rett[[i]]$Theta
  }

  theta_est = apply(Theta.array,c(1,2),mean)
  res$theta_est = theta_est

  W.array <- array(NA,c(dim(rett[[1]]$W),(niter-burnin)))
  for(i in 1:(niter-burnin)){
    W.array[,,i] <- rett[[i]]$W
  }

  W_est = apply(W.array,c(1,2),mean)
  res$W_est = W_est


  ZZip.array <- array(NA,c(dim(rett[[1]]$ZZip),(niter-burnin)))
  for(i in 1:(niter-burnin)){
    ZZip.array[,,i] <- rett[[i]]$ZZip
  }

  ZZip_est = apply(ZZip.array,c(1,2),mean)
  res$ZZip_est = ZZip_est
  ZZip_rowsums = rowSums(ZZip_est)
  res$ZZip_rowsums = ZZip_rowsums

  res$call <- match.call()
  class(res) <- "dpfa"
  return(res)
}


