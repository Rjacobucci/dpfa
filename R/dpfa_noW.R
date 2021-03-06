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
#' @param bias_init fill in
#' @param rk_init fill in
#' @param H_true fill in
#' @param phi_true fill in
#' @param theta_true fill in
#' @param psi_true fill in
#' @param verbose Print progress bar?
#' @return phi_est
#' @return psi_est
#' @return theta_est
#' @return ZZip_est
#' @return ZZip_rowsums
#' @keywords dynamic poisson factor analysis
#' @useDynLib dpfa, .registration=TRUE
#' @import Rcpp
#' @import RcppArmadillo
#' @import psych
#' @import extraDistr
#' @importFrom stats na.omit rgamma runif prcomp rbinom rbeta
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export


dpfa_noW <- function(data,
                 K,
                 N,
                 M,
                 numTime,
                 niter=5000,
                 burnin=floor(niter/2),
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
                 bias_init=1,
                 rk_init=1,
                 verbose=T,
                 H_true,
                 phi_true,
                 theta_true,
                 psi_true){ # inits



  Xmtot = data
  numSample = N
  numTotal = (numTime) * numSample
  timeSelect = 1:numTime
  timeSpan = c(diff(timeSelect),1)

  # hyperparameters
  rk = rep(rk_init,K)
  sk = rgamma(K,sk_a)*sk_b
  bias_0 = rep(bias_init,K)

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
                        init_con_H,
                        rk_vec=rk)#,
                       # phi_tru,
                       # psi_tru,
                        #w_tru,
                        #theta_tru)

  psi_init = inits[[1]]
  phi_init = inits[[2]]
  theta_init = inits[[3]]
  H_init = inits[[5]]


  Psi = psi_true#psi_init
  Phi = phi_true#phi_init
  Theta = theta_true#theta_init
  ZZip = H_true#H_init
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




   # print(summary(t(x_kn)))
    #print(head(t(x_kn)))

    if(verbose==TRUE){
      setTxtProgressBar(pb, b)
    }

     #W_time = sweep(W,2,timeSpan,'/') #-- only needed when non-equidistant time spacing


    Theta_3D = matrix_to_array(Theta, K, numSample,numTime)
    ZZip_3D = matrix_to_array(ZZip,K, numSample,numTime)
    #one.array = array(1,dim(Theta_3D))

   # C_kn <- calcC_kn(ZZip_3D, bias_0, Theta_3D, Phi)

  #  C_kn = matrix(C_kn, nrow=K,byrow=F) # was byrow=F -- checked, is correct

    C_kn = matrix(NA,K,ncol(ZZip))
    for(i in 1:ncol(ZZip)){
      for(k in 1:K){
        C_kn[k,i] = rpois(1,bias_0[k] + Phi[k,] %*% Theta[,i]) # not sure if this is the right truncated poisson
      }
    }


    #print(t(C_kn)[1:30,])
    #print(cor(t(C_kn)[1:30,],t(C_kn)[31:60,]))

    #print(C_kn[,1:10])
   # print(summary(t(C_kn)))

    out2 = mult_cpp(C_kn,Phi,Theta, ZZip)
    C_kk1 = out2[[1]]
    C_k1n = out2[[2]]

    #print(C_kk1/rowSums(C_kk1))

    #Pi_k = matrix(rbeta(K,a0,b0),K,1)
    Pi_k = rbeta(K,shape1 = a0 + rowSums(ZZip),shape2=b0 + numTotal*numTime - rowSums(ZZip));

    z_return=sample_Z(x_kn , p0, rk, Phi,Theta, sk, p1,C_k1n, numSample,Pi_k,ZZip)#Pi_k,
    #ZZip = z_return[[1]]
    ZZip = H_true
    #ZZip = matrix(1, nrow = K, ncol = numTotal) # this causes higher correlations

    Psi = matrix(NA,dim(x_pk)[1],dim(x_pk)[2])
    for (i in 1:dim(x_pk)[2]){
      Psi[,i] <-  rdirichlet(1,(alpha_psi+x_pk)[,i])
    }
    Psi = psi_true

    # chinese restaurant table distribution
    Lk = crt_cpp(x_kn,rk)

    sumbpi = rowSums(ZZip) * log(1-p0)
    rk = rgamma(K, rk_a + Lk)/( rk_b - sumbpi); # from code


    Theta = theta_true#calcTheta(rk, ZZip, x_kn, p0) # rk not the issue, tried rep(1,K) -- x_kn must be issue

    #Phi = matrix(NA,dim(x_pk)[2],dim(x_pk)[2])
    #for (i in 1:dim(x_pk)[2]){
     # Phi[,i] <-  rdirichlet(1,(alpha_psi+C_kk1)[,i])
      Phi <- rdirichlet(3,alpha_phi+C_kk1)
     # print(Phi)
    #}

    Lk = crt_cpp(C_k1n,sk)
    #print(C_k1n[,1:5])
    sumbpi = rowSums(ZZip) * log(1-p0)
    sk = rgamma(K, sk_a + Lk)/( 1/sk_b - sumbpi); # from code

    #print(sk);print(C_k1n[,1:10])

    # !!!!!
    #W = calcW(sk, ZZip, C_k1n, .5)
   # W = W_true


    # Pi_k = rbeta(K,shape1 = a0 + rowSums(ZZip),shape2=b0 + numTotal*numTime - rowSums(ZZip));

    #theta = matrix(rgamma(rk*ZZip+ x_kn),dim(x_kn)) * p0;# from code
    #theta = matrix(rgamma(matrix(rk,K,numTotal) + x_kn,1),dim(x_kn)) # from paper

    # nz = (Xmtot[,1:numTime]!=0)*1
    # lambda = Psi%*%Theta[,1:numTime]
    # RMSE = sum((Xmtot[,1:numTime][nz==1]-lambda[nz==1])^2)/numTotal
    # lambda_avg = lambda/colSums(lambda)
    # negLL = -  sum(Xmtot[nz==1]*log(lambda_avg[nz==1]))/sum(Xmtot[nz==1])

    if(b > burnin){
      count = count + 1
      rett[[count]] <- list(Theta=Theta,Lk=Lk,Psi=Psi,rk=rk, Phi = Phi,  sumbpi = sumbpi, x_kn = x_kn, ZZip = ZZip,#Xmtot = Xmtot,
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



  ZZip.array <- array(NA,c(dim(rett[[1]]$ZZip),(niter-burnin)))
  for(i in 1:(niter-burnin)){
    ZZip.array[,,i] <- rett[[i]]$ZZip
  }

  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }

  ZZip_est = apply(ZZip.array,c(1,2),Mode)
  res$ZZip_est = ZZip_est
  ZZip_rowsums = rowSums(ZZip_est)
  res$ZZip_rowsums = ZZip_rowsums

  res$call <- match.call()
  class(res) <- "dpfa_noW"
  return(res)
}


