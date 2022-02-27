init_vals <- function(M,
                      K,
                      N,
                      data,
                      numTime,
                      numTotal,
                      alpha_psi,
                      alpha_phi,
                      init_con_psi,
                      init_con_phi,
                      init_con_W,
                      init_con_theta,
                      init_con_H,
                      phi_tru,
                      psi_tru,
                      w_tru,
                      H_tru,
                      theta_tru){

  #------ psi ------
  if (init_con_psi == 1){
    psi_init = matrix(rgamma(M*K,shape=alpha_psi),M,K) #randg(alpha_psi,P,K)
    for (i in 1:K){
      psi_init[,i] = psi_init[,i]/sum(psi_init[,i])
    }
    psi_init = matrix(1/M,M,K)
  } else if (init_con_psi == 2){
    psi_init = psi_tru
  } else if (init_con_psi == 3){
    pca=prcomp(t(data))
    psi_init = abs(pca$rotation[,1:K])
    psi_init = psych::principal(t(data),K)$loadings
    for (i in 1:K){
      psi_init[,i] = as.numeric(psi_init[,i]/sum(psi_init[,i]))
    }
  }

  #------ phi ------
  if (init_con_phi ==1){
    #K1 = K
    #K2 = K
    #phi_init = matrix(rgamma(K1*K2,shape=alpha_phi),K1,K2) #randg(alpha_psi,P,K)
    oness = matrix(1,K,K)
    diag(oness) = 0
    matt = diag(K)*.8 + oness*(.2/(K-1))
    phi_init=rdirichlet(3,matt)
   # for (i in 1:K2){
   #   phi_init[,i] = phi_init[,i]/sum(phi_init[,i])
    #}
  } else if (init_con_phi ==2){
    phi_init = phi_tru
  }

  #------ W ------
  if (init_con_W ==1){
    W_init<-array(dim=c(K,numTotal))
    for (n in 1:numTotal){
      for (k in 1:K){
        W_init[k,n] = rgamma(1,3,0.5)
      }
    }
  } else if (init_con_W ==2){
    W_init = w_tru
  }

  #------ theta ------
  if (init_con_theta ==1){
    theta_init<-array(dim=c(K,numTotal))
    for (n in 1:numTotal){
      for (k in 1:K){
        theta_init[k,n] = rgamma(1,3,0.5)
      }
    }
  } else if (init_con_theta ==2){
    theta_init = theta_tru
  }

  #------ H ------
  if (init_con_H ==1){
    H_init<-array(rbinom(N*numTime,1, 0.5),dim=c(K,numTotal))
  } else if (init_con_H ==2){
    H_init = H_tru
  } else if (init_con_H ==3){
    H_init = array(1, dim=c(K,numTotal))
  }

  inits = list(psi_init, phi_init, theta_init, W_init, H_init)
}
