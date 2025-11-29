


library(gradfps)
library(fps)  # The modified version
library(RSpectra)
library(mvtnorm)
library(Matrix)
library(ggplot2)
library(astsa)
library(lpSolve)
library(lattice)
library(astsa)
library(waveslim)
library(complexplus)
library(MASS)
library(dplR)


Truncate_fun <- function(U, s){
  if(is.null(nrow(U))){
    p <- length(U)
    U_trunc <- matrix(0, nrow = p, ncol = 1)
    index_U_rowSorted <- sort(Mod(U)^2, decreasing = TRUE, index.return = TRUE )
    U_trunc[index_U_rowSorted$ix[1:s]] <- U[index_U_rowSorted$ix[1:s]]
  }else{
    p <- nrow(U)
    d <- ncol(U)
    U_trunc <- matrix(0, nrow = p, ncol = d)
    index_U_rowSorted <- sort(rowSums(Mod(U)^2), decreasing = TRUE, index.return = TRUE )
    U_trunc[index_U_rowSorted$ix[1:s],] <- U[index_U_rowSorted$ix[1:s],]
  }
  return(U_trunc)
}

SOAP_fun <- function(Sigma_hat,U0,s,n_iter){
  ## Initial values of SOAP
  U_trunc <- Truncate_fun(U0, s)
  U_n <- qr.Q(qr(U_trunc))

  ## Sparse orthogonal iteration
  for(t in 1:n_iter){
    if(n_iter == 0) return(U_n)
    V <- Sigma_hat%*%U_n
    V_n <- qr.Q(qr(V))
    rm(V)
    U_trunc <- Truncate_fun(V_n, s)
    rm(V_n)
    U_n <- qr.Q(qr(U_trunc))
    rm(U_trunc)
  }
  return(U_n)
}


LSPCA.f <- function(n,p, f_xx, d, eta, s, n_iter, theta){
  LSPCA0(n,p, f_xx, Pi = matrix(0, nrow = 2*p, ncol=2*p), d, lambda = 0.5 * sqrt(log(p) / n), lr = 0.02, maxiter = 60,
   control = list(fan_maxinc = 10, verbose = 0), eta, s, n_iter, nu=theta)
}


LSPCA <- function(X, d, eta, s, n_iter, theta, ntp = 10){

  p <- ncol(X)
  n <- nrow(X)

  ## Multitaper estimate of the spectral density matrix
  U <- sine.taper(n,ntp)
  X_tp <- apply(U, MARGIN = 2, function(u) u*X, simplify = FALSE)
  F_tp_list <- lapply(X_tp, FUN = function(Y) mvspec(Y,plot = FALSE) )

  len_freq <- n/2
  F_tp1 <- array(0, c(p, p, len_freq))
  for (ell in 1:len_freq) {
    for(j in 1:length(F_tp_list)){
      F_tp1[,,ell] <- F_tp1[,,ell] + F_tp_list[[j]]$fxx[,,ell]
    }
    F_tp1[,,ell] <- F_tp1[,,ell]/length(F_tp_list)
  }
  #plot(Re(F_tp1[1,1,])*(n), type = "l", ylab = " ", main = "Estimated spectral densities for the first ralization", ylim=c(1,400))
  #for(a in 2:p){
  #  lines(Re(F_tp1[a,a,])*n, col= a)
  #}
  f_xx <- F_tp1*n
  rm(U)
  rm(X_tp)
  rm(F_tp_list)
  gc()

  LSPCA0(n,p, f_xx, Pi = matrix(0, nrow = 2*p, ncol=2*p), d, lambda = 0.5 * sqrt(log(p) / n), lr = 0.02, maxiter = 60,
         control = list(fan_maxinc = 10, verbose = 0), eta, s, n_iter, nu=theta)
}


LSPCA0 <- function(n,p, f_xx, Pi, d, lambda, lr = 0.02, maxiter = 60,
                  control = list(fan_maxinc = 10, verbose = 0), eta, s, n_iter, nu){

  if(d == 1){
    return(LSDPC_ADMM_SOAP_PD_fun(n,p, f_xx, Pi,d=1, lambda, lr = 0.02, maxiter = 60,
                                  control = list(fan_maxinc = 10, verbose = 0), eta, s, n_iter, nu))
  }else if(d==2){
    return(LSDPC_ADMM_SOAP_PD2_fun(n,p, f_xx, Pi,d=2, lambda, lr = 0.02, maxiter = 60,
                                   control = list(fan_maxinc = 10, verbose = 0), eta, s, n_iter, nu))
  }else{
    return(LSDPC_ADMM_SOAP_PDd_fun(n,p, f_xx, Pi,ds=d, lambda, lr = 0.02, maxiter = 60,
                                   control = list(fan_maxinc = 10, verbose = 0), eta, s, n_iter, nu))
  }
}


Freq_par_selector <- function(X, f_xx1 = NULL, fxx_evec_s, d){
  if(d == 1){
    return(Freq_par_selector1(X, f_xx1 = NULL, fxx_evec_s))
  }else if(d == 2){
    return(Freq_par_selector2(X, f_xx1 = NULL, fxx_evec_s))
  }else{
    return(Freq_par_selector_d(X, f_xx1 = NULL, fxx_evec_s, ds=d))
  }
}


sparse_par_selector <- function(X, eta, nu=0, d){
  if(d == 1){
    return(sparse_par_selector1(X, eta, nu))
  }else if(d == 2){
    return(sparse_par_selector2(X, eta, nu))
  }else{
    return(sparse_par_selector_d(X, eta, nu, ds=d))
  }
}


##################################################################################
### d=1
##################################################################################




## **FPS+SOAP with initial value estimated from previous frequency component**
## Projection to previous estimate
LSDPC_ADMM_SOAP_PD_fun <- function(n,p, f_xx, Pi, d, lambda, lr = 0.02, maxiter = 60,
                                   control = list(fan_maxinc = 10, verbose = 0), eta, s, n_iter, nu) {


  len_freq <- n/2

  S <- array(0, c(p*2, p*2, len_freq))
  for (ell in 1:len_freq) {
    S[,,ell] <- cbind( rbind(Re(f_xx[,,ell]), -Im(f_xx[,,ell])) , rbind( Im(f_xx[,,ell]), Re(f_xx[,,ell])) )
  }
  gc()


  # Fantope projection and selection
  fxx_grad_evec <- rep(0,p)

  g_xx <- array(0,dim = c(p,p,len_freq)) # deflated spectral density matrix

  for(ell in 1:len_freq){

    ## Fantope projection and selection step
    # Initial value, should be noisy
    if(ell < len_freq){
      if(ell < 10){
        e_ell = RSpectra::eigs_sym(S[,,ell], k=2, which = "LA")
        x0_ell = tcrossprod(e_ell$vectors)

        tryCatch({ fxx_grad_ell <- gradfps::gradfps_prox_benchmark(
          S[,,ell], Pi = matrix(0, nrow = 2*p, ncol=2*p), d=2, lambda , x0 = x0_ell, lr = 0.02, maxiter = 60,
          control = list(fan_maxinc = 10, verbose = 0)
        )
        Re_H <- fxx_grad_ell$projection[1:p,1:p]
        Im_H <- fxx_grad_ell$projection[1:p,(p+1):(2*p)]
        H <- Re_H + Im_H*1i
        gc()
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
        )
      } else { }
    } else { }


    ## SOAP step
    U0 <- eigen(H)$vectors[,1]
    if(ell >= 10){
      U0 <- U_final
      rm(U_final)
    }

    ## Shrinking A
    if(ell < 10){
      f_shrunked <- (1-nu)*f_xx[,,ell] + nu*t(Conj(H))%*%f_xx[,,ell]%*%H
    }else{
      H00 <- fxx_grad_evec_ell%*%t(Conj(fxx_grad_evec_ell))
      f_shrunked <- (1-nu)*f_xx[,,ell] + nu*t(Conj(H00))%*%f_xx[,,ell]%*%H00
    }

    U_final <- SOAP_fun(f_shrunked, U0, s, n_iter)

    lam_ind <- 1

    ## first evec
    fxx_grad_evec_ell <- U_final[,lam_ind[1]]/Re(sqrt(Mod(t(Conj(U_final[,lam_ind[1]]))%*%U_final[,lam_ind[1]])))

    fxx_grad_evec <- cbind(fxx_grad_evec, fxx_grad_evec_ell)


    rm(U0)
    #rm(U_final)
    rm(e_ell)
    rm(x0_ell)
    gc()

    # deflation step
    if(ell < 10){
      g_xx[,,ell] <- f_xx[,,ell] - f_xx[,,ell]%*%H
    } else{
      g_xx[,,ell] <- f_xx[,,ell] - f_xx[,,ell]%*%fxx_grad_evec_ell%*%t(Conj(fxx_grad_evec_ell))
    }

  }

  ## optimization step
  # Coefficients of the linear objective function
  C <- sapply(1:len_freq, function(ell) Re(t(Conj(fxx_grad_evec[,ell+1]))%*%f_xx[,,ell]%*%fxx_grad_evec[,ell+1]))

  # Coefficients of the constraints
  A <- rbind(t(rep(1,len_freq)), diag(len_freq))

  # right hand side of the constraints
  B <- c(eta, rep(1,len_freq))

  # Direction of the constraints
  #constranints_direction  <- c("<=", "<=", "<=", ">=")
  constranints_direction  <- rep("<=", len_freq+1)


  # Find the optimal solution
  optimum <-  lpSolve::lp(direction="max",
                          objective.in = C,
                          const.mat = A,
                          const.dir = constranints_direction,
                          const.rhs = B,
                          all.int = FALSE)
  best_sol <- optimum$solution

  # Obtaining localized PCs
  freq_selector <- matrix(rep(best_sol, p), nrow = p, byrow = TRUE)
  LSDPCA <- freq_selector*(fxx_grad_evec[,-1])
  #LSDPCA_Re <- freq_selector*Re(fxx_grad_evec[,-1])
  #LSDPCA_Im <- freq_selector*Im(fxx_grad_evec[,-1])


  return(list(Evec_series_1 = fxx_grad_evec[,-1], LSDPCA, selected_freqs = best_sol, f_deflated = g_xx))
}

selector <- function(fxx_grad_evec, f_xx, len_freq ,eta,p){

  ## optimization step
  # Coefficients of the linear objective function
  C <- sapply(1:len_freq, function(ell) Re(t(Conj(fxx_grad_evec[,ell]))%*%f_xx[,,ell]%*%fxx_grad_evec[,ell]))

  # Coefficients of the constraints
  A <- rbind(t(rep(1,len_freq)), diag(len_freq))

  # right hand side of the constraints
  B <- c(eta, rep(1,len_freq))

  # Direction of the constraints
  #constranints_direction  <- c("<=", "<=", "<=", ">=")
  constranints_direction  <- rep("<=", len_freq+1)


  # Find the optimal solution
  optimum <-  lpSolve::lp(direction="max",
                 objective.in = C,
                 const.mat = A,
                 const.dir = constranints_direction,
                 const.rhs = B,
                 all.int = FALSE)
  best_sol <- optimum$solution

  # Obtaining localized PCs
  freq_selector <- matrix(rep(best_sol, p), nrow = p, byrow = TRUE)
  LSDPCA <- freq_selector*(fxx_grad_evec)
  #LSDPCA_Re <- freq_selector*Re(fxx_grad_evec)
  #LSDPCA_Im <- freq_selector*Im(fxx_grad_evec)

  return(list(LSDPCA, best_sol))
}



## Frequency parameter selection
Freq_par_selector1 <- function(X, f_xx1 = NULL, fxx_evec_s){
  n <- nrow(X)
  p <- ncol(X)

  if(is.null(f_xx1)){
    U <- sine.taper(n,10)

    X_tp <- apply(U, MARGIN = 2, function(u) u*X, simplify = FALSE)
    F_tp_list <- lapply(X_tp, FUN = function(Y) mvspec(Y,plot = FALSE) )

    len_freq <- n/2
    F_tp1 <- array(0, c(p, p, len_freq))
    for (ell in 1:len_freq) {
      for(j in 1:length(F_tp_list)){
        F_tp1[,,ell] <- F_tp1[,,ell] + F_tp_list[[j]]$fxx[,,ell]
      }
      F_tp1[,,ell] <- F_tp1[,,ell]/length(F_tp_list)
    }
    f_xx1 <- F_tp1*n
    rm(U)
    rm(X_tp)
    rm(F_tp_list)
    gc()
  }


  AIC <- rep(0, len_freq)
  AICc <- rep(0, len_freq)
  n2LL <- rep(0, len_freq)
  BIC <- rep(0, len_freq)

  # All Periodograms
  Periodograms <- mvspec(X, plot = FALSE)

  # sort e-vals
  # compute e-vals
  Evals <- rep(0,len_freq)
  for(ell in 1:len_freq){
    Evals_all <- eigen(f_xx1[,,ell])$values
    Evals[ell] <- Evals_all[1]
  }

  # sort e-vals
  Evals_sorted_index <- order(Evals, decreasing = TRUE)
  len_freq1 <- ceiling(len_freq*.9)
  LL <- rep(0,len_freq1)
  for(ell in 1:(len_freq1)){

    I_set <- Evals_sorted_index[1:ell]
    noise_ind <- setdiff(1:len_freq1, I_set)
    Sig_hat <- apply(Periodograms$fxx[,,noise_ind], c(1,2), mean)
    #plot(eigen(Sig_hat)$values)
    Indicator <- rep(0,len_freq1)
    Indicator[I_set] <- 1
    LL <- rep(0,len_freq1)
    for(j in 1:len_freq1){
      G0 <- f_xx1[,,j]%*%fxx_evec_s[,j]%*%t(Conj(fxx_evec_s[,j]))

      G <- Indicator[j]*G0 + Sig_hat
      G_e <- eigen(G)
      G_eval <- G_e$values[1]
      G_evec <- G_e$vectors[,1]

      G_inv <- (1/G_eval)*G_evec%*%t(Conj(G_evec))
      LL[j] <- Re(-p*log(pi) -sum(log(G_eval)) - sum(diag(G_inv%*%Periodograms$fxx[,,j])))

    }

    n2LL[ell] <- -2*sum(LL)

    AIC[ell] <- -2*sum(LL) + 2*ell

    AICc[ell] <- -2*sum(LL) + 2*ell +(2*(ell)^2 + 2*(ell))/(n-(ell)-1)

    BIC[ell] <- -2*sum(LL) + log(n)*ell
  }

  AICc_min_ind <- order(AICc[1:(len_freq1)], decreasing = FALSE)[1]
  AIC_min_ind <- order(AIC[1:(len_freq1)], decreasing = FALSE)[1]
  n2LL_min_ind <- order(n2LL[1:(len_freq1)], decreasing = FALSE)[1]
  BIC_min_ind <- order(BIC[1:(len_freq1)], decreasing = FALSE)[1]

  return(list(freq_par_AIC = AIC_min_ind, freq_par_AICc = AICc_min_ind, freq_par_BIC = BIC_min_ind))

}


## Sparsity parameter selection
sparse_par_selector1 <- function(X, eta, nu=0){
  ### sparsity selection second round
  n <- nrow(X)

  n0 <- n
  n <- n0/2
  X1 <- X[1:n,]

  U <- sine.taper(n,10)

  X_tp <- apply(U, MARGIN = 2, function(u) u*X1, simplify = FALSE)
  F_tp_list <- lapply(X_tp, FUN = function(Y) mvspec(Y,plot = FALSE) )

  len_freq <- n/2
  F_tp1 <- array(0, c(p, p, len_freq))
  for (ell in 1:len_freq) {
    for(j in 1:length(F_tp_list)){
      F_tp1[,,ell] <- F_tp1[,,ell] + F_tp_list[[j]]$fxx[,,ell]
    }
    F_tp1[,,ell] <- F_tp1[,,ell]/length(F_tp_list)
  }
  plot(Re(F_tp1[1,1,])*(n), type = "l", ylab = " ", main = "Estimated spectral densities for the first realization")
  for(a in 2:p){
    lines(Re(F_tp1[a,a,])*n, col= a)
  }
  f_xx1 <- F_tp1*n
  rm(U)
  rm(X_tp)
  rm(F_tp_list)
  gc()


  TMDist <- rep(0, 20)

  for(k in 1:20){
    print(paste("Round ", k))


    LSDPCA_ADMM_SOAP_Ex1 <- LSDPC_ADMM_SOAP_PD_fun(n,p, f_xx1, Pi <- matrix(0, nrow = 2*p, ncol=2*p), lambda = 0.5 * sqrt(log(p) / n), lr = 0.02, maxiter = 60,
                                                   control = list(fan_maxinc = 10, verbose = 0), eta=eta, s=k, n_iter = 20, nu=nu)



    #### AIC
    ## V3
    # All Periodograms
    X2 <- X[-c(1:n),]
    Periodograms <- mvspec(X2, plot = FALSE)
    Periodograms1 <- mvspec(X1, plot = FALSE)

    # sort e-vals
    # compute e-vals
    Evals <- rep(0,len_freq)
    for(ell in 1:len_freq){
      Evals_all <- eigen(f_xx1[,,ell])$values
      Evals[ell] <- Evals_all[1]
    }

    # sort e-vals
    Evals_sorted_index <- order(Evals, decreasing = TRUE)
    MD <- rep(0,len_freq)
    #ell <- eta
    ell <- ceiling(eta/2)

    I_set <- Evals_sorted_index[1:ell]
    Sig_hat <- apply(Periodograms1$fxx[,,-I_set], c(1,2), mean)
    #plot(eigen(Sig_hat)$values)
    Indicator <- rep(0,len_freq)
    Indicator[I_set] <- 1
    LL <- rep(0,len_freq)
    for(j in 1:len_freq){
      G0 <- f_xx1[,,j]%*%LSDPCA_ADMM_SOAP_Ex1[[1]][,j]%*%t(Conj(LSDPCA_ADMM_SOAP_Ex1[[1]][,j]))

      G <- Indicator[j]*G0 + Sig_hat
      MD[j] <- Mod(sum(diag(solve(G)%*%Periodograms$fxx[,,j])))

    }

    TMDist[k] <- sum(MD)

  }


  min_TMD_ind <- order(TMDist, decreasing = FALSE)[1]
  return(list(sparse_par = min_TMD_ind, TMD = TMDist))

}


##################################################################################
### d=2
##################################################################################

## **FPS+SOAP with initial value estimated from previous frequency component**

## FPS+SOAP


LSDPC_ADMM_SOAP_PD2_fun <- function(n,p, f_xx, Pi, d, lambda, lr = 0.02, maxiter = 60,
                                    control = list(fan_maxinc = 10, verbose = 0), eta, s, n_iter, nu) {
  len_freq <- n/2

  S <- array(0, c(p*2, p*2, len_freq))
  for (ell in 1:len_freq) {
    S[,,ell] <- cbind( rbind(Re(f_xx[,,ell]), -Im(f_xx[,,ell])) , rbind( Im(f_xx[,,ell]), Re(f_xx[,,ell])) )
  }
  gc()


  # Fantope projection and selection
  fxx_grad_evec <- rep(0,p)
  fxx_grad_evec_2 <- rep(0,p)
  fxx_grad_evec_3 <- rep(0,p)
  fxx_grad_evec_4 <- rep(0,p)

  for(ell in 1:len_freq){
    ## Fantope projection and selection step
    # Initial value, should be noisy
    if(ell < len_freq){
      if(ell < 10){
        e_ell = eigs_sym(S[,,ell], k=2, which = "LA")
        x0_ell = tcrossprod(e_ell$vectors)

        tryCatch({ fxx_grad_ell <- gradfps_prox_benchmark(
          S[,,ell], Pi = matrix(0, nrow = 2*p, ncol=2*p), d=4, lambda , x0 = x0_ell, lr = 0.02, maxiter = 60,
          control = list(fan_maxinc = 10, verbose = 0)
        )
        Re_H <- fxx_grad_ell$projection[1:p,1:p]
        Im_H <- fxx_grad_ell$projection[1:p,(p+1):(2*p)]
        H <- Re_H + Im_H*1i
        gc()
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
        )
      } else { }
    } else { }


    ## SOAP step
    U0 <- eigen(H)$vectors[,1:2]
    if(ell >= 10){
      U0 <- U_final
      rm(U_final)
    }
    ## Shrinking A
    if(ell < 10){
      f_shrunked <- (1-nu)*f_xx[,,ell] + nu*t(Conj(H))%*%f_xx[,,ell]%*%H
      #f_shrunked <- f_xx[,,ell] - nu*f_xx[,,ell]%*%H
    }else{
      H00 <- fxx_grad_evec_ell%*%t(Conj(fxx_grad_evec_ell))
      f_shrunked <- (1-nu)*f_xx[,,ell] + nu*t(Conj(H00))%*%f_xx[,,ell]%*%H00
      #f_shrunked <- f_xx[,,ell] - nu*f_xx[,,ell]%*%(diag(p)-fxx_grad_evec_ell%*%t(Conj(fxx_grad_evec_ell)))
    }

    U_final <- SOAP_fun(f_shrunked, U0, s, n_iter)

    lam1 <- Re(t(Conj(U_final[,1]))%*%f_xx[,,ell]%*%U_final[,1])
    lam2 <- Re(t(Conj(U_final[,2]))%*%f_xx[,,ell]%*%U_final[,2])
    lam_ind <- sort(c(lam1,lam2), decreasing = TRUE, index.return=TRUE)$ix
    #lam_ind <- 1:4

    ## first evec
    fxx_grad_evec_ell <- U_final[,lam_ind[1]]/Re(sqrt(Mod(t(Conj(U_final[,lam_ind[1]]))%*%U_final[,lam_ind[1]])))

    fxx_grad_evec <- cbind(fxx_grad_evec, fxx_grad_evec_ell)

    ## second evec
    fxx_grad_evec_ell_2 <- U_final[,lam_ind[2]]/Re(sqrt(Mod(t(Conj(U_final[,lam_ind[2]]))%*%U_final[,lam_ind[2]])))
    fxx_grad_evec_2 <- cbind(fxx_grad_evec_2, fxx_grad_evec_ell_2)


    rm(U0)
    #rm(U_final)
    rm(e_ell)
    rm(x0_ell)
    gc()
    #print(ell)
  }

  ## optimization step
  # Coefficients of the linear objective function
  C <- sapply(1:len_freq, function(ell) Re(t(Conj(fxx_grad_evec[,ell+1]))%*%f_xx[,,ell]%*%fxx_grad_evec[,ell+1]))

  # Coefficients of the constraints
  A <- rbind(t(rep(1,len_freq)), diag(len_freq))

  # right hand side of the constraints
  B <- c(eta, rep(1,len_freq))

  # Direction of the constraints
  #constranints_direction  <- c("<=", "<=", "<=", ">=")
  constranints_direction  <- rep("<=", len_freq+1)


  # Find the optimal solution
  optimum <-  lp(direction="max",
                 objective.in = C,
                 const.mat = A,
                 const.dir = constranints_direction,
                 const.rhs = B,
                 all.int = FALSE)
  best_sol <- optimum$solution

  # Obtaining localized PCs
  freq_selector <- matrix(rep(best_sol, p), nrow = p, byrow = TRUE)
  LSDPCA <- freq_selector*Mod(fxx_grad_evec[,-1])
  LSDPCA_Re <- freq_selector*Re(fxx_grad_evec[,-1])
  LSDPCA_Im <- freq_selector*Im(fxx_grad_evec[,-1])



  return(list(Evec_series_1 = fxx_grad_evec[,-1], Evec_series_2 = fxx_grad_evec_2[,-1], LSDPCA, selected_freqs = best_sol ))
}



## Frequency parameter selection
Freq_par_selector2 <- function(X, f_xx1 = NULL, fxx_evec_s){
  n <- nrow(X)
  p <- ncol(X)

  if(is.null(f_xx1)){
    U <- sine.taper(n,10)

    X_tp <- apply(U, MARGIN = 2, function(u) u*X, simplify = FALSE)
    F_tp_list <- lapply(X_tp, FUN = function(Y) mvspec(Y,plot = FALSE) )

    len_freq <- n/2
    F_tp1 <- array(0, c(p, p, len_freq))
    for (ell in 1:len_freq) {
      for(j in 1:length(F_tp_list)){
        F_tp1[,,ell] <- F_tp1[,,ell] + F_tp_list[[j]]$fxx[,,ell]
      }
      F_tp1[,,ell] <- F_tp1[,,ell]/length(F_tp_list)
    }
    f_xx1 <- F_tp1*n
    rm(U)
    rm(X_tp)
    rm(F_tp_list)
    gc()
  }


  AIC <- rep(0, len_freq)
  AICc <- rep(0, len_freq)
  n2LL <- rep(0, len_freq)
  BIC <- rep(0, len_freq)

  # All Periodograms
  Periodograms <- mvspec(X, plot = FALSE)

  # sort e-vals
  # compute e-vals
  Evals <- rep(0,len_freq)
  for(ell in 1:len_freq){
    Evals_all <- eigen(f_xx1[,,ell])$values
    Evals[ell] <- Evals_all[1]
  }

  # sort e-vals
  Evals_sorted_index <- order(Evals, decreasing = TRUE)
  len_freq1 <- ceiling(len_freq*.9)
  LL <- rep(0,len_freq1)
  for(ell in 1:(len_freq1)){

    I_set <- Evals_sorted_index[1:ell]
    noise_ind <- setdiff(1:len_freq1, I_set)
    Sig_hat <- apply(Periodograms$fxx[,,noise_ind], c(1,2), mean)
    #plot(eigen(Sig_hat)$values)
    Indicator <- rep(0,len_freq1)
    Indicator[I_set] <- 1
    LL <- rep(0,len_freq1)
    for(j in 1:len_freq1){
      G0 <- f_xx1[,,j]%*%fxx_evec_s[[1]][,j]%*%t(Conj(fxx_evec_s[[1]][,j])) +
        f_xx1[,,j]%*%fxx_evec_s[[2]][,j]%*%t(Conj(fxx_evec_s[[2]][,j]))

      G <- Indicator[j]*G0 + Sig_hat
      G_e <- eigen(G)
      G_eval <- G_e$values[1:4]
      G_evec <- G_e$vectors[,1:4]

      G_inv <- G_evec%*%diag(1/G_eval)%*%t(Conj(G_evec))
      LL[j] <- Re(-p*log(pi) -sum(log(G_eval)) - sum(diag(G_inv%*%Periodograms$fxx[,,j])))

    }

    n2LL[ell] <- -2*sum(LL)

    AIC[ell] <- -2*sum(LL) + 2*ell

    AICc[ell] <- -2*sum(LL) + 2*ell +(2*(ell)^2 + 2*(ell))/(n-(ell)-1)

    BIC[ell] <- -2*sum(LL) + log(n)*ell
  }

  AICc_min_ind <- order(AICc[1:(len_freq1)], decreasing = FALSE)[1]
  AIC_min_ind <- order(AIC[1:(len_freq1)], decreasing = FALSE)[1]
  n2LL_min_ind <- order(n2LL[1:(len_freq1)], decreasing = FALSE)[1]
  BIC_min_ind <- order(BIC[1:(len_freq1)], decreasing = FALSE)[1]

  return(list(freq_par_AIC = AIC_min_ind, freq_par_AICc = AICc_min_ind, freq_par_BIC = BIC_min_ind))

}


## Sparsity parameter selection
sparse_par_selector2 <- function(X, eta, nu=0){
  ### sparsity selection second round
  n <- nrow(X)

  n0 <- n
  n <- n0/2
  X1 <- X[1:n,]

  U <- sine.taper(n,10)

  X_tp <- apply(U, MARGIN = 2, function(u) u*X1, simplify = FALSE)
  F_tp_list <- lapply(X_tp, FUN = function(Y) mvspec(Y,plot = FALSE) )

  len_freq <- n/2
  F_tp1 <- array(0, c(p, p, len_freq))
  for (ell in 1:len_freq) {
    for(j in 1:length(F_tp_list)){
      F_tp1[,,ell] <- F_tp1[,,ell] + F_tp_list[[j]]$fxx[,,ell]
    }
    F_tp1[,,ell] <- F_tp1[,,ell]/length(F_tp_list)
  }
  plot(Re(F_tp1[1,1,])*(n), type = "l", ylab = " ", main = "Estimated spectral densities for the first realization")
  for(a in 2:p){
    lines(Re(F_tp1[a,a,])*n, col= a)
  }
  f_xx1 <- F_tp1*n
  rm(U)
  rm(X_tp)
  rm(F_tp_list)
  gc()


  TMDist <- rep(0, 20)

  for(k in 1:20){
    print(paste("Round ", k))


    LSDPCA_ADMM_SOAP_Ex1 <- LSDPC_ADMM_SOAP_PD2_fun(n,p, f_xx1, Pi <- matrix(0, nrow = 2*p, ncol=2*p), lambda = 0.5 * sqrt(log(p) / n), lr = 0.02, maxiter = 60,
                                                    control = list(fan_maxinc = 10, verbose = 0), eta=eta, s=k, n_iter = 20, nu=nu)



    #### AIC
    ## V3
    # All Periodograms
    X2 <- X[-c(1:n),]
    Periodograms <- mvspec(X2, plot = FALSE)
    Periodograms1 <- mvspec(X1, plot = FALSE)

    # sort e-vals
    # compute e-vals
    Evals <- rep(0,len_freq)
    for(ell in 1:len_freq){
      Evals_all <- eigen(f_xx1[,,ell])$values
      Evals[ell] <- Evals_all[1]
    }

    # sort e-vals
    Evals_sorted_index <- order(Evals, decreasing = TRUE)
    MD <- rep(0,len_freq)
    ell <- ceiling(eta/2)

    I_set <- Evals_sorted_index[1:ell]
    Sig_hat <- apply(Periodograms1$fxx[,,-I_set], c(1,2), mean)
    Indicator <- rep(0,len_freq)
    Indicator[I_set] <- 1
    LL <- rep(0,len_freq)
    for(j in 1:len_freq){
      G0 <- f_xx1[,,j]%*%LSDPCA_ADMM_SOAP_Ex1[[1]][,j]%*%t(Conj(LSDPCA_ADMM_SOAP_Ex1[[1]][,j])) +
        f_xx1[,,j]%*%LSDPCA_ADMM_SOAP_Ex1[[2]][,j]%*%t(Conj(LSDPCA_ADMM_SOAP_Ex1[[2]][,j]))

      G <- Indicator[j]*G0 + Sig_hat
      MD[j] <- Mod(sum(diag(solve(G)%*%Periodograms$fxx[,,j])))

    }

    TMDist[k] <- sum(MD)

  }


  min_TMD_ind <- order(TMDist, decreasing = FALSE)[1]
  return(list(sparse_par = min_TMD_ind, TMD = TMDist))

}



##################################################################################
### d=d
##################################################################################

LSDPC_ADMM_SOAP_PDd_fun <- function(n,p, f_xx, Pi, ds, lambda, lr = 0.02, maxiter = 60,
                                    control = list(fan_maxinc = 10, verbose = 0), eta, s, n_iter, nu) {
  len_freq <- n/2

  S <- array(0, c(p*2, p*2, len_freq))
  for (ell in 1:len_freq) {
    S[,,ell] <- cbind( rbind(Re(f_xx[,,ell]), -Im(f_xx[,,ell])) , rbind( Im(f_xx[,,ell]), Re(f_xx[,,ell])) )
  }
  gc()


  # Fantope projection and selection

  fxx_grad_evec_ALL <- as.list(rep(0,ds))
  for(k in 1:ds){
    fxx_grad_evec_ALL[[k]] <- matrix(0,nrow = p,ncol = len_freq)
  }

  for(ell in 1:len_freq){
    ## Fantope projection and selection step
    # Initial value, should be noisy
    if(ell < len_freq){
      if(ell < 10){
        e_ell = RSpectra::eigs_sym(S[,,ell], k=2, which = "LA")
        x0_ell = tcrossprod(e_ell$vectors)

        tryCatch({ fxx_grad_ell <- gradfps::gradfps_prox_benchmark(
          S[,,ell], Pi = matrix(0, nrow = 2*p, ncol=2*p), d=2*ds, lambda , x0 = x0_ell, lr = 0.02, maxiter = 60,
          control = list(fan_maxinc = 10, verbose = 0)
        )
        Re_H <- fxx_grad_ell$projection[1:p,1:p]
        Im_H <- fxx_grad_ell$projection[1:p,(p+1):(2*p)]
        H <- Re_H + Im_H*1i
        gc()
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
        )
      } else { }
    } else { }


    ## SOAP step
    U0 <- eigen(H)$vectors[,1:ds]
    if(ell >= 10){
      U0 <- U_final
      rm(U_final)
    }
    ## Shrinking A
    if(ell < 10){
      f_shrunked <- (1-nu)*f_xx[,,ell] + nu*t(Conj(H))%*%f_xx[,,ell]%*%H
    }else{
      H00 <- fxx_grad_evec_ell%*%t(Conj(fxx_grad_evec_ell))
      for(k in 2:ds){
        H00 <- H00 + fxx_grad_evec_ALL[[k]][,ell-1]%*%t(Conj(fxx_grad_evec_ALL[[k]][,ell-1]))
      }
      f_shrunked <- (1-nu)*f_xx[,,ell] + nu*t(Conj(H00))%*%f_xx[,,ell]%*%H00
    }

    U_final <- SOAP_fun(f_shrunked, U0, s, n_iter)

    lam <- sapply(1:ds, function(i) Re(t(Conj(U_final[,1]))%*%f_xx[,,ell]%*%U_final[,i]) )
    lam_ind <- sort(lam, decreasing = TRUE, index.return=TRUE)$ix

    ## first evec
    fxx_grad_evec_ell <- U_final[,lam_ind[1]]/Re(sqrt(Mod(t(Conj(U_final[,lam_ind[1]]))%*%U_final[,lam_ind[1]])))

    for (k in 1:ds) {
      fxx_grad_evec_ALL[[k]][,ell] <- U_final[,lam_ind[k]]/Re(sqrt(Mod(t(Conj(U_final[,lam_ind[k]]))%*%U_final[,lam_ind[k]])))
    }

    rm(U0)
    #rm(U_final)
    rm(e_ell)
    rm(x0_ell)
    gc()
  }

  ## optimization step
  # Coefficients of the linear objective function
  fxx_grad_evec <- fxx_grad_evec_ALL[[1]]
  C <- sapply(1:len_freq, function(ell) Re(t(Conj(fxx_grad_evec[,ell]))%*%f_xx[,,ell]%*%fxx_grad_evec[,ell]))

  # Coefficients of the constraints
  A <- rbind(t(rep(1,len_freq)), diag(len_freq))

  # right hand side of the constraints
  B <- c(eta, rep(1,len_freq))

  # Direction of the constraints
  #constranints_direction  <- c("<=", "<=", "<=", ">=")
  constranints_direction  <- rep("<=", len_freq)


  # Find the optimal solution
  optimum <-  lpSolve::lp(direction="max",
                          objective.in = C,
                          const.mat = A,
                          const.dir = constranints_direction,
                          const.rhs = B,
                          all.int = FALSE)
  best_sol <- optimum$solution

  # Obtaining localized PCs
  freq_selector <- matrix(rep(best_sol, p), nrow = p, byrow = TRUE)
  LSDPCA <- freq_selector*(fxx_grad_evec)


  return(c(Evec_series_ = fxx_grad_evec_ALL, LSDPCA, best_sol))
}



## Frequency parameter selection
Freq_par_selector_d <- function(X, f_xx1 = NULL, fxx_evec_s, ds){
  n <- nrow(X)
  p <- ncol(X)

  if(is.null(f_xx1)){
    U <- sine.taper(n,10)

    X_tp <- apply(U, MARGIN = 2, function(u) u*X, simplify = FALSE)
    F_tp_list <- lapply(X_tp, FUN = function(Y) mvspec(Y,plot = FALSE) )

    len_freq <- n/2
    F_tp1 <- array(0, c(p, p, len_freq))
    for (ell in 1:len_freq) {
      for(j in 1:length(F_tp_list)){
        F_tp1[,,ell] <- F_tp1[,,ell] + F_tp_list[[j]]$fxx[,,ell]
      }
      F_tp1[,,ell] <- F_tp1[,,ell]/length(F_tp_list)
    }
    f_xx1 <- F_tp1*n
    rm(U)
    rm(X_tp)
    rm(F_tp_list)
    gc()
  }


  AIC <- rep(0, len_freq)
  AICc <- rep(0, len_freq)
  n2LL <- rep(0, len_freq)
  BIC <- rep(0, len_freq)

  # All Periodograms
  Periodograms <- mvspec(X, plot = FALSE)

  # sort e-vals
  # compute e-vals
  Evals <- rep(0,len_freq)
  for(ell in 1:len_freq){
    Evals_all <- eigen(f_xx1[,,ell])$values
    Evals[ell] <- Evals_all[1]
  }

  # sort e-vals
  Evals_sorted_index <- order(Evals, decreasing = TRUE)
  len_freq1 <- ceiling(len_freq*.9)
  LL <- rep(0,len_freq1)
  for(ell in 1:(len_freq1)){

    I_set <- Evals_sorted_index[1:ell]
    noise_ind <- setdiff(1:len_freq1, I_set)
    Sig_hat <- apply(Periodograms$fxx[,,noise_ind], c(1,2), mean)
    #plot(eigen(Sig_hat)$values)
    Indicator <- rep(0,len_freq1)
    Indicator[I_set] <- 1
    LL <- rep(0,len_freq1)
    for(j in 1:len_freq1){
      G0 <- f_xx1[,,j]%*%fxx_evec_s[[1]][,j]%*%t(Conj(fxx_evec_s[[1]][,j]))
      for(k in 2:ds){
        G0 <- G0 + f_xx1[,,j]%*%fxx_evec_s[[k]][,j]%*%t(Conj(fxx_evec_s[[k]][,j]))
      }


      G <- Indicator[j]*G0 + Sig_hat
      G_e <- eigen(G)
      G_eval <- G_e$values[1:4]
      G_evec <- G_e$vectors[,1:4]

      G_inv <- G_evec%*%diag(1/G_eval)%*%t(Conj(G_evec))
      LL[j] <- Re(-p*log(pi) -sum(log(G_eval)) - sum(diag(G_inv%*%Periodograms$fxx[,,j])))

    }

    n2LL[ell] <- -2*sum(LL)

    AIC[ell] <- -2*sum(LL) + 2*ell

    AICc[ell] <- -2*sum(LL) + 2*ell +(2*(ell)^2 + 2*(ell))/(n-(ell)-1)

    BIC[ell] <- -2*sum(LL) + log(n)*ell
  }

  AICc_min_ind <- order(AICc[1:(len_freq1)], decreasing = FALSE)[1]
  AIC_min_ind <- order(AIC[1:(len_freq1)], decreasing = FALSE)[1]
  n2LL_min_ind <- order(n2LL[1:(len_freq1)], decreasing = FALSE)[1]
  BIC_min_ind <- order(BIC[1:(len_freq1)], decreasing = FALSE)[1]

  return(list(freq_par_AIC = AIC_min_ind, freq_par_AICc = AICc_min_ind, freq_par_BIC = BIC_min_ind))

}





## Sparsity parameter selection
sparse_par_selector_d <- function(X, eta, nu=0, ds){
  ### sparsity selection second round
  n <- nrow(X)

  n0 <- n
  n <- n0/2
  X1 <- X[1:n,]

  U <- sine.taper(n,10)

  X_tp <- apply(U, MARGIN = 2, function(u) u*X1, simplify = FALSE)
  F_tp_list <- lapply(X_tp, FUN = function(Y) mvspec(Y,plot = FALSE) )

  len_freq <- n/2
  F_tp1 <- array(0, c(p, p, len_freq))
  for (ell in 1:len_freq) {
    for(j in 1:length(F_tp_list)){
      F_tp1[,,ell] <- F_tp1[,,ell] + F_tp_list[[j]]$fxx[,,ell]
    }
    F_tp1[,,ell] <- F_tp1[,,ell]/length(F_tp_list)
  }
  plot(Re(F_tp1[1,1,])*(n), type = "l", ylab = " ", main = "Estimated spectral densities for the first realization")
  for(a in 2:p){
    lines(Re(F_tp1[a,a,])*n, col= a)
  }
  f_xx1 <- F_tp1*n
  rm(U)
  rm(X_tp)
  rm(F_tp_list)
  gc()


  TMDist <- rep(0, 20)

  for(k in 1:20){
    print(paste("Round ", k))


    LSDPCA_ADMM_SOAP_Ex1 <- LSDPC_ADMM_SOAP_PDd_fun(n,p, f_xx1, Pi <- matrix(0, nrow = 2*p, ncol=2*p), ds, lambda = 0.5 * sqrt(log(p) / n), lr = 0.02, maxiter = 60,
                                                    control = list(fan_maxinc = 10, verbose = 0), eta=eta, s=k, n_iter = 20, nu=nu)



    #### AIC
    ## V3
    # All Periodograms
    X2 <- X[-c(1:n),]
    Periodograms <- mvspec(X2, plot = FALSE)
    Periodograms1 <- mvspec(X1, plot = FALSE)

    # sort e-vals
    # compute e-vals
    Evals <- rep(0,len_freq)
    for(ell in 1:len_freq){
      Evals_all <- eigen(f_xx1[,,ell])$values
      Evals[ell] <- Evals_all[1]
    }

    # sort e-vals
    Evals_sorted_index <- order(Evals, decreasing = TRUE)
    MD <- rep(0,len_freq)
    ell <- ceiling(eta/2)

    I_set <- Evals_sorted_index[1:ell]
    Sig_hat <- apply(Periodograms1$fxx[,,-I_set], c(1,2), mean)
    Indicator <- rep(0,len_freq)
    Indicator[I_set] <- 1
    LL <- rep(0,len_freq)
    for(j in 1:len_freq){
      G0 <- f_xx1[,,j]%*%LSDPCA_ADMM_SOAP_Ex1[[1]][,j]%*%t(Conj(LSDPCA_ADMM_SOAP_Ex1[[1]][,j]))
      for(k in 2:ds){
        G0 <- G0 + f_xx1[,,j]%*%LSDPCA_ADMM_SOAP_Ex1[[k]][,j]%*%t(Conj(LSDPCA_ADMM_SOAP_Ex1[[k]][,j]))
      }


      G <- Indicator[j]*G0 + Sig_hat
      MD[j] <- Mod(sum(diag(solve(G)%*%Periodograms$fxx[,,j])))

    }

    TMDist[k] <- sum(MD)

  }


  min_TMD_ind <- order(TMDist, decreasing = FALSE)[1]
  return(list(sparse_par = min_TMD_ind, TMD = TMDist))

}





















#######################################################################################################################
#######################################################################################################################

selector0 <- function(fxx_grad_evec, f_xx, len_freq ,eta,p){

  ## optimization step
  # Coefficients of the linear objective function
  C <- sapply(1:len_freq, function(ell) Re(t(Conj(fxx_grad_evec[,ell]))%*%f_xx[,,ell]%*%fxx_grad_evec[,ell]))

  # Coefficients of the constraints
  A <- rbind(t(rep(1,len_freq)), diag(len_freq))

  # right hand side of the constraints
  B <- c(eta, rep(1,len_freq))

  # Direction of the constraints
  #constranints_direction  <- c("<=", "<=", "<=", ">=")
  constranints_direction  <- rep("<=", len_freq+1)


  # Find the optimal solution
  optimum <-  lpSolve::lp(direction="max",
                          objective.in = C,
                          const.mat = A,
                          const.dir = constranints_direction,
                          const.rhs = B,
                          all.int = FALSE)
  best_sol <- optimum$solution
  #print(best_sol)

  # Obtaining localized PCs
  freq_selector <- matrix(rep(best_sol, p), nrow = p, byrow = TRUE)
  LSDPCA <- freq_selector*Mod(fxx_grad_evec)
  LSDPCA_Re <- freq_selector*Re(fxx_grad_evec)
  LSDPCA_Im <- freq_selector*Im(fxx_grad_evec)


  g_xx <- array(0,dim = c(p,p,len_freq)) # deflated spectral density matrix

  for(ell in 1:len_freq){
    if(best_sol[ell] != 0){
      g_xx[,,ell] <- f_xx[,,ell] - f_xx[,,ell]%*%fxx_grad_evec[,ell]%*%t(Conj(fxx_grad_evec[,ell]))
    } else{
      g_xx[,,ell] <- f_xx[,,ell]
    }
  }
  return(list(LSDPCA, LSDPCA_Re, LSDPCA_Im, g_xx))
}

