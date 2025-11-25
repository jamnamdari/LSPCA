# Data Generation
We suggest to call the following libraries before running the examples.

```r
library(gradfps)
library(RSpectra)
library(mvtnorm)
library(Matrix)
library(ggplot2)
library(astsa)
library(lpSolve)
library(lattice)
library(astsa)
library(waveslim)
library(LSPCA)
library(dplR)

```

### Example

The R-code provided below generates the dataframe D in the R package. Detailed description of the process is provided in the manuscript. 


```r
set.seed(1000*3+12)

####################
## Data Generation
####################
p <- 64
n <- 1024
phi_1 <- 1.5
phi_2 <- -.75
a1 <- 1/20
a2 <- -1/1.15
omega <- seq(0,.5, length.out = n/2)
f_xx <- array(0,dim=c(p,p,length(omega)))
c1 <- 3
f_x_omega0 <- (1/c1^2)*arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), n.freq = n/2)$spec
phi_1 <- 1.5 +.05
f20 <- (1/c1^2)*arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), n.freq = n/2)$spec
phi_1 <- 1.5 -.05
f30 <- (1/c1^2)*arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), n.freq = n/2)$spec
phi_1 <- 1.5 +.15
f40 <- (1/c1^2)*arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), n.freq = n/2)$spec
phi_1 <- 1.5 -.15
f50 <- (1/c1^2)*arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), n.freq = n/2)$spec

f_x_omega <- f_x_omega0
for(ell in 1:length(omega)) {
  o_ell <- omega[ell]
  if( !(((o_ell >= .05) & (o_ell <= .25))  ) ) f_x_omega[ell] <- 0
}
f2 <- f20
for(ell in 1:length(omega)) {
  o_ell <- omega[ell]
  if( !(((o_ell >= .05) & (o_ell <= .25))  ) ) f2[ell] <- 0
}
f3 <- f30
for(ell in 1:length(omega)) {
  o_ell <- omega[ell]
  if( !(((o_ell >= .05) & (o_ell <= .25))  ) ) f3[ell] <- 0
}
f4 <- f40
for(ell in 1:length(omega)) {
  o_ell <- omega[ell]
  if( !(((o_ell >= .05) & (o_ell <= .25))  ) ) f4[ell] <- 0
}
f5 <- f50
for(ell in 1:length(omega)) {
  o_ell <- omega[ell]
  if( !(((o_ell >= .05) & (o_ell <= .25))  ) ) f5[ell] <- 0
}

f_e <- 1/4
for(ell in 1:length(omega)){
  f_xx[,,ell] <- diag(f_e, p)
  f_xx[1:5,1:5,ell] <- matrix(c(f_x_omega[ell],(2*1.1)*f_x_omega[ell],1.2*f_x_omega[ell],1.25*f_x_omega[ell], 3*.75*f_x_omega[ell],
                                2*1.1*f_x_omega[ell], 4*1.1*1.1*f_x_omega[ell]+f_e+1*f2[ell], 2*1.1*1.2*f_x_omega[ell],2*1.1*1.25*f_x_omega[ell],2*1.1*3*.75*f_x_omega[ell],
                                1.2*f_x_omega[ell], 2*1.1*1.2*f_x_omega[ell], 1.2*1.2*f_x_omega[ell]+f_e+f3[ell], 1.2*1.25*f_x_omega[ell],1.2*3*.75*f_x_omega[ell],
                                1.25*f_x_omega[ell],2*1.1*1.25*f_x_omega[ell],1.2*1.25*f_x_omega[ell],1.25*1.25*f_x_omega[ell]+f_e+1*f4[ell],1.25*3*.75*f_x_omega[ell],
                                3*.75*f_x_omega[ell],2*1.1*3*.75*f_x_omega[ell],1.2*3*.75*f_x_omega[ell],1.25*3*.75*f_x_omega[ell],9*.75*.75*f_x_omega[ell]+f_e+9*f5[ell]), nrow=5, byrow = TRUE)
}



f_evec11 <- matrix(0, nrow=p, ncol = length(omega))
for(ell in 1:length(omega)){
  f_evec_ell <- eigen(f_xx[,,ell])$vectors[,1]
  f_evec11[,ell] <- f_evec_ell
  #gc()
}


len_freq <- n/2
Coord10 <- matrix(0,nrow = 5, ncol = len_freq)
Coord20 <- matrix(0,nrow = 5, ncol = len_freq)
Coord30 <- matrix(0,nrow = 5, ncol = len_freq)
Coord40 <- matrix(0,nrow = 5, ncol = len_freq)
Coord50 <- matrix(0,nrow = 5, ncol = len_freq)

Coord1 <- matrix(0,nrow = 5, ncol = len_freq)
Coord2 <- matrix(0,nrow = 5, ncol = len_freq)
Coord3 <- matrix(0,nrow = 5, ncol = len_freq)
Coord4 <- matrix(0,nrow = 5, ncol = len_freq)
Coord5 <- matrix(0,nrow = 5, ncol = len_freq)


phi_1 <- 1.5
phi_2 <- -.75
a1 <- 1/20
a2 <- -1/1.15


Xt10 <- arima.sim(list(order=c(4,0,0), ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2)), n = n)*(1/c1)
Xt1 <- pass.filt(Xt10, W=c(0.05, 0.25), type="pass", method = "Butterworth")
temp1 <- arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), main="Autoregression")

phi_1 <- 1.5 +.05
yt20 <- arima.sim(list(order=c(4,0,0), ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2)), n = n)*(1/c1)
yt2 <- pass.filt(yt20, W=c(0.05, 0.25), type="pass", method = "Butterworth")
temp2 <- arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), main="Autoregression")

phi_1 <- 1.5 -.05
yt30 <- arima.sim(list(order=c(4,0,0), ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2)), n = n)*(1/c1)
yt3 <- pass.filt(yt30, W=c(0.05, 0.25), type="pass", method = "Butterworth")
temp3 <- arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), main="Autoregression")

phi_1 <- 1.5 +.15
yt40 <- arima.sim(list(order=c(4,0,0), ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2)), n = n)*(1/c1)
yt4 <- pass.filt(yt40, W=c(0.05, 0.25), type="pass", method = "Butterworth")
temp4 <- arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), main="Autoregression")

phi_1 <- 1.5 -.15
yt50 <- arima.sim(list(order=c(4,0,0), ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2)), n = n)*(1/c1)
yt5 <- pass.filt(yt50, W=c(0.05, 0.25), type="pass", method = "Butterworth")
temp5 <- arma.spec(ar=c(a1+phi_1,a2-a1*phi_1+phi_2,-(phi_1*a2+phi_2*a1), -phi_2*a2), main="Autoregression")

Xt2 <- 2*1.1*Xt1 + 1*yt2 + rnorm(n,0,1/2)
Xt3 <- 1.2*Xt1 + yt3 + rnorm(n,0,1/2)
Xt4 <- 1.25*Xt1 + 1*yt4 + rnorm(n,0,1/2)
Xt5 <- 3*0.75*Xt1 + 3*yt5 + rnorm(n,0,1/2)

X_omega_2 <- cbind(Xt1,Xt2,Xt3,Xt4,Xt5)
X_wn <- rmvnorm(n, sigma = diag(p-(ncol(X_omega_2))))
X_original <- cbind(X_omega_2, X_wn)
#shuffled_index <- sample.int(p)
shuffled_index <- 1:p
D <- X_original[,shuffled_index]

```
