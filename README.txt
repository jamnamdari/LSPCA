
# LSPCA

`LSPCA` is an R-package that provides an implementation of the localized and sparse principal components of multivariate time series in the frequency domain. The algorithm is based on the paper
[Localized Sparse Principal Component Analysis in Frequency Domain](https://arxiv.org/abs/2408.08177#)
by Jamshid Namdari, Amita Manatunga,Fabio Ferrarelli, and Robert Krafty.

## Installation

Currently the `LSPCA` has not been submitted to CRAN, but it can be installed dirrectly from GitHub:

```r
library(devtools)
install_github("jamnamdari/LSPCA")
```


The `LSPCA()` function requires the user to install and call the following libraries

+ library(astsa)
+ library(waveslim)
+ library(ggplot2)
+ library(RSpectra)
+ library(lpSolve)
+ library(gradfps)

Note that `gradfps` is not available on CRAN and can be installed from their authors` guithub page.

```r
library(devtools)
install_github("yixuan/gradfps")
```
The other libraries are aveilabe on CRAN and can be installed by calling `install.packages("name_of_the_package")`.

Note that `gradfps` would need a C++ compiler that supports the C++11 standard.
The authors of `gradfps` have suggested that "For best performance, it is strongly suggested linking your R to the
[OpenBLAS](https://www.openblas.net/) library for matrix computation.
You can achieve this with the help of the
[ropenblas](https://prdm0.github.io/ropenblas/) package".

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

The R-code provided below reproduces Figure 2 of the main manuscript.

We first generate a realization of a 64-dimensional time series.

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
X1 <- X_original[,shuffled_index]


U <- sine.taper(n,20)

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
f_xx1 <- F_tp1*n
rm(U)
rm(X_tp)
rm(F_tp_list)
#gc()

plot(omega, f_xx1[1,1,1:512], type="l", ylim = c(0,30), ylab="", xlab="Frequency")
for (k in 2:64) {
  points(omega, f_xx1[k,k,1:512], type="l", col=k)
}
abline(v=c(0.05,0.25))

```


Next we estimate the leading principal subspace of the underlying spectral density matrices over the frequency range [0,32] Hz.

```r
##################################################
## Estimate of 1-dimensional principal subspaces
##################################################
nu_v <- c(0,.2,.4,.6,.8,1)
k <- 1
nu <- nu_v[k]

LSDPCA_ADMM_SOAP_Ex1 <- LSPCA(n,p, f_xx1, lambda = 0.5 * sqrt(log(p) / n), d=1, lr = 0.02, maxiter = 60,
                                                 control = list(fan_maxinc = 10, verbose = 0), eta=200, s=5, n_iter = 20, nu=nu_v[k])


#gc()


k <- 4
nu <- nu_v[k]

LSDPCA_ADMM_SOAP_Ex2 <- LSPCA(n,p, f_xx1, lambda = 0.5 * sqrt(log(p) / n), d=1, lr = 0.02, maxiter = 60,
                                                 control = list(fan_maxinc = 10, verbose = 0), eta=200, s=5, n_iter = 20, nu=nu_v[k])


#gc()


f_MT_evec11 <- matrix(0, nrow=p, ncol = length(omega))
for(ell in 1:length(omega)){
  f_MT_evec_ell <- eigen(f_xx1[,,ell])$vectors[,1]
  f_MT_evec11[,ell] <- f_MT_evec_ell
  #gc()
}


```

## Plots

### Bottom Left panel of Figure 2

Bottom Left panel of Figure 2 of the main manuscript can be reproduced by the follwoing code.

```r
###########################
## Plots
###########################

library(latex2exp)


xlab <- "Hz"
ylab <- "Coordinate"
legend_title = ""
asp = .5 #0.2
bar_height = 10/2
font_size = 20/2



Localized_Est <- selector(LSDPCA_ADMM_SOAP_Ex1[[1]],f_xx1,n/2,(2/5)*512,p)
evecs <- t(Mod(Localized_Est[[1]]))
v = as.matrix(evecs)
lo = min(v)
hi = max(v)
#lo <- 0
#hi <- .9
r = max(abs(c(lo, hi)))
gdat = data.frame(x = as.integer(row(v)), y = as.integer(col(v)),
                  z = as.numeric(v))
ngrid = 1001
pal = c("#67000d", "#a50f15", "#cb181d", "#ef3b2c", "#fb6a4a",
        "#fc9272", "#fcbba1", "#fee0d2", "#ffffff", "#deebf7",
        "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5",
        "#08519c", "#08306b")
col_pal = colorRampPalette(pal)(ngrid)
col_val = seq(-r, r, length.out = ngrid)
lo_ind = findInterval(lo, col_val)
hi_ind = findInterval(hi, col_val)
colors = col_pal[lo_ind:hi_ind]

ggplot(gdat, aes(x = x, y = y, fill = z)) + geom_tile() +
  scale_x_continuous(xlab, expand = c(0, 0)) + scale_y_reverse(ylab,
                                                               breaks = 1:ncol(v), expand = c(0, 0)) + scale_fill_gradientn(legend_title,
                                                                                                                            colors = colors, limits = c(0, 1)) + guides(fill = guide_colorbar(barheight = bar_height)) +
  theme_bw(base_size = font_size) + theme(aspect.ratio = asp) + labs(title = TeX("$\\theta=0$")) +
  scale_y_continuous(name = "Coordinate", breaks=seq(0,p,10)) +
  scale_x_continuous(name = "Hz", breaks = seq(0,.5,by=.05)*(n), labels=function(x) x*64/n)

```

### Bottom Right panel of Figure 2

Bottom right panel of Figure 2 of the main manuscript can be reproduced by the following code.

```r
###########################
## Plots
###########################

library(latex2exp)


xlab <- "Hz"
ylab <- "Coordinate"
legend_title = ""
asp = .5 #0.2
bar_height = 10/2
font_size = 20/2



Localized_Est <- selector(LSDPCA_ADMM_SOAP_Ex2[[1]],f_xx1,n/2,(2/5)*512,p)
evecs <- t(Mod(Localized_Est[[1]]))
v = as.matrix(evecs)
lo = min(v)
hi = max(v)
#lo <- 0
#hi <- .9
r = max(abs(c(lo, hi)))
gdat = data.frame(x = as.integer(row(v)), y = as.integer(col(v)),
                  z = as.numeric(v))
ngrid = 1001
pal = c("#67000d", "#a50f15", "#cb181d", "#ef3b2c", "#fb6a4a",
        "#fc9272", "#fcbba1", "#fee0d2", "#ffffff", "#deebf7",
        "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5",
        "#08519c", "#08306b")
col_pal = colorRampPalette(pal)(ngrid)
col_val = seq(-r, r, length.out = ngrid)
lo_ind = findInterval(lo, col_val)
hi_ind = findInterval(hi, col_val)
colors = col_pal[lo_ind:hi_ind]

ggplot(gdat, aes(x = x, y = y, fill = z)) + geom_tile() +
  scale_x_continuous(xlab, expand = c(0, 0)) + scale_y_reverse(ylab,
                                                               breaks = 1:ncol(v), expand = c(0, 0)) + scale_fill_gradientn(legend_title,
                                                                                                                            colors = colors, limits = c(0, 1)) + guides(fill = guide_colorbar(barheight = bar_height)) +
  theme_bw(base_size = font_size) + theme(aspect.ratio = asp) + labs(title = TeX("$\\theta=0.6$")) +
  scale_y_continuous(name = "Coordinate", breaks=seq(0,p,10)) +
  scale_x_continuous(name = "Hz", breaks = seq(0,.5,by=.05)*(n), labels=function(x) x*64/n)

```

### Top Left Panel of Figure 2

```r
## Population
Localized_Est <- selector(f_evec11,f_xx1,n/2,(2/5)*512,p)
evecs <- t(Mod(Localized_Est[[1]]))
v = as.matrix(evecs)
lo = min(v)
hi = max(v)
r = max(abs(c(lo, hi)))
gdat = data.frame(x = as.integer(row(v)), y = as.integer(col(v)),
                  z = as.numeric(v))
ngrid = 1001
pal = c("#67000d", "#a50f15", "#cb181d", "#ef3b2c", "#fb6a4a",
        "#fc9272", "#fcbba1", "#fee0d2", "#ffffff", "#deebf7",
        "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5",
        "#08519c", "#08306b")
col_pal = colorRampPalette(pal)(ngrid)
col_val = seq(-r, r, length.out = ngrid)
lo_ind = findInterval(lo, col_val)
hi_ind = findInterval(hi, col_val)
colors = col_pal[lo_ind:hi_ind]

ggplot(gdat, aes(x = x, y = y, fill = z)) + geom_tile() +
  scale_x_continuous(xlab, expand = c(0, 0)) + scale_y_reverse(ylab,
                                                               breaks = 1:ncol(v), expand = c(0, 0)) + scale_fill_gradientn(legend_title,
                                                                                                                            colors = colors, limits = c(0, 1)) + guides(fill = guide_colorbar(barheight = bar_height)) +
  theme_bw(base_size = font_size) + theme(aspect.ratio = asp) + labs(title = TeX("Population")) +
  scale_y_continuous(name = "Coordinate", breaks=seq(0,p,10)) +
  scale_x_continuous(name = "Hz", breaks = seq(0,.5,by=.05)*(n), labels=function(x) x*64/n)

```

### Top Right Panel of Figure 2

```r
## Classic
Localized_Est <- selector(f_MT_evec11,f_xx1,n/2,512,p)
evecs <- t(Mod(Localized_Est[[1]]))
v = as.matrix(evecs)
lo = min(v)
hi = max(v)
r = max(abs(c(lo, hi)))
gdat = data.frame(x = as.integer(row(v)), y = as.integer(col(v)),
                  z = as.numeric(v))
ngrid = 1001
pal = c("#67000d", "#a50f15", "#cb181d", "#ef3b2c", "#fb6a4a",
        "#fc9272", "#fcbba1", "#fee0d2", "#ffffff", "#deebf7",
        "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5",
        "#08519c", "#08306b")
col_pal = colorRampPalette(pal)(ngrid)
col_val = seq(-r, r, length.out = ngrid)
lo_ind = findInterval(lo, col_val)
hi_ind = findInterval(hi, col_val)
colors = col_pal[lo_ind:hi_ind]

ggplot(gdat, aes(x = x, y = y, fill = z)) + geom_tile() +
  scale_x_continuous(xlab, expand = c(0, 0)) + scale_y_reverse(ylab,
                                                               breaks = 1:ncol(v), expand = c(0, 0)) + scale_fill_gradientn(legend_title,
                                                                                                                            colors = colors, limits = c(0, 1)) + guides(fill = guide_colorbar(barheight = bar_height)) +
  theme_bw(base_size = font_size) + theme(aspect.ratio = asp) + labs(title = TeX("Classic")) +
  scale_y_continuous(name = "Coordinate", breaks=seq(0,p,10)) +
  scale_x_continuous(name = "Hz", breaks = seq(0,.5,by=.05)*(n), labels=function(x) x*64/n)


```


# Data Analysis

## HC

You can use the following code to reproduce the top panels of Figure 5 of the main manuscript.

### Read in the Data

We first read in the data and estimate the spectral density matrices.

```r
## Data

X <- HC
p <- ncol(X)
n <- nrow(X)


## Multitaper estimate of the spectral density matrix
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
plot(Re(F_tp1[1,1,])*(n), type = "l", ylab = " ", main = "Estimated spectral density", ylim=c(1,2000))
for(a in 2:p){
  lines(Re(F_tp1[a,a,])*n, col= a)
}
f_xx1 <- F_tp1*n
rm(U)
rm(X_tp)
rm(F_tp_list)
#gc()

```

### Principal Subspace Estimation

Next we apply the LSPCA algorithm.

```r
## Localized and sparse PCA
LSDPCA_ADMM_SOAP_Ex2 <- LSPCA(n,p, f_xx1, lambda = 0.5 * sqrt(log(p) / n), d=2, lr = 0.02, maxiter = 60,
                                               control = list(fan_maxinc = 10, verbose = 0), eta=52, s=8, n_iter = 20, nu=0.6)

```

### Plots

#### First PC Loadings

The top left panel of Figure 5 of the main manuscript can be reproduced by the following code.

```r
## Plot

xlab <- "Hz"
ylab <- "Coordinate"
legend_title = ""
asp = .5 #0.2
bar_height = 10/2
font_size = 20/2



Localized_Est <- selector(LSDPCA_ADMM_SOAP_Ex2[[1]],f_xx1,n/2,52,p)
evecs <- t(Mod(Localized_Est[[1]]))
v = as.matrix(evecs)
lo = min(v)
hi = max(v)
r = max(abs(c(lo, hi)))
gdat = data.frame(x = as.integer(row(v)), y = as.integer(col(v)),
                  z = as.numeric(v))
ngrid = 1001
pal = c("#67000d", "#a50f15", "#cb181d", "#ef3b2c", "#fb6a4a",
        "#fc9272", "#fcbba1", "#fee0d2", "#ffffff", "#deebf7",
        "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5",
        "#08519c", "#08306b")
col_pal = colorRampPalette(pal)(ngrid)
col_val = seq(-r, r, length.out = ngrid)
lo_ind = findInterval(lo, col_val)
hi_ind = findInterval(hi, col_val)
colors = col_pal[lo_ind:hi_ind]

ggplot(gdat, aes(x = x, y = y, fill = z)) + geom_tile() +
  scale_x_continuous(xlab, expand = c(0, 0)) + scale_y_reverse(ylab,
                                                               breaks = 1:ncol(v), expand = c(0, 0)) + scale_fill_gradientn(legend_title,
                                                                                                                            colors = colors) + guides(fill = guide_colorbar(barheight = bar_height)) +
  theme_bw(base_size = font_size) + theme(aspect.ratio = asp) + labs(title = "HC - Modulus of the first PC loadings") +
  scale_y_continuous(name = "Coordinate", breaks=seq(0,p,10)) +
  scale_x_continuous(name = "Hz", breaks = seq(0,.5,by=.05)*(n), limits = c(0,320), labels=function(x) x*64/n)+
  #geom_vline(xintercept = 8*2, linetype="dotted", color = "red", size=.25)+
  geom_vline(xintercept = 64*2, linetype="dotted", color = "red", size=.25)+
  geom_vline(xintercept = 112*2, linetype="dotted", color = "red", size=.25)

```


#### Second PC Loadings

The top right panel of Figure 5 of the main manuscript can be reproduced by the following code.

```r
freq_s <- Localized_Est[[2]]
freq_selector <- matrix(rep(freq_s, p), nrow = p, byrow = TRUE)
evec2_Re <- freq_selector*Mod(LSDPCA_ADMM_SOAP_Ex2[[2]])
evecs <- t(evec2_Re)
v = as.matrix(evecs)
lo = min(v)
hi = max(v)
r = max(abs(c(lo, hi)))
gdat = data.frame(x = as.integer(row(v)), y = as.integer(col(v)),
                  z = as.numeric(v))
ngrid = 1001
pal = c("#67000d", "#a50f15", "#cb181d", "#ef3b2c", "#fb6a4a",
        "#fc9272", "#fcbba1", "#fee0d2", "#ffffff", "#deebf7",
        "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5",
        "#08519c", "#08306b")
col_pal = colorRampPalette(pal)(ngrid)
col_val = seq(-r, r, length.out = ngrid)
lo_ind = findInterval(lo, col_val)
hi_ind = findInterval(hi, col_val)
colors = col_pal[lo_ind:hi_ind]

ggplot(gdat, aes(x = x, y = y, fill = z)) + geom_tile() +
  scale_x_continuous(xlab, expand = c(0, 0)) + scale_y_reverse(ylab,
                                                               breaks = 1:ncol(v), expand = c(0, 0)) + scale_fill_gradientn(legend_title,
                                                                                                                            colors = colors) + guides(fill = guide_colorbar(barheight = bar_height)) +
  theme_bw(base_size = font_size) + theme(aspect.ratio = asp) + labs(title = "HC - Modulus of the second PC loadings") +
  scale_y_continuous(name = "Coordinate", breaks=seq(0,p,10)) +
  scale_x_continuous(name = "Hz", breaks = seq(0,.5,by=.05)*(n), limits = c(0,320), labels=function(x) x*64/n)+
  #geom_vline(xintercept = 8*2, linetype="dotted", color = "red", size=.25)+
  geom_vline(xintercept = 64*2, linetype="dotted", color = "red", size=.25)+
  geom_vline(xintercept = 112*2, linetype="dotted", color = "red", size=.25)

```

## FEP

You can use the following code to reproduce the bottom panels of Figure 5 of the main manuscript.

### Read in the Data

You can use the following code to reproduce bottom left panel of Figure 5 of the main manuscript.

We first read in the data and estimate the spectral density matrices.

```r
## Data

X <- FEP
p <- ncol(X)
n <- nrow(X)


## Multitaper estimate of the spectral density matrix
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
plot(Re(F_tp1[1,1,])*(n), type = "l", ylab = " ", main = "Estimated spectral density", ylim=c(1,400))
for(a in 2:p){
  lines(Re(F_tp1[a,a,])*n, col= a)
}
f_xx1 <- F_tp1*n
rm(U)
rm(X_tp)
rm(F_tp_list)
#gc()

```

### Principal Subspace Estimation

Next we apply the LSPCA algorithm.

```r
## Localized and sparse PCA
LSDPCA_ADMM_SOAP_Ex2 <- LSPCA(n,p, f_xx1, lambda = 0.5 * sqrt(log(p) / n), d=2, lr = 0.02, maxiter = 60,
                                               control = list(fan_maxinc = 10, verbose = 0), eta=41, s=8, n_iter = 20, nu=0.2)

```

### Plots

Finally, the bottom left panel of Figure 5 of the main manuscript can be reproduced by the following code.

```r
## Plot

xlab <- "Hz"
ylab <- "Coordinate"
legend_title = ""
asp = .5 #0.2
bar_height = 10/2
font_size = 20/2



Localized_Est <- selector(LSDPCA_ADMM_SOAP_Ex2[[1]],f_xx1,n/2,41,p)
evecs <- t(Mod(Localized_Est[[1]]))
v = as.matrix(evecs)
lo = min(v)
hi = max(v)
r = max(abs(c(lo, hi)))
gdat = data.frame(x = as.integer(row(v)), y = as.integer(col(v)),
                  z = as.numeric(v))
ngrid = 1001
pal = c("#67000d", "#a50f15", "#cb181d", "#ef3b2c", "#fb6a4a",
        "#fc9272", "#fcbba1", "#fee0d2", "#ffffff", "#deebf7",
        "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5",
        "#08519c", "#08306b")
col_pal = colorRampPalette(pal)(ngrid)
col_val = seq(-r, r, length.out = ngrid)
lo_ind = findInterval(lo, col_val)
hi_ind = findInterval(hi, col_val)
colors = col_pal[lo_ind:hi_ind]

ggplot(gdat, aes(x = x, y = y, fill = z)) + geom_tile() +
  scale_x_continuous(xlab, expand = c(0, 0)) + scale_y_reverse(ylab,
                                                               breaks = 1:ncol(v), expand = c(0, 0)) + scale_fill_gradientn(legend_title,
                                                                                                                            colors = colors) + guides(fill = guide_colorbar(barheight = bar_height)) +
  theme_bw(base_size = font_size) + theme(aspect.ratio = asp) + labs(title = "FEP - Modulus of the first PC loadings") +
  scale_y_continuous(name = "Coordinate", breaks=seq(0,p,10)) +
  scale_x_continuous(name = "Hz", breaks = seq(0,.5,by=.05)*(n), limits = c(0,320), labels=function(x) x*64/n)+
  #geom_vline(xintercept = 8*2, linetype="dotted", color = "red", size=.25)+
  geom_vline(xintercept = 64*2, linetype="dotted", color = "red", size=.25)+
  geom_vline(xintercept = 112*2, linetype="dotted", color = "red", size=.25)

```

#### Second PC Loadings

The bottom right panel of Figure 5 of the main manuscript can be reproduced by the following code.

```r
freq_s <- Localized_Est[[2]]
freq_selector <- matrix(rep(freq_s, p), nrow = p, byrow = TRUE)
evec2_Re <- freq_selector*Mod(LSDPCA_ADMM_SOAP_Ex2[[2]])
evecs <- t(evec2_Re)
v = as.matrix(evecs)
lo = min(v)
hi = max(v)
r = max(abs(c(lo, hi)))
gdat = data.frame(x = as.integer(row(v)), y = as.integer(col(v)),
                  z = as.numeric(v))
ngrid = 1001
pal = c("#67000d", "#a50f15", "#cb181d", "#ef3b2c", "#fb6a4a",
        "#fc9272", "#fcbba1", "#fee0d2", "#ffffff", "#deebf7",
        "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5",
        "#08519c", "#08306b")
col_pal = colorRampPalette(pal)(ngrid)
col_val = seq(-r, r, length.out = ngrid)
lo_ind = findInterval(lo, col_val)
hi_ind = findInterval(hi, col_val)
colors = col_pal[lo_ind:hi_ind]

ggplot(gdat, aes(x = x, y = y, fill = z)) + geom_tile() +
  scale_x_continuous(xlab, expand = c(0, 0)) + scale_y_reverse(ylab,
                                                               breaks = 1:ncol(v), expand = c(0, 0)) + scale_fill_gradientn(legend_title,
                                                                                                                            colors = colors) + guides(fill = guide_colorbar(barheight = bar_height)) +
  theme_bw(base_size = font_size) + theme(aspect.ratio = asp) + labs(title = "FEP - Modulus of the second PC loadings") +
  scale_y_continuous(name = "Coordinate", breaks=seq(0,p,10)) +
  scale_x_continuous(name = "Hz", breaks = seq(0,.5,by=.05)*(n), limits = c(0,320), labels=function(x) x*64/n)+
  #geom_vline(xintercept = 8*2, linetype="dotted", color = "red", size=.25)+
  geom_vline(xintercept = 64*2, linetype="dotted", color = "red", size=.25)+
  geom_vline(xintercept = 112*2, linetype="dotted", color = "red", size=.25)

```
