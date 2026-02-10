\name{selector}
\alias{selector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Frequency Slector
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
It selects the eta frequency components with the largest leading eigenvalue.
}
\usage{
selector(fxx_grad_evec, f_xx, len_freq, eta, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fxx_grad_evec}{
%%     ~~Describe \code{fxx_grad_evec} here~~
The estimated principal component loadings over the fundamental frequencies.
}
  \item{f_xx}{
%%     ~~Describe \code{f_xx} here~~
An estimate of the spectral density matrix.
}
  \item{len_freq}{
%%     ~~Describe \code{len_freq} here~~
len_freq = n/2 is the number of fundamental frequency components.
}
  \item{eta}{
%%     ~~Describe \code{eta} here~~
Number of frequency compoennts to be selected.
}
  \item{p}{
%%     ~~Describe \code{p} here~~
Dimension
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{LSDPCA}{The localized and sparse eigenvector series}
\item{best_sol}{0-1 valued vector representing the selected frequencies}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Jamshid Namdari
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

p <- ncol(D)
n <- nrow(D)

## Estimation
Ex1 <- LSPCA.f(n,p,f_D, d=1, eta=(2/5)*512+10, s=5, n_iter = 20, theta=0.8)

## Frequency Selection
LSPCA_EX1 <- selector(Ex1[[1]], f_D, n/2,(2/5)*512+10,p)

## Plot
plot(abs(LSPCA_EX1[[1]][1,]), main = "1st coordinate - Mod 1st sample evec over freq - eta = 500", pch=18, ylab="", ylim = c(0,1), type="l")

## We can compare it with the true eigenvectors
## We first construct the population Spectral density matrix, f_xx(.).
## We then compute the leading eigenvector of f_xx.

#############################################################################
## Below is the code to construct the population spectral density matrices
## Population Spectral Density Matrices

library(astsa)
library(mvtnorm)

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


## Leading eigenvector of the population spectral density matrices
f_evec11 <- matrix(0, nrow=p, ncol = length(omega))
for(ell in 1:length(omega)){
  f_evec_ell <- eigen(f_xx[,,ell])$vectors[,1]
  f_evec11[,ell] <- f_evec_ell
  #gc()
}

## Plot
plot(abs(LSPCA_EX1[[1]][1,]), main = "1st coordinate - Mod 1st sample evec over freq - eta = 500", pch=18, ylab="", ylim = c(0,1), type="l")
points(Mod(f_evec11[1,]), pch=18, ylab="", col = "green", type="l")
legend("topright", inset=.002, legend=c("Estimated", "Population"),
       col=c("black", "green"), lty=c(1,1))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
