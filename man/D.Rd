\name{D}
\alias{D}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
A Simulated Data
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
A simulated dataset used to produce Figure 2 in the main manuscript.
}
\usage{data("D")}
\format{
  The format is:
 num [1:1024, 1:64] -0.302 -0.182 -0.259 -0.514 -0.537 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:64] "Xt1" "Xt2" "Xt3" "Xt4" ...
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{

# Loading the data
data(D)

# Estimating the spectral density matrix of the underlying process
# Using multitaper estimator

library(waveslim)
library(astsa)

U <- sine.taper(n,20)
X_tp <- apply(U, MARGIN = 2, function(u) u*D, simplify = FALSE)
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

# Plot
omega <- seq(0,.5, length.out = n/2)
plot(omega, f_xx1[1,1,1:512], type="l", ylim = c(0,30), ylab="", xlab="Frequency")
for (k in 2:64) {
  points(omega, f_xx1[k,k,1:512], type="l", col=k)
}


}
\keyword{datasets}
