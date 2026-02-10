\name{Freq_par_selector}
\alias{Freq_par_selector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Frequency parameter selector
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Selects the localization (frequency) parameter
}
\usage{
Freq_par_selector(X, f_xx1 = NULL, fxx_evec_s, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
nxp data matrix
}
  \item{f_xx1}{
%%     ~~Describe \code{f_xx1} here~~
estimated spectral density matrix. If set to NULL, it will estimated it by the multitapper estimator using 10 sine tapper functions.
}
  \item{fxx_evec_s}{
%%     ~~Describe \code{fxx_evec_s} here~~
The sparse eigenvector series estimate
}
  \item{d}{
%%     ~~Describe \code{d} here~~
Dinesion of the principal subspace
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
\item{freq_par_AIC}{frequency parameter selected by AIC}
\item{freq_par_AICc}{frequency parameter selected by AICc}
\item{freq_par_BIC}{frequency parameter selected by BIC}
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


Ex1 <- LSPCA(D, d=1, eta=500, s=5, n_iter = 20, theta=0.8, ntp=20)

eta_hat <- Freq_par_selector(D,f_D,Ex1[[1]], d=1)

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
