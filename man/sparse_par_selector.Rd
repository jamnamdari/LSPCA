\name{sparse_par_selector}
\alias{sparse_par_selector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Sparsity paramter selector
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Selects the sparsity level by minimizing the Mahalanobis distance
}
\usage{
sparse_par_selector(X, eta, nu = 0, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
nxp data matrix
}
  \item{eta}{
%%     ~~Describe \code{eta} here~~
Localization paramter
}
  \item{nu}{
%%     ~~Describe \code{nu} here~~
Smoothing parameter
}
  \item{d}{
%%     ~~Describe \code{d} here~~
Dimension of the principal subspace
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
\item{sparse_par}{Selected sparsity level}
\item{TMD}{Vector of total Mahalanobis distances computed for s=1,...,20}
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
s_hat <- sparse_par_selector(D,eta_hat$freq_par_AIC,0.8, d=1)
plot(s_hat$TMD)
s_hat$sparse_par

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
