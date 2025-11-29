\name{LSPCA.f}
\alias{LSPCA.f}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Localized Sparse Principal Component Analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The `LSPCA.f` function estimates the localized and sparse principal components of a multivariate time series in the frequency domain. Users can pass the frequency transformation of the data, i.e. an estimate of the spectral density matrices, to the function.
}
\usage{
LSPCA.f(n, p, f_xx, d, eta, s, n_iter, theta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
Length of the time series
}
  \item{p}{
%%     ~~Describe \code{p} here~~
Dimension of the time series
}
  \item{f_xx}{
%%     ~~Describe \code{f_xx} here~~
An array of dimensions $p\times p\times n$ contining an estimate of the $p\times p$-dimensional spectral density matrices at frequencies $f=\frac{\pi\ell}{n}, \ell=1,\dots,n/2$.
}
  \item{d}{
%%     ~~Describe \code{d} here~~
Dimension of the principal subspaces to be estimated
}
  \item{eta}{
%%     ~~Describe \code{eta} here~~
Number of frequency components to be selected
}
  \item{s}{
%%     ~~Describe \code{s} here~~
Sparsity level of the principal subspaces
}
  \item{n_iter}{
%%     ~~Describe \code{n_iter} here~~
Number of iteration of the SOPA algorithem
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
Smoothing parameter
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
  \item{Evec_series_k}{The k-th eigenvector series, k=1,2,...}
  \item{LSDPCA}{The localized and sparse eigenvector series}
  \item{best_sol}{0-1 valued vector representing the selected frequencies}
  \item{f_deflated}{Deflated spectral denisty matrix (only returned when d = 1)}
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

  library(dplR)
  library(mvtnorm)
  library(dplR)
  library(waveslim)
  library(astsa)

  ##########################################################################################
  ## Example 1 : This example reproduces bottom right panel of Figure 2 in the manuscript
  ##########################################################################################


  ##################################################
  ## Estimate of 1-dimensional principal subspaces
  ##################################################

LSDPCA_Ex1 <- LSPCA.f(n,p,f_D, d=1, eta=(2/5)*512, s=5, n_iter = 20, theta=0.6)



## Plot
library(latex2exp)

xlab <- "Hz"
ylab <- "Coordinate"
legend_title = ""
asp = .5 #0.2
bar_height = 10/2
font_size = 20/2



freq_s <- LSDPCA_Ex1$selected_freqs
freq_selector <- matrix(rep(freq_s, p), nrow = p, byrow = TRUE)
evecs <- t(freq_selector*Mod(LSDPCA_Ex1$Evec_series_1))
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
  theme_bw(base_size = font_size) + theme(aspect.ratio = asp) + labs(title = TeX("$\\theta=0.6$")) +
  scale_y_continuous(name = "Coordinate", breaks=seq(0,p,10)) +
  scale_x_continuous(name = "frequency", breaks = seq(0,.5,by=.05)*(n), labels=function(x) x*64/n)


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
