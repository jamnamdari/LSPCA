\name{HC}
\alias{HC}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
EEG recording of a Healthy Control subject
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
The dataframe contains 64 chanell EEG recording of a Healthy Control subject collected as part of an NIH funded study (R01MH113827)). Data were collected while the participant sat in a chair and relaxed with their eyes open. The data are sampled at 64 Hz for one minutes.
}
\format{
 A data frame with 1028 rows and 64 columns. Columns represent chanells and rows represent time points.
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

## Localized and sparse PCA

LSPCA_HC <- LSPCA(HC, d=2, eta=52, s=8, n_iter = 20, theta=0.6)


## Plot
library(ggplot2)
library(latex2exp)


xlab <- "Hz"
ylab <- "Coordinate"
legend_title = ""
asp = .5 #0.2
bar_height = 10/2
font_size = 20/2
p <- ncol(FEP)
n <- nrow(FEP)


freq_s <- LSPCA_HC$selected_freqs
freq_selector <- matrix(rep(freq_s, p), nrow = p, byrow = TRUE)
evecs <- t(freq_selector*Mod(LSPCA_HC$Evec_series_1))
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


}
\keyword{datasets}
