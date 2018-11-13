\name{mmsc}
\alias{mmsc}
\title{Mark-Mark Scatterplot}
\description{
\code{mmsc} creates a mark-mark scatterplot of a planar marked point pattern, i.e. plots the marks of all point pairs with inter-point distances in a certain interval in a mark-mark coordinate system.
}
\usage{
mmsc(X, rmin, rmax, type=c("symmetric", "plain", "lower", "upper"), mono=FALSE, contour.col=c(4, 2), log=FALSE,
     mlim=range(X[,3]), mlab=c(expression(m[i]), expression(m[j])), main="Mark-mark scatterplot", maxs="r",
     pch=20, cex=1, col=1, contour.lwd=1, gos=0.1, extreme.ticks=TRUE, n.gray=12)}
\arguments{
  \item{X}{matrix representing the planar marked point pattern, i.e., the first two columns contain point coordinates and the third column contains the (numerical) marks associated with the points.}
  \item{rmin}{lower interval limit for inter-point distances.}
  \item{rmax}{upper interval limit for inter-point distances.}
  \item{type}{character string specifying the way in which the pairs \eqn{(m_i,m_j)}{(m[i],m[j])} of marks are arranged. With \code{"symmetric"} (the default) \eqn{(m_i,m_j)}{(m[i],m[j])} as well \eqn{(m_j,m_i)}{(m[j],m[i])} are plotted; with \code{"plain"} \eqn{(m_i,m_j)}{(m[i],m[j])} are plotted for \eqn{i<j} in the order given by \code{X}; with \code{"lower"} and \code{"upper"} \eqn{(m_i,m_j)}{(m[i],m[j])} are plotted such that \eqn{m_i>m_j)}{(m[i]>m[j])} and, respectively, \eqn{m_i<m_j)}{(m[i]<m[j])} is satisfied.}
  \item{mono}{logical: if \code{TRUE}, draw all mark-mark points with the same color \code{col}; if \code{FALSE} (the default), draw each mark-mark point with a gray tone between white or light-gray encoded by \code{gos} and black proportional to the inter-point distance between \code{rmin} and \code{rmax}.}
  \item{contour.col}{colors for two optional contour lines. See the 'Details' section. The default is \code{c(4, 2)}. The contour lines may be suppressed by choosing each \code{0} as color code.}
  \item{log}{logical: if \code{TRUE}, the axes are plotted with logarithmic scale. The default is \code{log = FALSE}.}
  \item{mlim}{limits for both axes. The default is \code{mlim = range(X[,3])}.}
  \item{mlab}{vector of two character strings for the names of the both axes.}
  \item{main}{plot title.}
  \item{maxs}{The style of axis interval calculation to be used for the both mark axes. Possible values are "r" (the default) and "i". Works like \code{xaxs} and \code{yaxs} in \code{\link[graphics]{par}}.}
  \item{pch}{plotting 'character' as in \code{\link[graphics]{points}}.}
  \item{cex}{character (or symbol) expansion. The default value is \code{cex = 1}.}
  \item{col}{color for the points drawn.}
  \item{contour.lwd}{line width for the contour lines drawn.}
  \item{gos}{lower limit for the range of possible gray tones. Must be not less than \code{0} and not larger than \code{1}. The default is \code{gos = 0.1}.}
  \item{extreme.ticks}{logical: if \code{TRUE} (the default), only ticks at \code{rmin} and \code{rmax} are displayed in the key on the right-hand side; if \code{FALSE}, ticks are chosen automatically.}
  \item{n.gray}{number of gray-level boxes displayed in the key on the right-hand side. The default is \code{n.gray = 12}.}

}
\details{The mark-mark scatterplot for a (planar) marked point pattern is obtained by plotting the marks of all point pairs with an inter-point distance \eqn{r} in an interval \eqn{[r_{\min},r_{\max}]}{[r_min,r_max]} in a mark-mark coordinate system.

There are two possible refinements.\cr
1. According to the inter-point distance of a point pair in the interval \eqn{[r_{\min},r_{\max}]}{[r_min,r_max]}, the points may be plotted in a corresponding gray tone between white/light-gray (= distance \eqn{r_{\min}}{r_min}) and black (= distance \eqn{r_{\max}}{r_max}).\cr
2. The mark-mark scatterplot may be complemented by contour lines. Currently the contour line corresponding to mark pairs \eqn{(m_i,m_j)}{(m[i],m[j])} which have a product \eqn{m_im_j}{m[i]*m[j]} equal to the square of the mean of marks as well as the contour line corresponding to mark pairs \eqn{(m_i,m_j)}{(m[i],m[j])} which have half squared differences \eqn{0.5(m_i-m_j)^2}{0.5*(m[i]-m[j])^2} equal to the variance of marks are implemented.
}
\author{Felix Ballani}
\examples{
## positions in a 100 x 100 window
x <- runif(200, 0, 100)
y <- runif(200, 0, 100)
## marks = distance to nearest neighbor
m <- c(sapply(1:200, function(i) min(sqrt((x[-i]-x[i])^2+(y[-i]-y[i])^2))))

## mark-mark scatterplot for inter-point distances from [0,6]

## with default values
mmsc(X=cbind(x,y,m), rmin=0.0, rmax=6.0)

## asymmetric
mmsc(X=cbind(x,y,m), rmin=0.0, rmax=6.0, type="lower")

## without supporting contour lines
mmsc(X=cbind(x,y,m), rmin=0.0, rmax=6.0, contour.col=c(0, 0))

## with more ticks in the key
mmsc(X=cbind(x,y,m), rmin=0.0, rmax=6.0, extreme.ticks=FALSE)

## without gray-tones
mmsc(X=cbind(x,y,m), rmin=0.0, rmax=6.0, mono=TRUE)
}

