\name{shoreacong}
\alias{shoreacong}
\docType{data}
\title{Tree pattern of the species Shorea congestiflora from Sri Lanka}
\description{The data record the locations of a tropical tree pattern of the species Shorea congestiflora from Sri Lanka in a 500 x 500 meter sampling region and their diameters at breast height (dbh).

The pattern was statistically analysed in Wiegand et al. (2007) by means of pair correlation functions for several size and age classes. In Wiegand and Moloney (2014), mark correlation function and mark variogram were estimated and briefly interpreted.
}
\usage{data(shoreacong)}
\format{
  A list with the following variables.
  \describe{
    \item{\code{x}}{a numeric vector containing the x-coordinates of the tree locations (in meters)}
    \item{\code{y}}{a numeric vector containing the y-coordinates of the tree locations (in meters)}
    \item{\code{dbh}}{a numeric vector containing the diameters at breast height (dbh) of the trees (in centimeters)}
    \item{\code{n}}{number of trees (= 850)}
    \item{\code{xrange}}{a numeric vector containing minimum and maximum x-coordinate of the rectangular observation window}
    \item{\code{yrange}}{a numeric vector containing minimum and maximum y-coordinate of the rectangular observation window}
  }
}
\source{
Prof T. Wiegand
}
\references{
Wiegand, T., Gunatilleke, C. V. S., Gunatilleke, I. A. U. N., Okuda, T., 2007. Analyzing the spatial structure of a Sri Lankan tree species with multiple scales of clustering. Ecology 88, 3088-3102.

Wiegand, T., Moloney, K.A., 2014. Handbook of spatial point pattern analysis, Chapman and Hall, Boca Raton.
}
\examples{
data(shoreacong)
}
\keyword{datasets}
