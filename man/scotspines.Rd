\name{scotspines}
\alias{scotspines}
\docType{data}
\title{Tree pattern of Scots pines from Siberia (Russia)}
\description{The data record the locations of a tree pattern of Scots pines (Pinus sylvestris L.) from Siberia (Russia) in a 56 x 56 meter sampling region, their total heights and their diameters at breast height (dbh).

The trees of the plot are 90 years old. The data were collected and ecologically interpreted by Busykin et al. (1985). The data only contain the inner region of the originally 60 x 60 meter observation window in order to reduce spatial inhomogeneity.

The data have been analysed e. g. by Gavrikov et al. (1993),  Gavrikov and Stoyan (1995) and Genet et al. (2014).
}
\usage{data(scotspines)}
\format{
  A list with the following variables.
  \describe{
    \item{\code{x}}{a numeric vector containing the x-coordinates of the tree locations (in meters)}
    \item{\code{y}}{a numeric vector containing the y-coordinates of the tree locations (in meters)}
    \item{\code{height}}{a numeric vector containing the total heights of the trees (in meters)}
    \item{\code{dbh}}{a numeric vector containing the diameters at breast height (dbh) of the trees (in centimeters)}
    \item{\code{n}}{number of trees (=484)}
    \item{\code{xrange}}{a numeric vector containing minimum and maximum x-coordinate of the rectangular observation window}
    \item{\code{yrange}}{a numeric vector containing minimum and maximum y-coordinate of the rectangular observation window}
  }
}
\source{
Busykin et al. (1985)
}
\references{
Busykin, A.I., Gavrikov, V.L., Sekretenko, O.P., Hlebopros, R.G. (1985). Analysis of forest cenoses (in Russian), Nauka, Novosibirsk.

Gavrikov, V.L., Grabarnik, P. Ya., Stoyan, D. (1993). Trunk-top relations in a Siberian pine forest. Biometrical Journal 35, 487-498.

Gavrikov, V., Stoyan, D. (1995). The use of marked point processes in ecological and environmental forest studies. Environmental and Ecological Statistics 2, 331-344.

Genet, A., Grabarnik, P., Sekretenko, O., Pothier, D. (2014). Incorporating the mechanisms underlying inter-tree competition into a random point process model to improve spatial tree pattern analysis in forestry. Ecological Modelling 288, 143-154.
}
\examples{
data(scotspines)
}
\keyword{datasets}
