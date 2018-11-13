mmsc <- function(X, rmin, rmax, type=c("symmetric", "plain", "lower", "upper"), mono=FALSE, contour.col=c(4, 2), log=FALSE, mlim=range(X[,3]), mlab=c(expression(m[i]), expression(m[j])), main="Mark-mark scatterplot", maxs="r", pch=20, cex=1, col=1, contour.lwd=1, gos=0.1, extreme.ticks=TRUE, n.gray=12) {
  if(!is.matrix(X))
    stop("Argument 'X' must be a matrix.")
  if(NCOL(X)<3)
    stop("Matrix 'X' has less than three columns.")
  X <- na.omit(X)
  if(NROW(X)<2)
    stop("Matrix 'X' has less than two rows after removing NAs.")
  n <- length(X[,1])
  if(missing(rmin)) {
    rmin <- 0.0
    warning("Unspecified argument 'rmin' has been set to 0.")
  }
  if(missing(rmax)) {
    R <- unlist(sapply(1:(n-1), function(i){rij <- sqrt((X[i,1]-X[(i+1):n,1])^2+(X[i,2]-X[(i+1):n,2])^2)}))
    rmax <- sort(R)[100]
    warning(paste("Unspecified argument 'rmax' has been set to ",rmax,", which corresponds to 100 mark pairs in the scatterplot.",sep=""))
    rm(R)
  }
  if((rmin<0.0) || (rmax<0.0) || !(rmin<rmax))
    stop("Argument 'rmin' must be nonnegative and less than argument 'rmax'.")
  if(log) log="xy"
  else log=""
  if(main=="") corr <- 0.0
  else corr <- 0.4
  oldpar <- par(no.readonly = TRUE)
  if(!mono) layout(matrix(c(1, 2), ncol = 2), widths = c(5, 1))
  else layout(1, widths = 5)
  par(mai=c(0.8, 0.8+corr, 0.1+corr, 0.1))
  plot(NA, NA, type="n", xlim=mlim, ylim=mlim, xlab=mlab[1], ylab=mlab[2], main=main, log=log, las=1, xaxs=maxs, yaxs=maxs)
  mlen <- 101
  if(length(contour.col)<2) contour.col <- c(4, 2)
  mm <- seq(mlim[1], mlim[2], len=mlen)
  kmm <- outer(mm, mm, function(m1,m2) m1*m2/mean(X[,3])^2)
  contour(mm, mm, kmm, add=TRUE, col=contour.col[1], levels=1, drawlabels=FALSE, lwd=contour.lwd)
  gam <- outer(mm, mm, function(m1,m2) 0.5*(m1-m2)^2)/var(X[,3])
  contour(mm, mm, gam, add=TRUE, col=contour.col[2], levels=1, drawlabels=FALSE, lwd=contour.lwd)
  if(log=="") lines(c(-1,mlim[2]*1.5),c(-1,mlim[2]*1.5))
  else lines(c(1e-9,mlim[2]*1.5),c(1e-9,mlim[2]*1.5))

  MMR<-unlist(sapply(1:(n-1), function(i){
    rij <- sqrt((X[i,1]-X[(i+1):n,1])^2+(X[i,2]-X[(i+1):n,2])^2)
    id <- ((i+1):n)[rij>rmin & rij<=rmax]
    if(length(id)>0) rbind(X[i,3], X[id,3], rij[rij>rmin & rij<=rmax])
  }))
  if(!is.null(MMR)) MMR <- matrix(MMR, ncol=3, byrow=TRUE)
  message(paste("Number of inter-point distances between rmin=",rmin," and rmax=",rmax," is ",NROW(MMR),".",sep=""))
  if(length(MMR)>0) {
    if(!mono) MMR <- MMR[order(MMR[,3],decreasing=TRUE),] # darker points (corresponding to larger distances) are plotted first
    if(!mono) col <- gray(1-(MMR[,3]+gos-rmin)/(rmax+gos-rmin))
  }
  type <- match.arg(type)
  switch(type,
         symmetric = for(i in 1:length(MMR[,3])) points(c(MMR[i,1],MMR[i,2]), c(MMR[i,2],MMR[i,1]), pch=pch, col=col[i], cex=cex),
         plain = points(MMR[,1], MMR[,2], pch=pch, col=col, cex=cex),
         lower = for(i in 1:length(MMR[,3])) points(max(MMR[i,1:2]), min(MMR[i,1:2]), pch=pch, col=col[i], cex=cex),
         upper = for(i in 1:length(MMR[,3])) points(min(MMR[i,1:2]), max(MMR[i,1:2]), pch=pch, col=col[i], cex=cex))
  rm(MMR)
  if(!mono) {
    par(mai=c(0.8, 0.2, 0.1+corr, 0.6))
    plot.new()
    plot.window(xlim=c(0, 1), ylim=c(rmin,rmax), xaxs="i", yaxs="i")
    levels <- seq(rmin, rmax, len=n.gray+1)
    rect(0, levels[-length(levels)], 1, levels[-1L], col=gray(1-(levels+gos-rmin)/(rmax+gos-rmin)))
    if(extreme.ticks) axis(4, las=1, at=c(rmin, rmax))
    else axis(4, las=1)
    axis(1, at=0.5, tick=FALSE, labels="r")
  }
  if(!mono) layout(1)
  par(oldpar)
}
