\name{panel.doublerug}
\alias{panel.doublerug}
\title{Lattice panel function for two rug plots}
\description{This panel function is defined to plot two rugs, one on top of the other in a multi-panel lattice display.}
\usage{
 panel.doublerug(xa = NULL, xb = NULL, 
	regular = TRUE, start = if (regular) 0 else 0.97, 
    end = if (regular) 0.03 else 1, x.units = rep("npc", 2), 
    lty = 1, lwd = 1) 
}
\arguments{
  \item{xa,xb}{Numeric vectors to be plotted.}
  \item{regular}{Logical flag indicating whether rug is to be drawn on the usual side (bottom/left) as opposed to the other side (top/right).}
  \item{start, end}{Start and end points for the rug ticks on the y-axis.}
  \item{x.units}{Character vectors, replicated to be of length two. Specifies the (grid) units associated with start and end above. x.units are for the rug on the x-axis and y-axis respectively (and thus are associated with start and end values on the y and x scales respectively). See \code{panel.rug} for more details.}
  \item{lty, lwd}{Line type and width arguments (see \code{par} for more details).}
}
\author{Dave Armstrong (UW-Milwaukee, Department of Political Science)}

