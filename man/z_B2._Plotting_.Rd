\name{B2. Plotting}
\alias{B2. Plotting}
\alias{plot}
\alias{plot,sssFit,missing-method}
\title{
Plot observed, fitted values or estimated spline from \code{scp} object.
}
\description{
Draw an image, perspective or levelplot of the observed, fitted values or estimated bivariate spline from the elements in a \code{scp} object.
}
\usage{
\S4method{plot}{sssFit,missing}(x, what, type, which, col.args, col.contour,
	    level.at, border, theta, phi, shade, ...)
}
\arguments{
	\item{x}{
		sssFit object from \link{scp}.
	}
	\item{what}{
		character. What to plot? One of \samp{obs} (for observed responses), \samp{fit} (for fitted values, the default) or \samp{g} (for the estimated bivariate spline).
	}
	\item{type}{
		character. Which type of plot? One of \samp{image} (the default), \samp{levelplot}, \samp{persp} or \samp{persp3d}.
	}
	\item{which}{
		character. Which color pattern? One of \samp{colorRampPalette} (default), \samp{colorRamp}, \samp{rainbow}, \samp{heat.colors}, \samp{terrain.colors}, \samp{topo.colors}, or \samp{cm.colors}.
	}
	\item{col.args}{
		named list. List with argument to pass to the color pattern function defined by \code{which}. See \link{colorRampPalette}, \link{colorRamp}, \link{rainbow}, \link{heat.colors}, \samp{terrain.colors}, \samp{topo.colors}, and \samp{cm.colors}.
	}
	\item{col.contour}{
		character. Only for \code{type="image"}. Color for the contours.
	}
	\item{level.at}{
		character or numeric vector. Only for \code{type="levelplot"}. Where to draw the levels at. If character, it is the name of the function to compute where to put the levels.
	}
	\item{border}{
		character. Color of the border.
	}
	\item{theta, phi, shade}{
		numeric. See \link{persp} or \link{persp3d}.
	}
	\item{\dots}{
		other arguments for \link{levelplot}, \link{image}, \link{persp} or \link{persp3d}.
	}

}
\author{
Mario A. Martinez Araya, \email{r@marioma.me}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~ 

