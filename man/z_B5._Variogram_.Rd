\name{B5. Variogram}
\alias{B5. Variogram}
\alias{Variogram}
\alias{Variogram,sssFit-method}
\title{
Compute and plot the semi-variogram from \code{scp} object.
}
\description{
Compute and plot the semi-variogram of the semiparametric model from a \code{scp} object.
}
\usage{
\S4method{Variogram}{sssFit}(object, distance, plot, ...)
}
\arguments{
	\item{object}{
		sssFit object from \link{scp}.
	}
	\item{distance}{
		numeric vector. The distances at which to compute the semi-variogram. By default is set to \code{NULL}.
	}
	\item{plot}{
		logical. \code{plot=TRUE} (the default) produce the semivariogram plot. \code{plot=FALSE} returns the values of the semivariogram at \code{distance}. See \emph{value}.
	}
	\item{\dots}{
		other graphical parameters to pass.
	}

}
\author{
Mario Martinez Araya, \email{mma@mariomartinezaraya.com}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~ 

