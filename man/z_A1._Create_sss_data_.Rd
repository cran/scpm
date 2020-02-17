\name{A1. Create sss data}
\alias{A1. Create sss data}
\alias{as.sss}
\alias{create.sss}
\alias{is.sss}
\alias{sss2df}
\title{
Convert an object to the class \samp{sss} for spatial smoothing splines
}
\description{
Create a matrix or data.frame to a valid dataset of class \samp{sss} for spatial smoothing splines. Those dataset can be used later by functions \link{s2D} for tensor product (natural) cubic splines or p-splines, and \link{scp} for estimating spatial smoothing splines models.
}
\usage{
as.sss(X, coords, coords.col, data.col, ...)
create.sss(coords, data, ...)
is.sss(x)
sss2df(x)
}
\arguments{
	\item{X}{
		a matrix or data-frame. Every row must correspond to a point location in a two-dimensional space (coordinates). Coordinates columns can be included in \code{X} or defined separately using the argument \code{coords}. Some columns can also correspond to variables measured at the different point locations.
	}
	\item{coords}{
		two-columns numeric matrix of coordinates (optional).
	}
	\item{data}{
		a data-frame containing the variables measured at the locations given by \code{coords}.
	}
	\item{coords.col}{
		numeric vector. The number of columns in \code{X} that contains the coordinates.
	}
	\item{data.col}{
		numeric vector. The number of columns in \code{X} that contains variables measured at the points locations.
	}
	\item{\dots}{
		slots elements to create a new \code{sss} dataset. Required slots are \code{data}, \code{coords}, \code{grid}, \code{knots}, \code{W}, \code{contract} (to be discarded in the future), and \code{regular}. See \emph{Value} for an explanation about each slot requirements.
	}
	\item{x}{
		an object to check validity as member of class \code{sss}.
	}

}
\value{
	\describe{
		\item{data}{a data-frame containing the variables measured at the locations given by \code{coords}.}
		\item{coords}{a matrix containing the two columns of observed coordinates for the data.}
		\item{grid}{a grid matrix containing the two columns of coordinates.}
		\item{knots}{a named list with the design points (knots) in every coordinate. Equivalent to a \code{grid.list} object.}
		\item{W}{a spatial incidence matrix. If \code{contract=TRUE} it is \eqn{W_{ij}}, otherwise \eqn{W_{ji}}.}
		\item{contract}{logical. The same value as the argument \code{contract}.}
		\item{regular}{logical. If the coordinates are observed at regular points it is \code{TRUE}, \code{FALSE} otherwise (missing coordinates in any direction).}
	}
}
\author{
Mario A. Martinez Araya, \email{r@marioma.me}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~ 

