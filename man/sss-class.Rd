\name{sss-class}
\Rdversion{1.1}
\docType{class}
\alias{sss-class}
\alias{sss}
\title{
Class \code{"sss"}}
\description{
Create a dataset for spatial smoothing splines. Those dataset can be used later by functions \link{s2D} for tensor product (natural) cubic splines or p-splines, and \link{scp} for estimating spatial smoothing splines models.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{sss(...)}.
%%  ~~ describe objects here ~~ 

}
\section{Slots}{
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
Mario Martinez Araya, \email{mma@mariomartinezaraya.com}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~ 
\usage{
sss(...)
}
\arguments{
	\item{\dots}{
		Slots elements to be included into the \code{sss} dataset. Allowed slots names are \code{data}, \code{coords}, \code{grid}, \code{knots}, \code{W}, \code{contract} (to be discarded in the future), and \code{regular}.
	}

}

