\name{sssFit-class}
\Rdversion{1.1}
\docType{class}
\alias{sssFit-class}
\alias{sssFit}
\title{
Class \code{"sssFit"}}
\description{
Output object from \code{scp} that can be used with methods for post-processing.
}
\section{Objects from the Class}{
Objects of this class are created by calls to the command \code{scp} (see \link{scp}). It is also possible to define an empty object of this class by calls of the form \code{sssFit(name)} however for further use this is subject to validity of the object.
}
\section{Slots}{
	\describe{
		\item{data}{an object of \code{sss} class containing the input data.}
		\item{zV}{numeric vector. Response variables measured at the locations given by \code{data@coords}.}
		\item{XL}{a named list with elements and covariates from \code{linear} command.}
		\item{XC}{a named list with elements and covariates from \code{cp} command.}
		\item{XF}{not implemented.}
		\item{XS}{a named list with elements from \code{s2D} command.}
		\item{fit}{a named list with different estimated parameters and summaries from the estimated model.}
		\item{call}{call to the fitted model.}
	}
}
\author{
Mario A. Martinez Araya, \email{r@marioma.me}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~ 
\usage{
sssFit(...)
}
\arguments{
	\item{\dots}{
		Slots elements to be included into the \code{sssFit} object. Allowed slots names are \code{data}, \code{zV}, \code{XL}, \code{XC}, \code{XF}, \code{XS}, \code{fit}, and \code{call}.
	}

}

