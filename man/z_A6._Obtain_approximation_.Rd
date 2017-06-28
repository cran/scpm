\name{A6. Obtain approximation}
\alias{A6. Obtain approximation}
\alias{scpApproximate}
\alias{scpApproximate,sssFit-method}
\title{
Linear approximation to a spline.
}
\description{
Approximation to a spatial semiparametric model based on a bivariate spline.
}
\usage{
\S4method{scpApproximate}{sssFit}(object, tol)
}
\arguments{
	\item{object}{
		an object of class \code{sssFit} from command \code{scp}.
	}
	\item{tol}{
		numeric. Numeric tolerance to use for some inversion of matrices. Default to \code{.Machine$double.neg.eps*1.0e-10}.
	}

}
\details{
\code{scpApproximate} compute an approximation to the spatial semiparametric model obtained from \code{scp}. This command update the fitted values and fitted spline in the input object of class \code{sssFit}. Then we can use the command \code{plot} for plotting the approximated semiparametric model.  
}
\value{
This command return an object of class \code{sssFit}.
}
\author{
Mario Martinez Araya, \email{mma@mariomartinezaraya.com}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~ 

