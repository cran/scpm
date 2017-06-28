\name{B1. Testing surface}
\alias{B1. Testing surface}
\alias{testSurface}
\alias{testSurface,sssFit-method}
\title{
Testing the surface model
}
\description{
Test the model for the surface of response. The null hypothesis is assumed as a linear model defined by the coordinates while the alternative hypothesis is assumed a bivariate spline (tensor product or thin-plate spline).
}
\usage{
\S4method{testSurface}{sssFit}(object, tol)
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
If we have defined a bivariate spline using \link{s2D} in the formula of \link{scp} then the model is an spatial semiparametric model based on splines (tensor products or thin-plate splines). In this case \code{testSurface} performs a test for the null hypothesis \eqn{H_0: g = X\beta} (linear model) against the alternative \eqn{H_1: g = X\beta + Zr} (spline model). When \eqn{g} is assumed as a thin-plate spline then this test is equivalent to test the null hypothesis  \eqn{H_0:} \emph{the pattern of response in the space is a plane} against the alternative \eqn{H_1:} \emph{the pattern of response in the space is a bivariate thin-plate spline}. In one dimension this test is equivalent to a test for linearity in the pattern of response.
}
\value{
Returns a table with the degrees of freedom, sum of squares and mean squares from different sources and the F test and its associated p-value.
}
\author{
Mario A. Martinez Araya <mma@mariomartinezaraya.com> [aut, cre]
}
%% ~Make other sections like Warning with \section{Warning }{....} ~ 

