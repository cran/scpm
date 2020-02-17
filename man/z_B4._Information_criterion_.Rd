\name{B4. Information criterion}
\alias{AIC}
\alias{AIC,sssFit-method}
\alias{AICc}
\alias{AICc,sssFit-method}
\alias{AICm}
\alias{AICm,sssFit-method}
\alias{B4. Information criterion}
\alias{BIC}
\alias{BIC,sssFit-method}
\alias{BICc}
\alias{BICc,sssFit-method}
\alias{BICj}
\alias{BICj,sssFit-method}
\alias{GIC}
\alias{GIC,sssFit-method}
\alias{GICb}
\alias{GICb,sssFit-method}
\alias{GIChq}
\alias{GIChq,sssFit-method}
\alias{GICpn}
\alias{GICpn,sssFit-method}
\title{
Information criterion of the estimated model from \code{scp} object.
}
\description{
Return the information criterion of the estimated model from a \code{scp} object.
}
\usage{
\S4method{AIC}{sssFit}(object, k, only.criterion)
\S4method{BIC}{sssFit}(object, only.criterion)
\S4method{AICm}{sssFit}(object, k, only.criterion)
\S4method{AICc}{sssFit}(object, k, only.criterion)
\S4method{BICc}{sssFit}(object, only.criterion)
\S4method{BICj}{sssFit}(object, k, tol, only.criterion)
\S4method{GIC}{sssFit}(object, k, only.criterion)
\S4method{GIChq}{sssFit}(object, k, only.criterion)
\S4method{GICpn}{sssFit}(object, only.criterion)
\S4method{GICb}{sssFit}(object, only.criterion)
}
\arguments{
	\item{object}{
		sssFit object from \link{scp}.
	}
	\item{k}{
		numeric. Factor multiplying the number of parameters in each criterion. Default to \code{k=2}.
	}
	\item{tol}{
		numeric. Value for the tolerance in some computation of inverse matrices. By default is set to \code{.Machine$double.neg.eps}.
	}
	\item{only.criterion}{
		logical. If \code{TRUE} (the default) returns only the value of the criterion.
	}

}
\details{
The information criterion for a mixed model is defined as \deqn{IC = -2\ell + penalty} where \eqn{\ell} is the log-likelihood \eqn{\ell(\vartheta)} or conditional log-likelihood \eqn{\ell(\vartheta|r)} (see \link{scp}). The penalty is expressed as \eqn{k\times  a_0\times \omega_{\mu_*,V}} where \eqn{\omega_{\mu_*,V} = \omega_{\mu_*} + \omega_V} is the (effective) number of parameters in the mean and variance and \eqn{k} and \eqn{a_0} are factors that depend on the criterion used. Thus the information criterion can be written as \deqn{IC = -2\ell + k\times  a_0\times \omega_{\mu_*,V}.} Note that \eqn{\mu_*} depends on the criterion being used so it can be \eqn{\mu_* = \mu_m} or \eqn{\mu_* = \mu}. See \link{scp}.
}
\value{
If \code{only.criterion=TRUE} returns the value of the criterion. If \code{only.criterion=FALSE} returns a list with the following elements:
	\describe{
		\item{\bold{logLik}}{numeric. The log-likelihood or conditional log-likelihood (given \eqn{r}) of the model depending of the criterion used.}
		\item{\bold{criterion}}{numeric. The value of the information criterion.}
		\item{\bold{ka0}}{numeric. Factors \eqn{ka_0} multiplying the number of parameters. Depends on the criterion selected.}
		\item{\bold{numpar}}{numeric. The (effective) number of parameters. Depends on the criterion selected.}
		\item{\bold{penalty}}{numeric. The value of the penalty.}
	}
}
\references{
	\itemize{
		\item{\enc{Mueller}{Muller}, Samuel; Scealy, J. L. and Welsh, A. H. (2013) Model Selection in Linear Mixed Models. Statist. Sci. \bold{28}, no. 2, 135--167. \href{http://doi.org/10.1214/12-STS410}{doi:10.1214/12-STS410}. \url{http://projecteuclid.org/euclid.ss/1369147909}.}
	}
}
\author{
Mario A. Martinez Araya, \email{r@marioma.me}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~ 

