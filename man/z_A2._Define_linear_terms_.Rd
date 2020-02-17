\name{A2. Define linear terms}
\alias{A2. Define linear terms}
\alias{linear}
\title{
Linear components of the mean of the model
}
\description{
Define parametric components of the mean as linear terms.
}
\usage{
linear(formula, data = NULL, contrasts = NULL, intercept = FALSE)
}
\arguments{
	\item{formula}{
		formula expression. A formula expression as described in \link{formula}.
	}
	\item{data}{
		data frame. Where to search for the covariates?
	}
	\item{contrasts}{
		character. A contrast method for factor covariates. Default to \samp{contr.treatment}.
	}
	\item{intercept}{
		logical. \code{TRUE} to include an intercept term, \code{FALSE} otherwise (default).
	}

}
\author{
Mario A. Martinez Araya, \email{r@marioma.me}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~ 

