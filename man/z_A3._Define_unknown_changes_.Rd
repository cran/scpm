\name{A3. Define unknown changes}
\alias{A3. Define unknown changes}
\alias{cp}
\title{
Changes in the pattern of response
}
\description{
Define unknown changes in the pattern of response to be estimated.
}
\usage{
cp(x, psi, data = NULL, groups = NULL, contrasts = NULL,
        only.UV = FALSE)
}
\arguments{
	\item{x}{
		numeric vector. Covariate over which range define unknown change-points.
	}
	\item{psi}{
		numeric vector. Starting values for the change-points.
	}
	\item{data}{
		data frame. Where to search for the covariate?
	}
	\item{groups}{
		not used at the moment. To be implemented.
	}
	\item{contrasts}{
		character. A contrast method for factor covariates in groups. Default to \samp{contr.treatment}.
	}
	\item{only.UV}{
		logical. Not required.
	}

}
\author{
Mario A. Martinez Araya, \email{r@marioma.me}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~ 

