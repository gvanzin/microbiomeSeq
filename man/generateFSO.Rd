% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_fso.R
\name{generateFSO}
\alias{generateFSO}
\title{Fuzzy Set Ordination}
\usage{
generateFSO(physeq, grouping_column, method, indices = NULL,
  filename = NULL, type = 1, step_across = F)
}
\arguments{
\item{grouping_column}{(Required). Character string specifying name of a categorical variable that is preffered for grouping the information.
information.}

\item{type}{(optional). An integer (1 or 2), for fso or mfso.}

\item{step_across}{(optional). logical variable setting it to TRUE is for step-across correction which might improve the ordination}

\item{physeq(Required).}{A \code{phyloseq} object containing merged information of abundance,
taxonomic assignment, sample data including the measured variables and categorical information
of the samples, and / or phylogenetic tree if available.}

\item{method.}{An integer specifying method for computing similarity idices.}
}
\description{
This function uses fuzzy set ordination to test effects of pertubation in
environmental variables to community structure.
It returns  an ordination plot which is annotated with a correlation between the original values of
the variable and the fuzzy values together with corresponding significance label.
}
\references{
\url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015

\url{http://www.nku.edu/~boycer/fso/}
}
\author{
Alfred Ssekagiri \email{assekagiri@gmail.com}, Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
}
