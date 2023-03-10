% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_lat_dim.R
\name{gen_lat_dim}
\alias{gen_lat_dim}
\title{Generate MDS Space}
\usage{
gen_lat_dim(ndims = 2, n = 16, min = -1, max = 1)
}
\arguments{
\item{ndims}{number of dimensions}

\item{n}{number of stimuli}

\item{min}{minimum value the coordinates can have (defualt -1)}

\item{max}{maximum value the coordinates can have (default 1)}
}
\value{
matrix with \code{ndims} columns and \code{n} rows
}
\description{
Generate MDS Space
}
\examples{
gen_lat_dim()
}
\author{
David Izydorczyk
}
