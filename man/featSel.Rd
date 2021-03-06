% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hypeRimg.R
\name{featSel}
\alias{featSel}
\title{Feature selection for final characterization}
\usage{
featSel(obj, k = 2)
}
\arguments{
\item{obj}{hypeRimg object}

\item{k}{parameter k from mixtools::normalmixEM()}
}
\description{
Utilizes normalmixE; implemented in the
mixtools package to test for the presence of > 1 population
and identifies minima. Cells are then assigned to high/low
groups.
}
