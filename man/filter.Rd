% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.R
\name{filter}
\alias{filter}
\title{Filter.}
\usage{
filter(data, threshold = 1)
}
\arguments{
\item{data}{list consisting of count data and sample table.}

\item{threshold}{minimum mean-log2-CPM to keep the data.}
}
\value{
input list containing filtered data.
}
\description{
Filter low expressing counts by mean logarithm counts-per-million.
}
