% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/motif_plot.R
\name{motif_plot}
\alias{motif_plot}
\title{Code for plotting the Motif based on a specific CDR3 length and V gene (see \link{netplot_ClusTCR2} for ).}
\usage{
motif_plot(
  ClusTCR,
  Clust_column_name = "Clust_size_order",
  Clust_selected = NULL
)
}
\arguments{
\item{ClusTCR}{Matrix file produce from \link{mcl_cluster}}

\item{Clust_column_name}{Name of clustering column from mcl_cluster file e.g. cluster}

\item{Clust_selected}{Select which cluster to display. Only one at a time.}
}
\value{
A ggplot object representing the motif.
}
\description{
Code for plotting the Motif based on a specific CDR3 length and V gene (see \link{netplot_ClusTCR2} for ).
}
\examples{
# Example usage of mcl_cluster function with a stored file
example_file <- read.csv(system.file("extdata", "my_data.csv",package = "ClusTCR2"))
# Perform clustering using mcl_cluster function
step1 <- ClusTCR(example_file,allele = FALSE)
# perform mcl
step2 <- mcl_cluster(step1)
# print the motif plot for the simple clustering
print(motif_plot(step2,Clust_selected = 1))
}
