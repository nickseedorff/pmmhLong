% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{gen_data_mats}
\alias{gen_data_mats}
\title{Build data matrices for ingestion by gen_pmll_estim()}
\usage{
gen_data_mats(distance_mat, sigma, phi, mean_vec)
}
\arguments{
\item{distance_mat}{list of objects to make relevant data matrices}

\item{sigma}{scalar}

\item{phi}{scalar}

\item{mean_vec}{vector, linear predictor without the GP}
}
\value{
list of matrices
}
\description{
Build data matrices for ingestion by gen_pmll_estim()
}
