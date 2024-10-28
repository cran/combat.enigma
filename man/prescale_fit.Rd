\name{prescale_fit}
\alias{prescale_fit}
\alias{prescale_apply}
\title{Fit and apply brain imaging prescaling}
\description{Fit and apply a prescaling of magnetic resonance imaging (MRI) data from different sites, especially relevant for fMRI, where different devices can have varying units of measurement.}
\usage{
prescale_fit(dat, site, cov = NULL, n.min = 10, impute_missing_cov = FALSE,
             verbose = TRUE)
prescale_apply(combat_parameters, dat, site, cov = NULL, verbose = TRUE)
}
\arguments{
  \item{dat}{matrix or data.frame with the MRI data (e.g., ROIs, voxels, vertexes, etc.).}
  \item{site}{factor specifying the site of each individual.}
  \item{cov}{matrix or data.frame with the covariates.}
  \item{combat_parameters}{a list of combat parameters returned by \code{prescale_fit}.}
  \item{n.min}{(optional) number specifying the minimum size of a site to be analyzed.}
  \item{impute_missing_cov}{(optional) logical, whether to impute missing covariates.}
  \item{verbose}{(optional) logical, whether to print some messages during execution.}
}
\details{\code{prescale_fit} function finds a prescaling parameter for each site so that, after calling \code{prescale_apply}, the voxelwise-median standard deviation after removing the effects of covariates is 1 (in the training data, but it might not be the case in new data).}
\value{\code{prescale_fit} returns a list of parameters for \code{prescale_apply}; \code{prescale_apply} returns the list of parameters plus the prescaled data (item \code{dat.combat})}
\references{
  Neural correlates of human fear conditioning and sources of variability: A mega-analysis and normative modeling study of fMRI data from 2,199 individuals, to be submitted.
}
\author{Joaquim Radua}
\seealso{\code{\link{combat_fit}}, and \code{\link{lmm_fit}}}
\examples{
raw_mri = combat_example[,6:19]
site = factor(combat_example$site)

# Fit and apply prescale to prescale mri data across sites
mod = as.matrix(combat_example[,c("disorder", "age", "sex")])
prescaling = prescale_fit(raw_mri, site, mod)
prescaled_mri = prescale_apply(prescaling, raw_mri, site, mod)$dat.combat
}
\keyword{enigma}
\keyword{harmonization}
\keyword{mri}
