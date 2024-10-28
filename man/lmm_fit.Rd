\name{lmm_fit}
\alias{lmm_fit}
\title{Fit LMM harmonization and obtain model coefficients}
\description{Fit linear mixed-effects models (LMM) to test magnetic resonance imaging (MRI) data from different sites. Briefly, this function calls \code{\link[nlme]{lme}} (from the \code{nlme} package) massively to fit LMM for each ROI/voxel/vertex with the variables in \code{cov} as fixed-effects factors and \code{site} as a random intercept. Importantly, it also addresses specific brain imaging details (e.g., prescaling fMRI measures if needed or accounting for varying collinearity due to differing brain coverage; see "Neural correlates of human fear conditioning and sources of variability: A mega-analysis and normative modeling study of fMRI data from 2,199 individuals").}
\usage{
lmm_fit(dat, site, cov = NULL, n.min = 10, impute_missing_cov = FALSE,
        prescaling = FALSE, lin.hyp = NULL, verbose = TRUE)
}
\arguments{
  \item{dat}{matrix or data.frame with the MRI data (e.g., ROIs, voxels, vertexes, etc.).}
  \item{site}{factor specifying the site of each individual.}
  \item{cov}{matrix or data.frame with the fixed-effects  covariates.}
  \item{n.min}{(optional) number specifying the minimum size of a site to be analyzed.}
  \item{impute_missing_cov}{(optional) logical, whether to impute missing covariates.}
  \item{prescaling}{(optional) logical, whether to prescale the sites' data before conducting ComBat. See 'Details'.}
  \item{lin.hyp}{A hypothesis vector or matrix giving linear combinations of coefficients by rows as described in \code{\link[car]{linearHypothesis}} (from the \code{car} package).}
  \item{verbose}{(optional) logical, whether to print some messages during execution.}
}
\details{In several situations, the use of ComBat is problematic, such as when all participants of some sites have missing data for some ROIs or voxels. In these cases, a linear mixed effects model (LMM) may be preferable. On another note, setting \code{prescaling} equal to \code{TRUE} may be especially beneficial when the sites use different scales, as may be easily the case in fMRI mega-analyses.}
\value{A list of parameters plus the results of the LMM.}
\references{
  Pinheiro JC, Bates DM (2000).
  \emph{Mixed-Effects Models in S and S-PLUS}.
  Springer, New York, doi:10.1007/b98882.
  
  Neural correlates of human fear conditioning and sources of variability: A mega-analysis and normative modeling study of fMRI data from 2,199 individuals, to be submitted.
}
\author{Joaquim Radua}
\seealso{\code{\link{combat_fit}}, \code{\link{prescale_fit}}, and \code{\link[car]{linearHypothesis}}}
\examples{
raw_mri = combat_example[,6:19]
site = factor(combat_example$site)
mod = as.matrix(combat_example[,c("disorder", "age", "sex")])

# Estimate the effects of disorder with simple linear models
# WITHOUT harmonizing MRI data across sites
Out_raw = t(apply(raw_mri, 2, function (y) {
  m = summary(lm(y ~ mod, data = combat_example))
  coef(m)[2,c(1,3,4)]
}))

# Estimate the effects of disorder with simple linear models
# AFTER ComBat harmonizing MRI data across sites
combat = combat_fit(raw_mri, site, mod)
harmonized_mri = combat_apply(combat, raw_mri, site, mod)$dat.combat
Out_combat = t(apply(harmonized_mri, 2, function (y) {
  m = summary(lm(y ~ mod, data = combat_example))
  coef(m)[2,c(1,3,4)]
}))

# Estimate the effects of disorder with LMM harmonizing MRI data across sites
lmm = lmm_fit(raw_mri, site, mod)
Out_lmm = data.frame(
  b = lmm$b_map[[2]],
  t = lmm$t_map[[2]],
  p = 2 * pnorm(-abs(lmm$z_map[[2]]))
)
rownames(Out_lmm) = colnames(raw_mri)

# Results without harmonizing, with combat_fit and with lmm_fit:
cat("\nRaw results:\n")
print(round(Out_raw, 3))
cat("\nComBat results:\n")
print(round(Out_combat, 3))
cat("\nLMM results:\n")
print(round(Out_lmm, 3))

# Correlations between the three methods:
cat("\nCorrelation in coefficients:\n")
print(round(cor(cbind(Out_raw[,1], Out_combat[,1], Out_lmm[,1])), 3))
cat("\nCorrelation in p-values:\n")
print(round(cor(cbind(Out_raw[,3], Out_combat[,3], Out_lmm[,3])), 3))
}
\keyword{enigma}
\keyword{harmonization}
\keyword{mri}
