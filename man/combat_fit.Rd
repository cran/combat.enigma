\name{combat_fit}
\alias{combat_fit}
\alias{combat_apply}
\title{Fit and apply ComBat harmonization}
\description{Fit and apply ComBat to harmonize magnetic resonance imaging (MRI) data from different sites. Briefly, ComBat is a batch adjustment method that removes additive and multiplicative differences between sites due to the use of different scanning devices. \code{combat_fit} fits the ComBat model, while \code{combat_apply} applies the harmonization using the model. As detailed below, the original \code{combat} function from the \code{sva} package was first modified by Jean-Philippe Fortin for the harmonization of MRI data and then modified by Joaquim Radua to create separate functions for fitting and applying the harmonization and allow missing values and constant rows for its use within the Enhancing Neuro Imaging Genetics through Meta-Analysis (ENIGMA) Consortium. Finding the effects of the site in one sample and then removing them from another sample is interesting in several scenarios, such as case-control studies (where the effects of the site are better found in controls) or in machine-learning studies (where the effects of the site must be found in the training sample).}
\usage{
combat_fit(dat, site, cov = NULL, n.min = 10, impute_missing_cov = FALSE,
           prescaling = FALSE, impute_empty_sites = FALSE, eb = TRUE, verbose = TRUE)
combat_apply(combat_parameters, dat, site, cov = NULL, verbose = TRUE)
}
\arguments{
  \item{dat}{matrix or data.frame with the MRI data (e.g., ROIs, voxels, vertexes, etc.).}
  \item{site}{factor specifying the site of each individual.}
  \item{cov}{matrix or data.frame with the covariates.}
  \item{combat_parameters}{a list of combat parameters returned by \code{combat_fit}.}
  \item{n.min}{(optional) number specifying the minimum size of a site to be analyzed.}
  \item{impute_missing_cov}{(optional) logical, whether to impute missing covariates.}
  \item{prescaling}{(optional) logical, whether to prescale the sites' data before conducting ComBat. See 'Details'.}
  \item{impute_empty_sites}{(optional) logical, whether to impute data when all participants of some sites have missing data for some ROIs or voxels. See 'Details'.}
  \item{eb}{(optional) logical, whether to use empirical Bayes.}
  \item{verbose}{(optional) logical, whether to print some messages during execution.}
}
\details{The original ComBat function in the \code{sva} package was first modified by Jean-Philippe Fortin for the harmonization of MRI data and modified again by Joaquim Radua to create separate functions for fitting and applying the harmonization, allowing missings and constant rows and minor changes in the arguments of the functions to facilitate their use. The current version is thus based on the one described in "Increased power by harmonizing structural MRI site differences with the ComBat batch adjustment method in ENIGMA" below, but please also acknowledge the previous versions.
In situations where the use of ComBat is problematic, such as when all participants of some sites have missing data for some ROIs or voxels, a linear mixed effects model (LMM) may be preferable; please see \code{\link{lmm_fit}}. Finally, setting \code{prescaling} equal to \code{TRUE} may be especially beneficial when the sites use different scales, as it may be easily the case in fMRI mega-analyses (see "Neural correlates of human fear conditioning and sources of variability: A mega-analysis and normative modeling study of fMRI data from 2,199 individuals").}
\value{\code{combat_fit} returns a list of ComBat parameters for \code{combat_apply}; \code{combat_apply} returns the list of parameters plus the harmonized data (item \code{dat.combat})}
\references{
  Fortin JP, Cullen N, Sheline YI, Taylor WD, Aselcioglu I, Cook PA, Adams P, Cooper C, Fava M, McGrath PJ, McInnis M, Phillips ML, Trivedi MH, Weissman MM, Shinohara RT. (2017)
  Harmonization of cortical thickness measurements across scanners and sites.
  \emph{Neuroimage}, \bold{167}, 104--120, doi:10.1016/j.neuroimage.2017.11.024.
  
  Radua J, Vieta E, Shinohara R, Kochunov P, Quide Y, Green MJ, Weickert CS, Weickert T, Bruggemann J, Kircher T, Nenadic I, Cairns MJ, Seal M, Schall U, Henskens F, Fullerton JM, Mowry B, Pantelis C, Lenroot R, Cropley V, Loughland C, Scott R, Wolf D, Satterthwaite TD, Tan Y, Sim K, Piras F, Spalletta G, Banaj N, Pomarol-Clotet E, Solanes A, Albajes-Eizagirre A, Canales-Rodriguez EJ, Sarro S, Di Giorgio A, Bertolino A, Stablein M, Oertel V, Knochel C, Borgwardt S, du Plessis S, Yun JY, Kwon JS, Dannlowski U, Hahn T, Grotegerd D, Alloza C, Arango C, Janssen J, Diaz-Caneja C, Jiang W, Calhoun V, Ehrlich S, Yang K, Cascella NG, Takayanagi Y, Sawa A, Tomyshev A, Lebedeva I, Kaleda V, Kirschner M, Hoschl C, Tomecek D, Skoch A, van Amelsvoort T, Bakker G, James A, Preda A, Weideman A, Stein DJ, Howells F, Uhlmann A, Temmingh H, Lopez-Jaramillo C, Diaz-Zuluaga A, Fortea L, Martinez-Heras E, Solana E, Llufriu S, Jahanshad N, Thompson P, Turner J, van Erp T; ENIGMA Consortium collaborators. (2020)
  Increased power by harmonizing structural MRI site differences with the ComBat batch adjustment method in ENIGMA.
  \emph{Neuroimage}, \bold{218}, 116956, doi:10.1016/j.neuroimage.2020.116956.

  Neural correlates of human fear conditioning and sources of variability: A mega-analysis and normative modeling study of fMRI data from 2,199 individuals, submitted.
}
\author{Joaquim Radua, based on the previous work of others (see Details)}
\seealso{\code{\link{prescale_fit}} and \code{\link{lmm_fit}}}
\examples{
raw_mri = combat_example[,6:19]
site = factor(combat_example$site)

# Fit and apply ComBat to harmonize mri data across sites
mod = as.matrix(combat_example[,c("disorder", "age", "sex")])
combat = combat_fit(raw_mri, site, mod)
harmonized_mri = combat_apply(combat, raw_mri, site, mod)$dat.combat

# Run ANOVAs to check the effects of site
Out = NULL
for (j in 1:ncol(raw_mri)) {
  raw_mri_anova = anova(lm(raw_mri[,j] ~ mod), lm(raw_mri[,j] ~ mod + site))
  harmonized_mri_anova = anova(lm(harmonized_mri[,j] ~ mod), lm(harmonized_mri[,j] ~ mod + site))
  Out = rbind(Out, data.frame(
    roi = colnames(raw_mri)[j],
    raw_mri.F = round(raw_mri_anova$F[2], 1),
    raw_mri.P = round(raw_mri_anova$`Pr(>F)`[2], 3),
    harmonized_mri.F = round(harmonized_mri_anova$F[2], 1),
    harmonized_mri.P = round(harmonized_mri_anova$`Pr(>F)`[2], 3)
  ))
}
Out[,-1] = apply(Out[,-1], 2, function (x) {ifelse(x == 0, "<0.001", x)})
Out
}
\keyword{enigma}
\keyword{harmonization}
\keyword{mri}
