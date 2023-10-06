\name{combat_fit}
\alias{combat_fit}
\title{Fit ComBat harmonization}
\description{Fit ComBat to harmonize magnetic resonance imaging (MRI) data from different sites. Briefly, ComBat is a batch adjustment method that removes additive and multiplicative differences between sites due to the use of different scanning devices. \code{combat_fit} fits the ComBat model, while \code{\link{combat_apply}} applies the harmonization using the model. As detailed below, the original \code{combat} function from the \code{sva} package was first modified by Jean-Philippe Fortin for the harmonization of MRI data and then modified by Joaquim Radua to create separate functions for fitting and applying the harmonization and allow missing values and constant rows for its use within the Enhancing Neuro Imaging Genetics through Meta-Analysis (ENIGMA) Consortium.}
\usage{
combat_fit(dat, site, cov = NULL, eb = TRUE, verbose = TRUE)
}
\arguments{
  \item{dat}{matrix or data.frame with the MRI data (e.g., ROIs, voxels, vertexes, etc.).}
  \item{site}{factor specifying the site of each individual.}
  \item{cov}{matrix or data.frame with the covariates.}
  \item{eb}{(optional) logical, whether to use empirical Bayes.}
  \item{verbose}{(optional) logical, whether to print some messages during execution.}
}
\details{The original function in the \code{sva} package was first modified by Jean-Philippe Fortin for the harmonization of MRI data and modified again by Joaquim Radua to create separate functions for fitting and applying the harmonization, allowing missings and constant rows and minor changes in the arguments of the functions to facilitate their use. The current version is thus the one described in "Increased power by harmonizing structural MRI site differences with the ComBat batch adjustment method in ENIGMA" below, but please also acknowledge the previous versions.}
\value{A list of ComBat parameters for \code{\link{combat_apply}}}
\references{
  Fortin JP, Cullen N, Sheline YI, Taylor WD, Aselcioglu I, Cook PA, Adams P, Cooper C, Fava M, McGrath PJ, McInnis M, Phillips ML, Trivedi MH, Weissman MM, Shinohara RT. (2017)
  Harmonization of cortical thickness measurements across scanners and sites.
  \emph{Neuroimage}, \bold{167}, 104--120, doi:10.1016/j.neuroimage.2017.11.024.
  
  Radua J, Vieta E, Shinohara R, Kochunov P, Quide Y, Green MJ, Weickert CS, Weickert T, Bruggemann J, Kircher T, Nenadic I, Cairns MJ, Seal M, Schall U, Henskens F, Fullerton JM, Mowry B, Pantelis C, Lenroot R, Cropley V, Loughland C, Scott R, Wolf D, Satterthwaite TD, Tan Y, Sim K, Piras F, Spalletta G, Banaj N, Pomarol-Clotet E, Solanes A, Albajes-Eizagirre A, Canales-Rodriguez EJ, Sarro S, Di Giorgio A, Bertolino A, Stablein M, Oertel V, Knochel C, Borgwardt S, du Plessis S, Yun JY, Kwon JS, Dannlowski U, Hahn T, Grotegerd D, Alloza C, Arango C, Janssen J, Diaz-Caneja C, Jiang W, Calhoun V, Ehrlich S, Yang K, Cascella NG, Takayanagi Y, Sawa A, Tomyshev A, Lebedeva I, Kaleda V, Kirschner M, Hoschl C, Tomecek D, Skoch A, van Amelsvoort T, Bakker G, James A, Preda A, Weideman A, Stein DJ, Howells F, Uhlmann A, Temmingh H, Lopez-Jaramillo C, Diaz-Zuluaga A, Fortea L, Martinez-Heras E, Solana E, Llufriu S, Jahanshad N, Thompson P, Turner J, van Erp T; ENIGMA Consortium collaborators. (2020)
  Increased power by harmonizing structural MRI site differences with the ComBat batch adjustment method in ENIGMA.
  \emph{Neuroimage}, \bold{218}, 116956, doi:10.1016/j.neuroimage.2020.116956.
}
\author{Joaquim Radua, based on the previous work of others (see Details)}
\seealso{\code{\link{combat_apply}}}
\examples{
raw_mri = combat_example[,6:19]
site = factor(combat_example$site)

# Fit and apply combat to harmonize mri data across sites
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
