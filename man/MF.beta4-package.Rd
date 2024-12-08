\name{MF.beta4-package}
\alias{MF.beta4-package}
\title{Measuring ecosystem multifunctionality and assessing BEF relationships}

\description{

\code{MF.beta4} is an R package for measuring ecosystem multifunctionality and assessing biodiversity–ecosystem function (BEF) relationships. The measures are illustrated using ecosystem function and biodiversity data collected from a total of 209 plots across six European countries (the FunDivEUROPE dataset). All data are available in the Dryad repository; see Ratcliffe et al. (2017a, b) and Scherer-Lorenzen et al. (2023) for details. The software was originally developed for the Beta4 project (Müller et al. 2022), which studied the effect of enhancing beta diversity between forest patches on ecosystem multifunctionality and forest resilience across spatial scales.    \cr


Based on a framework of Hill-Chao numbers of orders q = 0, 1 and 2, \code{MF.beta4} features the following multifunctionality measures for a single and multiple ecosystems; see Chao et al. (2024) for
pertinent methodology and decomposition theory. \cr

(1)	Multifunctionality measures in a single ecosystem: \code{MF.beta4} computes a class of weighted multifunctionality measures for given function weights. Multifunctionality measures that correct for strong correlations between ecosystem functions, in order to avoid redundancy, are also provided. When biodiversity data are available, \code{MF.beta4} also provides graphics for assessing biodiversity-ecosystem functioning (BEF) relationships between within-ecosystem multifunctionality and species diversity for orders q = 0, 1, and 2. \cr

(2)	Multifunctionality measures in multiple ecosystems: for given function weights, \code{MF.beta4} computes the gamma multifunctionality of pooled ecosystems, the within-ecosystem component (alpha multifunctionality) and the among-ecosystem component (beta multifunctionality). The correlation between functions can also be corrected for. 

When biodiversity data are available, \code{MF.beta4} also provides graphics to assess biodiversity-ecosystem functioning (BEF) relationships between gamma/alpha/beta multifunctionality and species diversity for orders q = 0, 1, and 2, comparing all pairs of ecosystems/plots or multiple ecosystems/plots.  \cr


This package includes four main functions:   \cr

1. \code{function_normalization} transforms ecosystem functions data to values between 0 and 1. \cr 

2. \code{MF1_single} computes a class of weighted multifuctionality measures in a single ecosystem for given individual function weights separately for two cases: (i) correlations between functions are not corrected for, and (ii) correlations between functions are corrected for. \cr

3. \code{MF2_multiple} computes alpha, beta and gamma multifuctionality measures of orders q = 0, 1 and 2 for given function weights for all pairs of ecosystems (if \code{by_pair = TRUE}) or for multiple ecosystems (if \code{by_pair = FALSE}) separately for two cases (i) correlations between functions are not corrected for, and (ii) correlations between functions are corrected for.\cr

4. \code{MFggplot} provides the graphical BEF relationships based on the output from the function \code{MF1_single} or \code{MF2_multiple}. \cr\cr



NOTE: Latest updates as of Dec. 8, 2024: In earlier versions, multifunctionality decomposition (alpha, beta and gamma) was performed only for pairs of plots/ecosystems. In the updated version, we have added a logical argument \code{by_pair} in the main function \code{MF2_multiple} to specify whether multifunctionality decomposition will be performed for all pairs of ecosystems or not.  If \code{by_pair =TRUE}, alpha/beta/gamma multifunctionality will be computed for all pairs of ecosystems/plots in the input data; if \code{by_pair = FALSE}, alpha/beta/gamma multifunctionality will be computed for K plots (i.e., K can be greater than two) when data for K plots are provided in the input data. Default is \code{by_pair =TRUE}.

}

\author{
Anne Chao, Chun-Yu Liu, K. H. Hu  \cr
Maintainer: Anne Chao <chao@stat.nthu.edu.tw> 
}

\references{

Chao, A., Chiu, C. H., Hu, K. H., van der Plas, F., Cadotte, M. W., Mitesser, O., et al. (2024). Hill-Chao numbers in multifunctionality allows decomposing gamma multifunctionality into alpha and beta components. Ecology Letters. 27, e14336. Available from: https://doi.org/10.1111/ele.14336 

Müller, J., Mitesser, O. Cadotte, M. W., van der Plas, F., Mori, A, Ammer, C., Eisenhauer N. (2022). Enhancing the structural diversity between forest patches - a concept and real-world experiment to study biodiversity and multifunctionality across spatial scales. Global Change Biology, 29, 1437-1450.

Ratcliffe, S. Wirth, C., Jucker, T. van der Plas, F., Scherer-Lorenzen, M. Verheyen, K. et al. (2017a). Biodiversity and ecosystem functioning relations in European forests depend on environmental context. Ecology Letters, 20, 1414-1426.

Ratcliffe, S. Wirth, C., Jucker, T., van der Plas, F., Scherer-Lorenzen, M., Verheyen, K. et al. (2017b). Data for Biodiversity and ecosystem functioning relations in European forests depend on environmental context. https://doi.org/10.6084/m9.figshare.5368846.v1

Scherer-Lorenzen, M. et al. (2023). The functional significance of tree species diversity in European forests - the FunDivEUROPE dataset [Dataset]. Dryad. https://doi.org/10.5061/dryad.9ghx3ffpz

}
