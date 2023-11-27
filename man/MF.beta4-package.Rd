\name{MF.beta4-package}
\alias{MF.beta4-package}
\title{Measuring ecosystem multifunctionality and assessing BEF relationships}

\description{

\code{MF.beta4} is an R package for measuring ecosystem multifunctionality and assessing BEF relationships. The measures are illustrated by using ecosystem function and biodiversity data collected in a total of 209 plots in six European countries (the FunDivEUROPE dataset). All data are available from the Dryad repository; see Ratcliffe et al. (2017a, b) and Scherer-Lorenzen et al. (2023) for data details. The software was originally developed for the Beta4 project (Müller et al. 2022) on the effect of enhancing the beta diversity between forest patches on ecosystem multifunctionality and forest resilience across spatial scales.    \cr


Based on a framework of Hill-Chao numbers of orders q = 0, 1 and 2, \code{MF.beta4} features the following multifunctionality measures for a single and multiple ecosystems; see Chao et al. (2023) for pertinent methodology and decomposition theory. \cr

(1)	Multifunctionality measures in a single ecosystem: \code{MF.beta4} computes a class of weighted multifunctionality measures for given function weights. Multifunctionality measures that correct for strong correlations between ecosystem functions to avoid redundancy are also provided. \cr

(2)	Multifunctionality measures in multiple ecosystems: for given function weights, \code{MF.beta4} computes the gamma multifunctionality of pooled ecosystems, the within-ecosystem component (alpha multifunctionality) and the among-ecosystem component (beta multifunctionality). The correlation between functions can also be corrected for.   \cr

Based on biodiversity and function data from ecosystems, this package also provides graphics for assessing biodiversity-ecosystem functioning (BEF) relationships across scales.   \cr


This package includes four main functions:   \cr

1. \code{function_normalization} transforms ecosystem functions data to values between 0 and 1. \cr 

2. \code{MF1_single} computes a class of weighted multifuctionality measures in a single ecosystem for given individual function weights separately for two cases: (i) correlations between functions are not corrected for, and (ii) correlations between functions are corrected for. \cr

3. \code{MF2_multiple} computes alpha, beta and gamma multifuctionality measures in multiple ecosystems for given function weights separately for two cases (i) correlations between functions are not corrected for, and (ii) correlations between functions are corrected for.\cr

4. \code{MFggplot} provides the graphical BEF relationships based on the output from the function \code{MF1_single} or \code{MF2_multiple}.

}

\author{
Anne Chao, Chun-Yu Liu, K. H. Hu  \cr
Maintainer: Anne Chao <chao@stat.nthu.edu.tw> 
}

\references{

Chao, A., Chiu, C. H., Hu, K. H., van der Plas, F., Cadotte, M. W., Mitesser, O., et al. (2023). Hill-Chao numbers in multifunctionality allows decomposing gamma multifunctionality into alpha and beta components. To appear in Ecology Letters.  

Müller, J., Mitesser, O. Cadotte, M. W., van der Plas, F., Mori, A, Ammer, C., Eisenhauer N. (2023). Enhancing the structural diversity between forest patches – a concept and real-world experiment to study biodiversity and multifunctionality across spatial scales. Global Change Biology, 29, 1437–1450.

Ratcliffe, S. Wirth, C., Jucker, T. van der Plas, F., Scherer-Lorenzen, M. Verheyen, K. et al. (2017a). Biodiversity and ecosystem functioning relations in European forests depend on environmental context. Ecology Letters, 20, 1414–1426.

Ratcliffe, S. Wirth, C., Jucker, T., van der Plas, F., Scherer-Lorenzen, M., Verheyen, K. et al. (2017b). Data for Biodiversity and ecosystem functioning relations in European forests depend on environmental context. https://doi.org/10.6084/m9.figshare.5368846.v1

Scherer-Lorenzen, M. et al. (2023). The functional significance of tree species diversity in European forests - the FunDivEUROPE dataset [Dataset]. Dryad. https://doi.org/10.5061/dryad.9ghx3ffpz

}
