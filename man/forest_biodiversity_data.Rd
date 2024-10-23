\name{forest_biodiversity_data}
\alias{forest_biodiversity_data}
\title{Biodiversity data for six European forests}

\description{
In addition to row and column names, this dataset consists of four columns: the “plotID” column indicates the names of plots, the “species” column includes species names, the column “abundance” (basal area as a proxy of species abundance), and the “country” column includes the corresponding stratifying variable. Because missing values of “basal area” in the original dataset were imputed by the mean of the same species within the country, and basal areas were combined for two species (“Betula pendula” and “Betula pubescens”), the dataset provided with the package is slightly different from the original dataset provided in Scherer-Lorenzen et al. (2023). 
}

\usage{
data("forest_biodiversity_data")
}

\format{
a data.frame with 481 rows (the total number of combinations of plot and tree species in 209 plots) and 4 columns (plotID, species name, the corresponding basal area as species abundance, and country).
}

\references{

Scherer-Lorenzen, M. et al. (2023). The functional significance of tree species diversity in European forests - the FunDivEUROPE dataset [Dataset]. Dryad. https://doi.org/10.5061/dryad.9ghx3ffpz

}


\keyword{datasets}
