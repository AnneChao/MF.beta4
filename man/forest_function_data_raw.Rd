\name{forest_function_data_raw}
\alias{forest_function_data_raw}
\title{Ecosystem function data for six European forests}

\description{
In addition to plot information, this dataset includes raw values of 26 ecosystem functions collected from 209 plots (each with 30 m × 30 m) in six European countries, representing six major European forest types: boreal forest (Finland, 28 plots); hemi-boreal (Poland, 43 plots); temperate deciduous (Germany, 38 plots); mountainous deciduous (Romania, 28 plots); thermophilous deciduous (Italy, 36 plots); and Mediterranean mixed (Spain, 36 plots). See Table 1 of Ratcliffe et al. (2017a) for a description of the 26 functions. Each plot is designated as an ecosystem in assessing BEF relationships. See Ratcliffe et al. (2017b) and Scherer-Lorenzen et al. (2023) for the original dataset. For each missing value of functions in the original dataset, the mean of the given function within the country was imputed. An additional column “country” for each plot is added (as the last column) because function normalization and relevant analyses will be performed within each country. Thus, the dataset provided with the package is slightly different from the original one. 
}

\usage{
data("forest_function_data_raw")
}

\format{
a data.frame with 209 plots (rows) and 32 columns, in addition to row and column names; the first 5 columns show the relevant plot information, followed by 26 raw ecosystem functions (columns 6 to 31). The last column shows the corresponding country for each plot.
}


\references{

Ratcliffe, S. Wirth, C., Jucker, T. van der Plas, F., Scherer-Lorenzen, M. Verheyen, K. et al. (2017a). Biodiversity and ecosystem functioning relations in European forests depend on environmental context. Ecology Letters, 20, 1414–1426. \cr

Ratcliffe, S. Wirth, C., Jucker, T., van der Plas, F., Scherer-Lorenzen, M., Verheyen, K. et al. (2017b). Data for Biodiversity and ecosystem functioning relations in European forests depend on environmental con- text. https://doi.org/10.6084/m9.figshare.5368846.v1 \cr

Scherer-Lorenzen, M. et al. (2023). The functional significance of tree species diversity in European forests - the FunDivEUROPE dataset [Dataset]. Dryad. https://doi.org/10.5061/dryad.9ghx3ffpz

}

\keyword{datasets}
