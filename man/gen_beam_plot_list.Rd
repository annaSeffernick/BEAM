% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_beam_plot_list.R
\name{gen_beam_plot_list}
\alias{gen_beam_plot_list}
\title{Generate BEAM Plot List}
\usage{
gen_beam_plot_list(
  beam.result,
  beam.specs,
  beam.feat.pvals,
  number.pairs = 1,
  set.id,
  feat.id = NULL,
  title.size = 10,
  pair.order = "both",
  endpt.order = NULL
)
}
\arguments{
\item{beam.result}{Result of prep.beam.data}

\item{beam.specs}{A data.frame of strings with columns name, mtx, mdl, plot}

\item{beam.feat.pvals}{List of feature-level p-values from compute_feature_pvalues}

\item{number.pairs}{Numeric; number of features to display in clinical plots, ordered by significance}

\item{set.id}{A character with set name; must be in beam.result$beam.data$set.data$set.id}

\item{feat.id}{Default NULL; a character with feature name; must be in beam.result$beam.data$set.data$row.id}

\item{title.size}{A numeric. Specify the size of individual plot titles. Default is 10.}

\item{pair.order}{One of c("both", "omic", "endpoint"). Default is "both." Specify how to choose feature-endpoint plots to include. If "both", find the best (based on q, p, effect size) feature-omic pair for each type of omic and each endpoint separately. If "omic", within each omic, find the best feature-endpoint pair and then plot this feature with all endpoints. If "endpoint", need to specify endpt.order as the name of chosen endpoint. Then, within each omic, find the feature with best association with the selected endpoint, and plot this feature for all endpoints.}

\item{endpt.order}{Default NULL. If pair.order="endpoint", specify character with endpoint name (from beam.specs$name, after the period).}
}
\value{
A list of plots for the specified set and/or feature.
}
\description{
Internal function: generate a list of clinical feature plots.
}
\examples{
test.beam.data <- prep_beam_data(main.data=clinf, mtx.data=omicdat,
                                 mtx.anns=omicann, set.data=setdat,
                                 set.anns=NULL, n.boot=10, seed=123)
specs <- prep_beam_specs(beam.data=test.beam.data, endpts=c("MRD29", "EFS"), firth=TRUE)
test.beam.stats <- compute_beam_stats(beam.data=test.beam.data, beam.specs=specs)
test.feat.pvals <- compute_feature_pvalues(beam.stats=test.beam.stats)
plot.specs <- prep_beam_plot(beam.data=test.beam.data, beam.specs=specs)
plot.list <- gen_beam_plot_list(beam.result=test.beam.stats, beam.specs=plot.specs,
                                beam.feat.pvals=test.feat.pvals, number.pairs=1,
                                set.id="ENSG00000099810", feat.id=NULL, title.size=11,
                                pair.order="omic", endpt.order=NULL)
}
