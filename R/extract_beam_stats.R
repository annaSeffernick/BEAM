#' Extract beam stats for a specific set
#'
#' @param beam.stats A beam.stats object, which is a list with beam.stats (the association matrices), the beam.specs, and the beam.data
#' @param set.id A character of a set id name (an entry in in beam.data$set.data$set.id)
#'
#' @return A matrix with with estimated associations for each endpoint and each omic feature linked to the set
#' @export
#'
#' @examples
#' test.beam.data <- prep_beam_data(main.data=clinf, mtx.data=omicdat,
#'                                  mtx.anns=omicann, set.data=setdat,
#'                                  set.anns=NULL, n.boot=10, seed=123)
#' specs <- prep_beam_specs(beam.data=test.beam.data, endpts=c("MRD29", "EFS"), firth=TRUE)
#' test.beam.stats <- compute_beam_stats(beam.data=test.beam.data, beam.specs=specs)
#' test.stats <- extract_beam_stats(test.beam.stats, set.id="ENSG00000099810")
extract_beam_stats=function(beam.stats,set.id)
{
  set.mtch=which(beam.stats$beam.data$set.data$set.id==set.id) # extract rows corresponding to the set.id
  if (length(set.mtch)==0)
    stop(paste0("set.id ", set.id," not found.")) # stop if set.id not found
  mtch.set=beam.stats$beam.data$set.data[set.mtch,] # extract relevant set.data with set.id, mtx.id, row.id columns
  n.set.rows=length(mtch.set$set.id) # number of rows in set (number of omics types available?)
  beam.stat.mtx=NULL
  beam.stat.rownames=NULL
  for (i in 1:n.set.rows)
  {
    stat.mtch=which(substring(names(beam.stats$beam.stats),1,nchar(mtch.set$mtx.id[i]))==mtch.set$mtx.id[i]) # extract indices to match each omics type
    for (j in stat.mtch)
    {
      row.mtch=which(rownames(beam.stats$beam.stats[[j]])==mtch.set$row.id[i])
      beam.stat.mtx=rbind(beam.stat.mtx,beam.stats$beam.stats[[j]][row.mtch,]) # Extract submatrix
      beam.stat.rownames=c(beam.stat.rownames,
                           paste(rownames(beam.stats$beam.stats[[j]])[row.mtch],
                                 names(beam.stats$beam.stats[j])))
      #mtch.set$mtx.id[j])) # Extract rownames
    }
  }
  rownames(beam.stat.mtx)=beam.stat.rownames
  return(beam.stat.mtx)
}
