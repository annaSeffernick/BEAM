#' Compute feature level p-values from BEAM statistics
#'
#' @param beam.stats A beam.stats object, which is a list with beam.stats (the association matrices), the beam.specs, and the beam.data
#'
#' @return A list of feature level p-values, with each entry a data frame for a different omics/endpoint associaiton, with columns id, gene, beta, p, q
#' @importFrom stats p.adjust
#' @export
#'
#' @examples
#' test.beam.data <- prep_beam_data(main.data=clinf, mtx.data=omicdat,
#'                                  mtx.anns=omicann, set.data=setdat,
#'                                  set.anns=NULL, n.boot=10, seed=123)
#' specs <- prep_beam_specs(beam.data=test.beam.data, endpts=c("MRD29", "EFS"), firth=TRUE)
#' test.beam.stats <- compute_beam_stats(beam.data=test.beam.data, beam.specs=specs)
#' test.feat.pvals <- compute_feature_pvalues(beam.stats=test.beam.stats)
compute_feature_pvalues=function(beam.stats)
{
  n.mtx=length(beam.stats$beam.stats)
  pvals=vector("list",n.mtx)
  for (i in 1:n.mtx)
  {
    message(paste0("Computing feature p-values for stat matrix ",i," of ",n.mtx,": ",date()))
    beta=beam.stats$beam.stats[[i]]
    message(paste0("  This matrix has ",nrow(beta)," features."))
    B=rowSums(!is.na(beta))
    geq0=rowSums(beta[,-1]>=0,na.rm=T)
    leq0=rowSums(beta[,-1]<=0,na.rm=T)
    p=2*(pmin(geq0,leq0)+1)/(B+1)
    p[p>1]=1
    pi.hat=min(1,2*mean(p,na.rm=T))
    q=stats::p.adjust(p,"fdr")
    q=pi.hat*q

    clean.beta=clean_Bmtx(beta)
    beta.hat=rowMeans(clean.beta,na.rm=T)
    pvl.mtx=cbind(beta=beta.hat,p=p,q=q)
    rownames(pvl.mtx)=rownames(beta)

    mtx.name=beam.stats$beam.specs[i,"mtx"]
    ann.indx=which(beam.stats$beam.data$anns.mtch[,"mtx.data"]==mtx.name)
    ann.name=beam.stats$beam.data$anns.mtch[ann.indx,"mtx.anns"]
    ann.id=beam.stats$beam.data$anns.mtch[ann.indx,"id.clm"]
    ann.data=beam.stats$beam.data$mtx.anns[[ann.name]]

    ann.pvl=merge(ann.data,pvl.mtx,
                  by.x=ann.id,by.y="row.names",
                  all=T)

    pvals[[i]]=ann.pvl


  }
  names(pvals)=names(beam.stats$beam.stats)
  return(pvals)
}
