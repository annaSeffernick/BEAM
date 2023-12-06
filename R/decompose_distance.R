#' Decompose the distance between observed and null
#'
#' Use distance decomposition to identify relative contribution of each omic feature/endpoint associaiton
#'
#' @param beam.stats A beam.stats object from compute_beam_stats function
#' @param z Logical indicating whether to z-scale each vector of one coefficient estimate across bootstraps before analysis
#' @param mess.freq  Message frequency; default 25.
#'
#' @return A list of dataframes. List names correspond to set names. Each data.frame has Association column with feature name and association type and Contribution column with relative contribution, ordered greatest to least contribution.
#' @importFrom stats prcomp
#' @importFrom stats predict
#' @importFrom stats sd
#' @export
#'
#' @examples
#' test.beam.data <- prep_beam_data(main.data=clinf, mtx.data=omicdat,
#'                                  mtx.anns=omicann, set.data=setdat,
#'                                  set.anns=NULL, n.boot=10, seed=123)
#' specs <- prep_beam_specs(beam.data=test.beam.data, endpts=c("MRD29", "EFS", "OS"), firth=TRUE)
#' test.beam.stats <- compute_beam_stats(beam.data=test.beam.data, beam.specs=specs)
#' test.dist <- decompose_distance(test.beam.stats, z=TRUE)
decompose_distance <- function(beam.stats, z=TRUE, mess.freq=25)
{
  ######################################################
  # Identify the matrices with analyses specified
  beam.specs=beam.stats$beam.specs
  specs.mtx=beam.specs[,"mtx"]
  stat.ids=names(beam.stats$beam.stats)

  #################################################
  # Extract the set data with analyses specified
  message(paste0("Preparing bootstrap results for calculating feature set p-values: ",date()))
  set.data=beam.stats$beam.data$set.data

  keep.set.data=is.element(set.data$mtx.id,specs.mtx)
  set.data=set.data[keep.set.data,]


  ord=order(set.data$mtx.id)
  set.data=set.data[ord,]
  a=nrow(set.data)
  new.sect=which((set.data$mtx.id[-1]!=set.data$mtx.id[-a]))
  row.start=c(1,new.sect+1)
  row.end=c(new.sect,a)
  mtx.index=cbind.data.frame(row.start=row.start,
                             row.end=row.end,
                             mtx.id=set.data$mtx.id[row.end])

  ###############################################
  # Merge set data with analysis stat matrices

  message(paste0("Finding stats for each data matrix:",date()))
  set.mtch=NULL
  for (i in 1:nrow(mtx.index))
  {
    message(paste0("  Finding stats for data matrix ",mtx.index$mtx.id[i],": ",date()))
    spec.rows=which(specs.mtx==mtx.index$mtx.id[i])
    for (j in spec.rows)
    {
      message(paste0("   Finding features with with ",stat.ids[j]," stats: ",date()))
      row.ids=row.names(beam.stats$beam.stats[[j]])
      stat.rows=cbind(mtx.row=paste0(mtx.index$mtx.id[i],"_",row.ids),
                      stat.id=stat.ids[j],
                      row.id=row.ids)
      set.mtch=rbind.data.frame(set.mtch,stat.rows)
    }
  }
  message(paste0("Found ",nrow(set.mtch)," rows of stats: ",date()))

  message(paste0("Merging stats with feature-sets: ",date()))
  set.data$mtx.row=paste0(set.data$mtx.id,"_",set.data$row.id)
  set.mtch=merge(set.data[c("mtx.row","set.id")],
                 set.mtch,
                 by="mtx.row")
  message(paste0("Merged feature-set stat rows: ",nrow(set.mtch)))

  ##########################################
  # Index feature-set stat rows
  message(paste0("Ordering and indexing feature sets: ",date()))
  ord=order(set.mtch$set.id,
            set.mtch$stat.id)
  set.mtch=set.mtch[ord,]

  m=nrow(set.mtch)
  new.sect=which((set.mtch$set.id[-1]!=set.mtch$set.id[-m])|
                   (set.mtch$stat.id[-1]!=set.mtch$stat.id[-m]))
  row.start=c(1,new.sect+1)
  row.end=c(new.sect,m)

  set.stat.index=cbind.data.frame(row.start=row.start,
                                  row.end=row.end,
                                  set.id=set.mtch$set.id[row.start],
                                  stat.id=set.mtch$stat.id[row.start])

  m=nrow(set.stat.index)
  new.sect=which(set.stat.index$set.id[-1]!=set.stat.index$set.id[-m])
  row.start=c(1,new.sect+1)
  row.end=c(new.sect,m)
  set.index=cbind.data.frame(row.start=row.start,
                             row.end=row.end,
                             set.id=set.stat.index$set.id[row.start])

  ########################################
  # clean up bootstrap coefficient matrices
  message(paste0("  Cleaning up beam.stat matrices:",date()))
  mtx.list=names(beam.stats$beam.stats)
  n.mtx=length(mtx.list)
  for (i in 1:n.mtx)
  {
    message(paste0("   Working on matrix ",i," of ",n.mtx,": ",date()))
    B.mtx=beam.stats$beam.stats[[i]]
    B.cln=clean_Bmtx(B.mtx)
    beam.stats$beam.stats[[i]]=B.cln
  }


  print(head(set.index))
  print(head(set.stat.index))
  print(head(set.mtch))

  ##############################################
  # Extract bootstrap coefs and compute p-value for each set
  r=nrow(set.index)
  #p.set=rep(NA,r)
  features=rep("",r)
  #origin.to.center=rep(NA,r)
  #mean.from.center=rep(NA,r)
  #distance.ratio=rep(NA,r)
  dist.decomp.list <- list()

  for (i in 1:r)
  {
    if (((i-1)%%mess.freq)==0)
    {
      message(paste0("Computing p-value for feature set ",i," of ",r,": ",date()))
      message(set.index[i,])
    }

    ind1.start=set.index$row.start[i]
    ind1.end=set.index$row.end[i]
    ind2.start=set.stat.index$row.start[ind1.start]
    ind2.end=set.stat.index$row.end[ind1.end]
    ft.str=paste(set.stat.index$row.end[ind1.start:ind1.end]-
                   set.stat.index$row.start[ind1.start:ind1.end]+1,
                 set.stat.index$stat.id[ind1.start:ind1.end])
    features[i]=paste(ft.str,collapse=" | ")
    m=ind2.end-ind2.start+1
    B=NULL
    #obs.vec <- NULL
    for (j in ind1.start:ind1.end)
    {
      stat.id=set.stat.index$stat.id[j]
      row.start=set.stat.index$row.start[j]
      row.end=set.stat.index$row.end[j]
      row.ids=set.mtch$row.id[row.start:row.end]
      #b.mtx=beam.stats$beam.stats[[stat.id]][row.ids,-1]
      #obs.vec <- c(obs.vec, beam.stats$beam.stats[[stat.id]][row.ids, 1])
      b.mtx=beam.stats$beam.stats[[stat.id]][row.ids,] # keep observed value as first column
      names.b.mtx <- c(rownames(B), paste(row.ids, stat.id, sep="_"))
      B=rbind(B,b.mtx)
      rownames(B) <- names.b.mtx

    }
    B=t(B)                          # transpose for distance calculations
    #clm.sd=apply(B,2,sd)
    #clm.sd=apply(B[-1,],2,sd)
    #if (((i-1)%%mess.freq)==0) message(summary(clm.sd))
    dist.res <- decomp.dis(B, z=z)
    dist.decomp.list[[i]] <- dist.res
    names(dist.decomp.list)[i] <- set.index[i,"set.id"]
  }
  message(paste0("Finished computing distance decomposition at: ",date()))
  res <- dist.decomp.list
  return(res)
}

#B.mtx <- B
decomp.dis <- function(B.mtx,     # bootstrap matrix with observed result in row 1 - make sure it has names
                       z=TRUE)       # indicates whether z-scale each vector of one coefficient estimate across bootstraps before analysis
{
  cent=B.mtx[1,]                             # observed result is first row, it will be the center for PCA calculations
  B0=B.mtx[-1,]                              # remove observed result
  B.mtx=B0                                   # bootstrap association matrix
  b0=b=nrow(B.mtx)                           # number of bootstraps

  # Check column variability in B.mtx
  col.sd <- apply(B.mtx, 2, stats::sd)
  # If there is a constant column, scale the other columns by hand to avoid prcomp error
  if(any(col.sd==0)&z){
    B.mtx.scale <- B.mtx
    for(i in 1:ncol(B.mtx.scale)){
      if(col.sd[i]!=0){B.mtx.scale[,i] <- (B.mtx[,i]-cent[i])/f(B.mtx[,i]-cent[i])}
    }
    pca.res=stats::prcomp(B.mtx.scale,
                          center=F,
                          scale.=F)
  }
  else{
    pca.res=stats::prcomp(B.mtx,                      # principal components for all bootstraps
                          center=cent,                # use the observed result as the center
                          scale.=z)               # user option for rescaling PCs before computing distance

  }

  null=matrix(0,1,ncol(B.mtx))               # define vector for null (origin: all coefs = 0)
  colnames(null)=colnames(B.mtx)             # assign column names
  pca.null=stats::predict(pca.res,null)             # project original null into PC space
  pca.mtx=rbind(pca.res$x,pca.null)          # matrix of PCs (bootstraps and null)

  pca.dist=rowSums(pca.mtx^2)                # distance of null and each bootstrap from observed
  p.indx=which(pca.dist[-(b+1)]>=pca.dist[(b+1)])    # index of bootstraps farther from observed than null
  dist.to.null <- pca.dist[(b+1)]
  mean.from.obs <- mean(pca.dist[-(b+1)])
  # Distance decomposition
  pca.rel <- pca.mtx[(b+1),]^2/pca.dist[(b+1)]
  pca.load <- pca.res$rotation
  pca.wt <- t(abs(pca.load)/colSums(abs(pca.load)))
  dist.decomp <- pca.rel %*% pca.wt
  dist.decomp.df <- data.frame(Association=colnames(dist.decomp), Contribution=t(dist.decomp))
  dist.decomp.df.or <- dist.decomp.df[order(dist.decomp.df$Contribution, decreasing=TRUE),]
  return(dist.decomp.df.or)
}
