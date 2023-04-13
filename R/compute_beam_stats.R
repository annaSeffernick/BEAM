#' Compute bootstrap model coefficients for BEAM
#'
#' @param beam.data Result of prep.beam.data
#' @param beam.specs A data.frame of strings with columns name, mtx, mdl (string with R model with mtx.row)
#'
#' @return A beam.stats object, which is a list with beam.stats (the association matrices), the beam.specs, and the beam.data
#' @importFrom stats coef
#' @importFrom stats sd
#' @importFrom stats lm
#' @importFrom survival coxph
#' @importFrom coxphf coxphf
#' @export
#'
#' @examples
#' data(clinf)
#' data(omicdat)
#' data(omicann)
#' data(setdat)
#' test.beam.data <- prep_beam_data(main.data=clinf, mtx.data=omicdat,
#'                                  mtx.anns=omicann, set.data=setdat,
#'                                  set.anns=NULL, n.boot=100, seed=123)
#' specs <- prep_beam_specs(beam.data=test.beam.data, endpts=c("MRD29", "EFS", "OS"))
#' test.beam.stats <- compute_beam_stats(beam.data=test.beam.data, beam.specs=specs)
compute_beam_stats=function(beam.data, beam.specs)
{
  ##############################
  # Check data
  if(!inherits(beam.data, "beam.data"))
    stop("beam.data must be the result of prep.beam.data.")

  print(beam.data)
  print(beam.specs)

  ##################################
  # check the specs
  mtx.names=names(beam.data$mtx.data)
  main.clms=colnames(beam.data$main.data)
  spec.check=check_beam_specs(beam.specs,
                              mtx.names)

  n.spec=nrow(beam.specs)
  boot.index=beam.data$boot.index
  print(class(boot.index))
  main.data=beam.data$main.data

  beam.stats=vector("list",n.spec)
  for (i in 1:n.spec)                    # loop over specifications
  {
    message(paste("Working on analysis",
                  i,"of",n.spec,":",date()))
    print(beam.specs[i,c("mtx","mdl")])
    mtx.name=beam.specs[i,"mtx"]
    mtx=beam.data$mtx.data[[mtx.name]]
    mdl=beam.specs[i,"mdl"]
    x.clm=paste0(mtx.name,".clm")

    print(boot.index)
    beam.stats[[i]]=boot_model_coefs(boot.index,
                                     mtx,
                                     main.data,
                                     x.clm,
                                     mdl)
  }

  names(beam.stats)=beam.specs[,"name"]

  res=list(beam.stats=beam.stats,
           beam.specs=beam.specs,
           beam.data=beam.data)

  class(res)="beam.stats"

  return(res)
}

###############################################
# compute bootstrap model coefficient matrix

boot_model_coefs=function(boot.index,
                          X,
                          main.data,
                          x.clm,
                          mdl,
                          mess.freq=10)

{
  #print(dim(X))
  m=nrow(X)
  b=nrow(boot.index)
  #print(b)
  #print(head(boot.index))
  #print(dim(boot.index))
  res=matrix(NA,m,b)
  rownames(res)=rownames(X)
  colnames(res)=paste0("boot_",0:(b-1))
  for (i in 1:b)
  {
    if (((i-1)%%mess.freq)==0)
      message(paste("  Working on bootstrap",
                    i-1,"of",b,":",date()))
    boot.ind=boot.index[i,]
    #print(boot.ind)
    boot.data=main.data[boot.ind,]
    x.index=boot.data[,x.clm]
    print(x.index)
    boot.X=X[,x.index]
    print(dim(boot.X))
    res[,i]=row_model_coefs(boot.X,boot.data,mdl)
  }
  return(res)
}

#########################################
# For each row of X, fit a model and get the first coefficient

row_model_coefs=function(X,main.data,mdl)

{
  res=apply(X,1,model_coef,
            main.data=main.data,
            mdl=mdl)
  return(res)
}

############################################
# Get the first non-intercept coefficient for a fitted model

model_coef=function(x,main.data,mdl)

{
  res=NA
  sdx=stats::sd(x,na.rm=T)
  if (is.na(sdx)) return(0)
  if (sdx==0) return(0)

  mnx=mean(x,na.rm=T)
  main.data$mtx.row=(x-mnx)/sdx

  model.fit=try(eval(parse(text=mdl)))


  if (!inherits(model.fit[1],"try.error"))
  {
    beta=stats::coef(model.fit)
    int=(names(beta)=="(Intercept)")
    beta=beta[!int]
    res=beta[1]
  }

  #################################
  # special cases in coxph
  if(inherits(model.fit,"coxph"))
  {
    model.terms=get_model_terms(model.fit)

    evtm=model.fit$model[,model.terms[1]][,1]
    event=model.fit$model[,model.terms[1]][,2]
    x.temp=model.fit$model[,model.terms[2]]

    # no events
    if (is.null(event)) return(0)
    if (mean(event,na.rm=T)==0) return(0)

    # no variability in event time
    if ((stats::sd(evtm,na.rm=T)==0)&(mean(event,na.rm=T)==1)) return(0)

    # no variability in x
    if (stats::sd(x.temp,na.rm=T)==0) return(0)

    # check for monotone likelihood
    x.event=range(x.temp[event>0],na.rm=T)
    x.none=range(x.temp[event==0],na.rm=T)

    if (all(x.none>=x.event[2])) return(-Inf) # if x for cases with no event is always less than x for cases with event
    if (all(x.none<=x.event[1])) return(Inf)  #
  }

  #######################################
  # special cases in lm
  if (inherits(model.fit,"lm"))
  {
    model.terms=get_model_terms(model.fit)
    y.temp=model.fit$model[,model.terms[1]]
    x.temp=model.fit$model[,model.terms[2]]

    # no variability in x or y after removing NAs
    sdx=stats::sd(x.temp,na.rm=T)
    if (sdx==0) res=0
    sdy=stats::sd(y.temp,na.rm=T)
    if (sdy==0) res=0
  }


  if (inherits(model.fit[1],"try-error"))
  {
    message("Error fitting model, R call shown below:")
    print(mdl)
    message("Error fitting model, mtx.row summarized below:")
    print(summary(main.data$mtx.row))
    message("Error fitting model, summary of data shown below: ")
    print(summary(main.data))
  }

  return(res)
}

#######################################
# Get the terms of a model

get_model_terms=function(model.fit)

{
  model.terms=deparse(attr(model.fit$terms,"variables"))
  model.terms=unlist(strsplit(model.terms,split=","))
  model.terms=gsub("list(","",model.terms,fixed=T)
  model.terms=gsub(")","",model.terms,fixed=T)
  model.terms=gsub(" ","",model.terms)
  return(model.terms)
}

#######################################
# Modify coxphf

coxphf2=function(form,data)
{
  cox.res=coxph(form,data=data,x=T,model=T)
  all.zero=(sd(cox.res$model[,2])==0)
  if (all.zero) return(0)
  coxphf.res=try(coxphf(cox.res$formula,
                        data=cox.res$model))
  if (inherits(coxphf.res, "try-error")) return(0)
  else return(coxphf.res)
}


