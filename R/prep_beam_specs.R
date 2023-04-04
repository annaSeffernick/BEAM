#' Prepare beam.specs data.frame for BEAM model fitting
#'
#' @param beam.data A beam.data object from prep_beam_data
#' @param endpts A vector of endpoint variable names in main.data
#' @param adjvars Default NULL, optional vector of adjustment variable names in main.data
#'
#' @return beam.specs, a data.frame specifying the omics-endpoint association models to be fit
#' @export
#'
#' @examples
#' data(clinf)
#' data(omicdat)
#' data(omicann)
#' data(setdat)
#' data(specs)
#' test.beam.data <- prep_beam_data(main.data=clinf, mtx.data=omicdat,
#'                                  mtx.anns=omicann, set.data=setdat,
#'                                  set.anns=NULL, n.boot=100, seed=123)
#' #Without adjustment
#' prep_beam_specs(beam.data=test.beam.data, endpts=c("MRD29", "OS", "EFS"))
#' # With adjustment
#' prep_beam_specs(beam.data=test.beam.data, endpts=c("OS", "EFS"), adjvars=c("MRD29"))
prep_beam_specs <- function(beam.data, endpts, adjvars=NULL){
  if(!inherits(beam.data, "beam.data"))
    stop("beam.data must be a beam.data object from prep_beam_data")
  mainData <- beam.data$main.data
  name.mtx <- names(beam.data$mtx.data)
  if(!all(endpts %in% colnames(mainData)))
    stop("endpts must be in column names of beam.data$main.data")
  if(is.null(adjvars)){
    beam.specs <- data.frame()
    for(j in 1:length(endpts)){
      temp.ep.name <- endpts[j]
      temp.ep.vec <- mainData[,grep(temp.ep.name, colnames(mainData))]
      temp.ep.type <- class(temp.ep.vec)
      #if(inherits(temp.ep.type, "Surv")){temp.mdl <- paste0("coxphf(", temp.ep.name, "~mtx.row,data=main.data)")}
      #if(inherits(temp.ep.type, "numeric")){temp.mdl <-paste0("lm(", temp.ep.name, "~mtx.row,data=main.data,model=T)")}
      temp.mdl <- ifelse(inherits(temp.ep.vec, "Surv"), paste0("coxphf(", temp.ep.name, "~mtx.row,data=main.data)"),
                         ifelse(inherits(temp.ep.vec, "numeric"), paste0("lm(", temp.ep.name, "~mtx.row,data=main.data,model=T)"),
                                ifelse(inherits(temp.ep.vec, "factor"), paste0("glm(", temp.ep.name,"~mtx.row,data=main.data,family='binomial')"), NULL)))
      for(k in 1:length(name.mtx)){
        temp.vec <- c(paste(name.mtx[k], temp.ep.name, sep="."), name.mtx[k], temp.mdl)
        beam.specs <- rbind.data.frame(beam.specs, temp.vec)
      }
    }
    colnames(beam.specs) <- c("name", "mtx", "mdl")
    return(beam.specs)
  }
  else{
    if(!all(adjvars %in% colnames(mainData)))
      stop("adjvars must be in column names of beam.data$main.data")
    adjvars.add <- paste(adjvars, collapse="+")
    beam.specs <- data.frame()
    for(j in 1:length(endpts)){
      temp.ep.name <- endpts[j]
      temp.ep.vec <- mainData[,grep(temp.ep.name, colnames(mainData))]
      temp.ep.type <- class(temp.ep.vec)
      temp.mdl <- ifelse(inherits(temp.ep.vec, "Surv"), paste0("coxphf(", temp.ep.name, "~mtx.row+",adjvars.add,",data=main.data)"),
                         ifelse(inherits(temp.ep.vec, "numeric"), paste0("lm(", temp.ep.name, "~mtx.row+", adjvars.add,",data=main.data,model=T)"),
                                ifelse(inherits(temp.ep.vec, "factor"), paste0("glm(", temp.ep.name,"~mtx.row+",adjvars.add,",data=main.data,family='binomial', model=T)"), NULL)))
      for(k in 1:length(name.mtx)){
        temp.vec <- c(paste(name.mtx[k], temp.ep.name, sep="."), name.mtx[k], temp.mdl)
        beam.specs <- rbind.data.frame(beam.specs, temp.vec)
      }
    }
    colnames(beam.specs) <- c("name", "mtx", "mdl")
    return(beam.specs)
  }
}
