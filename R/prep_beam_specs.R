#' Prepare beam.specs data.frame for BEAM model fitting
#'
#' @param beam.data A beam.data object from prep_beam_data
#' @param endpts A vector of endpoint variable names in main.data
#' @param adjvars Default NULL, optional vector of adjustment variable names in main.data
#' @param endptmdl Optional model specification data.frame with endpoint name column called "endpt" and model string column called "mdl"
#'
#' @return beam.specs, a data.frame specifying the omics-endpoint association models to be fit
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
#' #Without adjustment
#' prep_beam_specs(beam.data=test.beam.data, endpts=c("MRD29", "OS", "EFS"))
#' # With adjustment
#' prep_beam_specs(beam.data=test.beam.data, endpts=c("OS", "EFS"), adjvars=c("MRD29"))
prep_beam_specs <- function(beam.data, endpts, adjvars=NULL, endptmdl=NULL){
  # check that beam.data is a beam.data object
  if(!inherits(beam.data, "beam.data"))
    stop("beam.data must be  beam.data object from prep_beam_data")
  # check that endpts is a character vector
  if(!inherits(endpts, "character"))
    stop("endpts must be a character vector of endpoint names")
  # Extract main data (clinical data)
  mainData <- beam.data$main.data
  name.mtx <- names(beam.data$mtx.data)
  # Check that endpoints are present in clinical data
  if(!all(endpts %in% colnames(mainData)))
    stop("endpts must be in column names of beam.data$main.data")
  if(is.null(adjvars)){
    beam.specs <- data.frame()
    for(j in 1:length(endpts)){
      temp.ep.name <- endpts[j]
      temp.ep.vec <- mainData[,grep(temp.ep.name, colnames(mainData))]
      temp.ep.type <- class(temp.ep.vec)
      if(!is.null(endptmdl)){
        temp.mdl <- endptmdl[grep(temp.ep.name, endptmdl$endpt),]$mdl
        message(paste0("Fitting ", temp.mdl, " for ", temp.ep.name))
      }
      else if(inherits(temp.ep.vec, "Surv")){
        temp.mdl <- paste0("coxph(", temp.ep.name,"~mtx.row,data=main.data)")
        message(paste0(temp.ep.name," is survival endpoint, fitting coxph"))
      }
      else if(inherits(temp.ep.vec, "numeric")){
        uni.len <- length(unique(temp.ep.vec))
        if(uni.len==2){
          temp.mdl <- paste0("glm(", temp.ep.name,"~mtx.row,data=main.data,family='binomial')")
          message(paste0(temp.ep.name," is binary, fitting glm with family='binomial'"))
        }
        else{
          temp.mdl <- paste0("lm(", temp.ep.name,"~mtx.row,data=main.data,model=T)")
          message(paste0(temp.ep.name," is continuous, fitting lm"))
        }
      }
      else if(inherits(temp.ep.vec, "int")){
        uni.len <- length(unique(temp.ep.vec))
        if(uni.len==2){
          temp.mdl <- paste0("glm(", temp.ep.name,"~mtx.row,data=main.data,family='binomial')")
          message(paste0(temp.ep.name," is binary, fitting glm with family='binomial'"))
        }
        else{
          temp.mdl <- paste0("glm", temp.ep.name,"~mtx.row,data=main.data,family=poisson())")
          message(paste0(temp.ep.name," is count, fitting glm with family='poisson'"))
        }
      }
      else if(inherits(temp.ep.vec, "factor")){
        uni.lev <- length(levels(temp.ep.vec))
        if(uni.lev==2){
          temp.mdl <- paste0("glm(", temp.ep.name,"~mtx.row,data=main.data,family='binomial')")
          message(paste0(temp.ep.name," is binary, fitting glm with family='binomial'"))
        }
        else{
          stop(paste0(temp.ep.name," is a factor variable. Please input appropriate model using endptmdl."))
        }
      }
      else{
        stop(paste0(temp.ep.name, " is not a recognized variable type. Please input approriate model using endptmdl."))
      }
      for(k in 1:length(name.mtx)){
        temp.vec <- c(paste(name.mtx[k], temp.ep.name, sep="."), name.mtx[k], temp.mdl)
        beam.specs <- rbind.data.frame(beam.specs, temp.vec)
      }
    }
    colnames(beam.specs) <- c("name", "mtx", "mdl")
    return(beam.specs)
  }
  else{
    if(!inherits(adjvars, "character"))
      stop("adjvars must be a character vector of adjustment variable names")
    if(!all(endpts %in% colnames(mainData)))
      stop("adjvars must be in column names of beam.data$main.data")
    adjvars.add <- paste(adjvars, collapse="+")
    beam.specs <- data.frame()
    for(j in 1:length(endpts)){
      temp.ep.name <- endpts[j]
      temp.ep.vec <- mainData[,grep(temp.ep.name, colnames(mainData))]
      temp.ep.type <- class(temp.ep.vec)
      if(!is.null(endptmdl)){
        temp.mdl <- endptmdl[grep(temp.ep.name, endptmdl$endpt),]$mdl
        message(paste0("Fitting ", temp.mdl, " for ", temp.ep.name,"; Note that you must specify adjustment variables in mdl string."))
      }
      else if(inherits(temp.ep.vec, "Surv")){
        temp.mdl <- paste0("coxph(", temp.ep.name,"~mtx.row+",adjvars.add,",data=main.data)")
        message(paste0(temp.ep.name," is survival endpoint, fitting coxph"))
      }
      else if(inherits(temp.ep.vec, "numeric")){
        uni.len <- length(unique(temp.ep.vec))
        if(uni.len==2){
          temp.mdl <- paste0("glm(", temp.ep.name,"~mtx.row+",adjvars.add,",data=main.data,family='binomial')")
          message(paste0(temp.ep.name," is binary, fitting glm with family='binomial'"))
        }
        else{
          temp.mdl <- paste0("lm(", temp.ep.name,"~mtx.row+",adjvars.add,",data=main.data,model=T)")
          message(paste0(temp.ep.name," is continuous, fitting lm"))
        }
      }
      else if(inherits(temp.ep.vec, "int")){
        uni.len <- length(unique(temp.ep.vec))
        if(uni.len==2){
          temp.mdl <- paste0("glm(", temp.ep.name,"~mtx.row+",adjvars.add,",data=main.data,family='binomial')")
          message(paste0(temp.ep.name," is binary, fitting glm with family='binomial'"))
        }
        else{
          temp.mdl <- paste0("glm", temp.ep.name,"~mtx.row+",adjvars.add,"data=main.data,family=poisson())")
          message(paste0(temp.ep.name," is count, fitting glm with family='poisson'"))
        }
      }
      else if(inherits(temp.ep.vec, "factor")){
        uni.lev <- length(levels(temp.ep.vec))
        if(uni.lev==2){
          temp.mdl <- paste0("glm(", temp.ep.name,"~mtx.row+",adjvars.add,",data=main.data,family='binomial')")
          message(paste0(temp.ep.name," is binary, fitting glm with family='binomial'"))
        }
        else{
          stop(paste0(temp.ep.name," is a factor variable. Please input appropriate model using endptmdl."))
        }
      }
      else{
        stop(paste0(temp.ep.name, " is not a recognized variable type. Please input approriate model using endptmdl."))
      }
      for(k in 1:length(name.mtx)){
        temp.vec <- c(paste(name.mtx[k], temp.ep.name, sep="."), name.mtx[k], temp.mdl)
        beam.specs <- rbind.data.frame(beam.specs, temp.vec)
      }
    }
    colnames(beam.specs) <- c("name", "mtx", "mdl")
    return(beam.specs)
  }

}
