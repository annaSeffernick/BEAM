#' Print summary information about a beam.data object
#'
#' @param beam.data A beam.data object from prep_beam_data()
#'
#' @return Messages about the beam.data object
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
#' print_beam_data(test.beam.data)
print_beam_data=function(beam.data)

{
  cat(paste0("main.data: ",
             nrow(beam.data$main.data)," rows and ",
             ncol(beam.data$main.data)," columns. \n \n"))

  n.mtx=length(beam.data$mtx.data)

  cat(paste0("mtx.data: \n"))

  mtx.names=names(beam.data$mtx.data)
  for (i in 1:n.mtx)
  {
    cat(paste0("  mtx.data ",mtx.names[i],": ",
               ncol(beam.data$mtx.data[[i]])," columns linked to ",
               sum(!is.na(beam.data$main.data[,paste0(mtx.names[i],".clm")])),
               " rows of main.data. \n"))
  }

  cat(paste0("\n","mtx.anns: \n"))
  n.anns=length(beam.data$mtx.anns)
  if (n.anns>0)
  {
    ann.names=names(beam.data$mtx.anns)
    for (i in 1:n.anns)
    {
      cat(paste0("  ",ann.names[i],": ",
                 nrow(beam.data$mtx.anns[[i]])," rows and ",
                 ncol(beam.data$mtx.anns[[i]])," columns. \n"))
    }
  }

  cat(paste0("\n","anns.mtch: \n"))
  print(beam.data$anns.mtch)

  cat(paste0("\n","set.data: ",
             nrow(beam.data$set.data)," rows assigning sets to data.mtx rows."))


  cat(paste0("\n \n","set.anns: ",
             nrow(beam.data$set.anns)," rows of set annotations."))

  cat(paste0("\n \n","boot.index: ",
             nrow(beam.data$boot.index)," rows and ",
             ncol(beam.data$boot.index)," columns of bootstrap indices. \n"))


}
