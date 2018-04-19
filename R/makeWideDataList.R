#' Convert Long Form Data to List of Wide Form Data
#'
#' The function takes a \code{data.frame} and \code{list} consisting of short
#' and long format right-censored failure times. The function reshapes the long
#' format into the wide format needed for calls to \code{mean_tmle}. The list
#' returned by the function will have number of entries equal to
#' \code{length(trtOfInterest) + 1}. The first will contain the observed
#' \code{trt} columns and will set \code{C.t} (the censoring counting process)
#' equal to the observed value of censoring. The subsequent entries will set
#' \code{trt} equal to each level of \code{trtOfInterest} and set \code{C.t} to
#' zero for everyone.
#'
#' @param dat The short form \code{data.frame}
#' @param allJ Numeric vector indicating the labels of all causes of failure.
#' @param uniqtrt The values of \code{trtOfInterest} passed to \code{mean_tmle}.
#' @param adjustVars A data.frame of adjustment variables that will be used in
#'        estimating the conditional treatment, censoring, and failure (hazard
#'        or conditional mean) probabilities.
#' @param dataList A list of long format \code{data.frame} objects. See
#'        \code{?makeDataList} for more details on formatting.
#' @param msm.formula A valid right-hand-side of a formula that can include 
#'        variables \code{trt} and \code{colnames(adjustVars)}
#' @param msmWeightList A list of weights in same format as this list will be (
#'        i.e., first entry corresponding to observed trt, latter to set values
#'        of trt)
#' @param t0 The timepoint at which \code{survtmle} was called to evaluate.
#' @param ... Other arguments. Not currently used.
#'
#' @importFrom stats reshape
#'
#' @return A list of \code{data.frame} objects as described above.
#'

makeWideDataList <- function(dat,
                             allJ,
                             uniqtrt,
                             adjustVars,
                             dataList,
                             msm.formula = NULL, 
                             msmWeightList = NULL, 
                             t0, ...) {
  wideDataList <- vector(mode = "list", length = length(dataList))
  ftype.ind <- grepl("ftype", msm.formula)
  dlNames <- colnames(dataList[[1]])
  g_names <- dlNames[grepl("g_.*", dlNames)]
  # g_names <- g_names[-which(g_names == "g_obsz")]
  
  dropVars <- c("trt", names(adjustVars),"ftime","ftype", g_names)

  wideDataList[[1]] <- data.frame(dat$trt, dat[, names(adjustVars)], g_obsz = dat$g_obsz,
                                  stats::reshape(dataList[[2]][, !(names(dataList[[2]]) %in%
                                                             dropVars)],
                                                 direction = "wide",
                                                 timevar = "t", idvar = "id"))
  colnames(wideDataList[[1]])[1] <- c("trt")
  colnames(wideDataList[[1]])[2:(1 + ncol(adjustVars))] <- names(adjustVars)
  
  
  # set Nj0=0 for all j -- makes things easier to run in a loop later
  eval(parse(text = paste0(paste0("wideDataList[[1]]$N", allJ, ".0",
                                 collapse = "<-"),
                          "<- wideDataList[[1]]$C.0 <- 0")))
  
  
 
  wideDataList[2:length(dataList)] <- lapply(dataList[2:length(dataList)],
                                             function(x) {
   if(!msm.ftype){
     out <- data.frame(dat[, names(adjustVars)],
                       stats::reshape(x[, !(names(x) %in% dropVars)],
                                      direction = "wide", timevar = "t", idvar = "id")
                       ,row.names = NULL)
     }else { 
     out <- data.frame(x$ftype , dat[,  names(adjustVars)],
                         stats::reshape(x[, !(names(x) %in% dropVars)],
                                        direction = "wide", timevar = "t", idvar = "id")
                         ,row.names = NULL)
       }
    
     if(is.null(msm.formula)){
        out[, paste0("C.", 1:t0)] <- 0
      }
    names(out)[1:(ncol(adjustVars)+1)] <- c("ftype", names(adjustVars))
    # set Nj0=0 for all j -- makes things easier to run in a loop later
    eval(parse(text = paste0(paste0("out$N", allJ, ".0", collapse = "<-"),
                            "<- out$C.0 <- 0")))
    out
  })
  
  if(!ftype.ind){
    names(wideDataList) <- c("obs", uniqtrt)  
  } else {
    names(wideDataList) <- c("obs",paste0("Z", rep(uniqtrt, each = length(allJ)),
                                            "J", rep(allJ, length(uniqtrt))))}

  if(!ftype.ind){
  for(z in uniqtrt){
    wideDataList[[paste0(z)]]$trt <- z
    wideDataList[[paste0(z)]]$g_obsz <- dat[[paste0("g_",z)]]
    }
  } else {
    for(z in uniqtrt){
      for (j in allJ){
        wideDataList[[paste0("Z", z, "J", j)]]$trt <- z
        wideDataList[[paste0("Z", z, "J", j)]]$g_obsz <- dat[[paste0("g_",z)]] 
      }
    }
  }
  
  if(is.null(msm.formula)){
    wideDataList <- lapply(wideDataList, function(x){
      # make clever covariates
      for(z in uniqtrt) {
        for(t in 1:t0) {
          x[[paste0("H",z,".",t)]] <- 
            (x$trt==z & x[[paste0("C.",t-1)]]==0) / (x[[paste0("G_dC.",t)]]*x$g_obsz)
        }
          x[[paste0("H",z,".0")]] <- (x$trt==z) / x$g_obsz
      }
      x
    })
  }else{
    if(!msm.ftype){
      wideDataList <- mapply(wdl = wideDataList, mw = msmWeightList, FUN = function(wdl, mw){
      # browser()
      msmModelMatrix <- model.matrix(as.formula(paste0("N1.0 ~ ",msm.formula)), data = wdl)
      msm.p <- dim(msmModelMatrix)[2]
      for(t in 1:t0){
        for(p in 1:msm.p){
          wdl[[paste0("H",p,".",t,".obs")]] <- 
            as.numeric(msmModelMatrix[,p] * mw * as.numeric(wdl[,paste0("C.",t-1)]==0) / (wdl[[paste0("G_dC.",t)]] * wdl$g_obsz))
          wdl[[paste0("H",p,".",t,".pred")]] <- 
            as.numeric(msmModelMatrix[,p] * mw / (wdl[[paste0("G_dC.",t)]] * wdl$g_obsz))
          
        }
      }
      for(p in 1:msm.p){
        wdl[[paste0("H",p,".",0,".obs")]] <- as.numeric(msmModelMatrix[,p] * mw / wdl$g_obsz)
        wdl[[paste0("H",p,".",0,".pred")]] <- as.numeric(msmModelMatrix[,p] * mw / wdl$g_obsz)
      }
      wdl
      }, SIMPLIFY = FALSE)} else {  
        
        
        ObswideDataList <- wideDataList[[1]]
        for(j in allJ){
          ObswideDataList$ftype <- j
          Obs.msmModelMatrix <- model.matrix(as.formula(paste0("N1.0 ~ ",msm.formula)), data = ObswideDataList)
          msm.p <- dim(Obs.msmModelMatrix)[2]
          for(t in 1:t0){
          for(p in 1:msm.p){
            ObswideDataList[[paste0("H", j, ".", p,".",t,".obs")]] <- 
              as.numeric(Obs.msmModelMatrix[,p] * mw * as.numeric(ObswideDataList[,paste0("C.",t-1)]==0) / (ObswideDataList[[paste0("G_dC.",t)]] * ObswideDataList$g_obsz))
            ObswideDataList[[paste0("H", j, ".", p,".",t,".pred")]] <- 
              as.numeric(Obs.msmModelMatrix[,p] * mw / (ObswideDataList[[paste0("G_dC.",t)]] * ObswideDataList$g_obsz))
          }
          }    
          ObswideDataList[-which(colnames(ObswideDataList) == "ftype")]
        }
        wideDataList[[1]] <- ObswideDataList
        wideDataList[-1] <- mapply(wdl = wideDataList[-1], mw = msmWeightList[-1], FUN = function(wdl, mw){
        msmModelMatrix <- model.matrix(as.formula(paste0("N1.0 ~ ",msm.formula)), data = wdl)
        msm.p <- dim(msmModelMatrix)[2]
        for(t in 1:t0){
          for(p in 1:msm.p){
            j <- (wdl$ftype)[1]
            wdl[[paste0("H", j, ".", p,".",t,".obs")]] <- 
              as.numeric(msmModelMatrix[,p] * mw * as.numeric(wdl[,paste0("C.",t-1)]==0) / (wdl[[paste0("G_dC.",t)]] * wdl$g_obsz))
            wdl[[paste0("H", j, ".", p,".",t,".pred")]] <- 
              as.numeric(msmModelMatrix[,p] * mw / (wdl[[paste0("G_dC.",t)]] * wdl$g_obsz))
            
          }
        }
        for(p in 1:msm.p){
          j <- (wdl$ftype)[1]
          wdl[[paste0("H", j, ".", p,".",0,".obs")]] <- as.numeric(msmModelMatrix[,p] * mw / wdl$g_obsz)
          wdl[[paste0("H", j, ".", p,".",0,".pred")]] <- as.numeric(msmModelMatrix[,p] * mw / wdl$g_obsz)
        }
        wdl
      }, SIMPLIFY = FALSE)
      
    }
    
    
  }

  return(wideDataList)
}
