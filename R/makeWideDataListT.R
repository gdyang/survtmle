#' Convert Long Form Data to List of Wide Form Data (TVC verison)
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

makeWideDataListT <- function(dat,
                             allJ,
                             uniqtrt,
                             trt,
                             trtOfInterest,
                             adjustVars,
                             wideData,
                             trtofTime,
                             msm.formula = NULL,
                             msmWeightList = NULL,
                             t0, ...){

  ### create widelist based on the numer of interested treatment regeime
  wideDataList <- vector(mode = "list", length = ncol(trtOfInterest))
  # if(is.null(msm.formula)){
  #   ind.ftype <- FALSE
  # } else {
  #   ind.ftype <- grepl("ftype", msm.formula)
  #   ind.ftime <- grepl("ftime", msm.formula)
  # }
  s.list <- t0
  dlNames <- colnames(dat)
  g_names <- dlNames[grepl("g_.*", dlNames)]
  # g_names <- g_names[-which(g_names == "g_obsz")]

  #dropVars <- c("trt", names(adjustVars),"ftime","ftype", g_names)
  ### extract baseline covariates
  baseVar <- adjustVars[,grepl("*_t0"  ,names(adjustVars))]


  wideDataList[[1]] <- data.frame(dat, trt[,-1],  baseVar, wideData[,-1])

  #colnames(wideDataList[[1]])[1] <- c("trt")
  #colnames(wideDataList[[1]])[2:(1 + ncol(adjustVars))] <- names(adjustVars)


  # set Nj0=0 for all j -- makes things easier to run in a loop later
  eval(parse(text = paste0(paste0("wideDataList[[1]]$N", allJ, ".0",
                                  collapse = "<-"),
                           "<- wideDataList[[1]]$C.0 <- 0")))



  wideDataList[2:length(wideDataList)] <- lapply( wideDataList[2:length(wideDataList)],
                                             function(x) {

                                               out <- data.frame(dat, baseVar, wideData[, -1])


                                               if(is.null(msm.formula)){
                                                 out[, paste0("C.", 1:max(t0))] <- 0
                                               }
                                               # set Nj0=0 for all j -- makes things easier to run in a loop later
                                               eval(parse(text = paste0(paste0("out$N", allJ, ".0", collapse = "<-"),
                                                                        "<- out$C.0 <- 0")))
                                               out
                                             })

  names(wideDataList) <- c("obs", colnames(trtOfInterest)[-1])


  # if(!ind.ftype){
  for(z in 2:ncol(trtOfInterest)){
    wideDataList[[z]] <- cbind(wideDataList[[z]], trt[,-1], dat[, grepl("g_obsz_*", colnames(dat))])
    for (trt.t in trtofTime){
    wideDataList[[z]][, paste0("trt_t", trt.t)] <- ifelse(!is.na(wideDataList[[z]][, paste0("trt_t", trt.t)]),
                                                              trtOfInterest[which(trtofTime == trt.t), z],
                                                              NA)
    }
  }
  # } else {
  #   for(z in uniqtrt){
  #     for (j in allJ){
  #       wideDataList[[paste0("Z", z, "J", j)]]$trt <- z
  #       wideDataList[[paste0("Z", z, "J", j)]]$g_obsz <- dat[[paste0("g_",z)]]
  #     }
  #   }
  # }

  if(is.null(msm.formula)){
    wideDataList <- lapply(wideDataList, function(x){
      # make clever covariates


###### do we need to create clever covariates for  [[paste0("H",z,".",t)]] for each regime?
      for(z in uniqtrt) {
        for(t in 1:max(t0)) {

          x[[paste0("H",z,".",t)]] <-
            (x$trt==z & x[[paste0("C.",t-1)]]==0) / (x[[paste0("G_dC.",t)]]*x[[paste0("g_obsz_t.",t)]])
        }
        x[[paste0("H",z,".0")]] <- (x$trt==z) / x[[paste0("g_obsz_t.",t)]]
      }
      x
    })
  }else{



    #create indicator for patient whether follow r regimen
    #### missing treated as not following throgh
    wideDataList <- lapply(wideDataList, function(wdl){
      for (r.ind in seq_len(ncol(trtOfInterest)-1)){
        wdl[[paste0("regimen", r.ind)]] <-  apply(wdl[,grepl("trt_t*", colnames(wdl))], 1, function(x){
          as.numeric(sum( x == trtOfInterest[[paste0("regimen", r.ind)]]) == length(trtofTime)) } )

        wdl[[paste0("regimen", r.ind)]][is.na(wdl[[paste0("regimen", r.ind)]])] <- 0
      }
      wdl
    })

      for (r.ind in seq_len(ncol(trtOfInterest)-1)){
        wideDataList[[r.ind+1]][[paste0("regimen", r.ind)]] <- 1
      }



    obs <- wideDataList[[1]]
    temp.obs.fill <- obs[seq_len(length(allJ)),]
    temp.obs.fill$ftype <- c(allJ)
    temp.obs.fill$ftime <- 0

    cfact <- wideDataList[[2]]
    temp.cfact.fill <- cfact[seq_len(length(allJ)),]
    temp.cfact.fill$ftype <-  c(allJ)
    temp.cfact.fill$ftime <- 0

    # TO DO: This breaks with factor(ftype)
    wideDataList <- mapply(wdl = wideDataList, mw = msmWeightList, FUN = function(wdl, mw){
      wdl.temp <- wdl
    seq_regmen <- seq_len(ncol(trtOfInterest)-1)
    for (r in seq_regmen){
      for(j in allJ){
        for (s in s.list){
          wdl$ftype <- j
          wdl$ftime <- s
          wdl.temp$ftype <- j
          wdl.temp$ftime <- s
          if( identical(colnames(wdl.temp), colnames(temp.obs.fill))){
            temp.fill <- temp.obs.fill} else{temp.fill <-temp.cfact.fill}
          wdl.new <- rbind(temp.fill, wdl.temp)


###### Do we want to make msm.formula into a list
        wdl.new[[paste0("regimen", r)]] <- 1
        wdl.new[[paste0("regimen", seq_regmen[which(seq_regmen != r)])]] <- 0
        msmModelMatrix <- stats::model.matrix(as.formula(paste0("N1.0 ~ ", msm.formula)), data = wdl.new)[-seq_len(nrow(temp.fill)),]
        msm.p <- dim(msmModelMatrix)[2]


          for(t in 1:s){

##### create indicator for following the treatment regimen till t-1
            ## check !
           past_t <- trtofTime[trtofTime <= t-1]
           ind.r.t_minus1 <-  apply(wdl[,paste0("trt_t", past_t), drop = FALSE], 1, function(x){
              as.numeric(sum((x == trtOfInterest[(trtofTime <= t-1), paste0("regimen", r)] ) ) == length(past_t)) } )

           ind.r.t_minus1[is.na(ind.r.t_minus1)] <- 0

            for(p in 1:msm.p){
              if ( (t-1) %in% trtofTime){
                ##### j: failure type, p: number of covariates, t: time, s: t0, r: treatment regime
                  wdl[[paste0("H", j, ".", p,".",t, ".", s, ".", r ,".obs")]] <-
                  as.numeric(msmModelMatrix[,p] * mw * as.numeric(wdl[,paste0("C.",t)]==0)* ind.r.t_minus1/
                               (wdl[[paste0("G_",r,"dC_t",t)]] * wdl[[paste0("g_",r,"_t",t-1)]]))
             #     as.numeric(msmModelMatrix[,p] * mw * ind.r.t_minus1/ (wdl[[paste0("G_",r,"dC_t",t)]] * wdl[[paste0("g_",r,"_t",t-1)]]))
                  wdl[[paste0("H", j, ".", p,".",t, ".", s, ".", r , ".pred")]] <-
                  as.numeric(msmModelMatrix[,p] * mw / (wdl[[paste0("G_",r,"dC_t",t)]] * wdl[[paste0("g_",r,"_t",t-1)]]))

              }else{
              #### if the time if not in trtofTime
                prev_time_of_trt <- max(trtofTime[trtofTime < t-1])

              #### H for each regimen
              wdl[[paste0("H", j, ".", p,".",t, ".", s, ".", r ,".obs")]] <-
                as.numeric(msmModelMatrix[,p] * mw * as.numeric(wdl[,paste0("C.",t)]==0)*ind.r.t_minus1/
                             (wdl[[paste0("G_",r,"dC_t",t)]]* wdl[[paste0("g_",r,"_t",prev_time_of_trt)]]))
              wdl[[paste0("H", j, ".", p,".",t, ".", s, ".", r , ".pred")]] <-
                as.numeric(msmModelMatrix[,p] * mw / (wdl[[paste0("G_",r,"dC_t",t)]]*wdl[[paste0("g_",r,"_t",prev_time_of_trt)]]))
              }
            }

          }

        } # end of t0 loop

        wdl <- wdl[,-which(colnames(wdl) == "ftype")]
      } # end of failure type
      } # end of treatment regimen

      wdl[,-which(colnames(wdl)=="ftime")]
    }, SIMPLIFY = FALSE)


  }

  return(wideDataList)
}
