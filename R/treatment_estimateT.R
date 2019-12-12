
estimateTreatmentT <- function(dat, adjustVars,trt, glm.trt = NULL, SL.trt = NULL,
                              returnModels = FALSE, verbose = FALSE,
                              gtol = 1e-3, trtofTime = NULL, trtOfInterest = NULL, ...) {
n <- length(dat[,1])



##### Time varing confounding
uniqtrt <- unique(unlist(trtOfInterest[,-1]))
ntrt <- length(uniqtrt)

regime <- trtOfInterest[,-1]
trtofTime <- trtOfInterest[,1]
n_regime <- ncol(regime)

for (t.ind in trtofTime){


  if(ntrt == 1) {
    eval(parse(text = paste0("dat$g_", c(1:n_regime), "<- 1")))
  } else if (ntrt == 2){
    ## create data to fit the model
    id_include <- !is.na(trt[[paste0("trt_t", t.ind)]])
    trt_outcome <- trt[id_include, paste0("trt_t", t.ind)]


    # binarize the outcome
    thisY <- as.numeric(trt_outcome == max(trt_outcome))

    if (t.ind != 0){
      trt_include <-c("id", paste0("trt_t", trtofTime[trtofTime < t.ind]))
      trt_past <- trt[id_include, trt_include]

    } else { trt_past = data.frame(id = trt[id_include, "id"])}

    var_include <- NULL
    for (var.ind in trtofTime[trtofTime <= t.ind]){
      var_include <- c(var_include, colnames(adjustVars)[grepl(paste0("t",var.ind), colnames(adjustVars))])
    }

    var_past <- adjustVars[id_include,c("id", var_include)]

    trt_var_past <- merge(trt_past, var_past, by = "id")[,-1]

    # fit Super Learner
    #  NEED TO FIX SUPER LEARNER PART
    if(!is.null(SL.trt)) {
      if(class(SL.trt) != "SuperLearner") {
        trtMod <- SuperLearner::SuperLearner(Y = thisY, X = trt_var_past,
                                             newX = trt_var_past,
                                             SL.library = SL.trt,
                                             id = id_include, verbose = verbose,
                                             family = "binomial")
      } else {
        trtMod <- SL.trt
      }
      #dat[[paste0("g_",max(dat$trt))]] <- trtMod$SL.predict
      #dat[[paste0("g_",min(dat$trt))]] <- 1 - trtMod$SL.predict

    } else if(!is.null(glm.trt) & is.null(SL.trt)) {
      # set up model formula and data for the treatment regression
      #trt_form <- paste("thisY", "~", glm.trt, sep = " ")

      # fit GLM if Super Learner not requested
      if(!("glm" %in% class(glm.trt)) & !("speedglm" %in% class(glm.trt))) {
        # fit the treatment model

        if (is.list(glm.trt)){
          if ( length(glm.trt)  ==  length(trtofTime) ){
            trt_form <- paste("thisY", "~", glm.trt[[paste0("t", t.ind)]], sep = " ")
          } else{
            stop("need to specified regression form for treatment
                 at all the time points of assigning treatment")
          }
        } else {
          trt_form <- paste("thisY", "~", glm.trt, sep = " ")
        }

        #trt_form <- paste("thisY", "~", glm.trt, sep = " ")
        trt_data_in <- as.data.frame(cbind(trt_var_past, thisY))

        trtMod <- fast_glm(reg_form = stats::as.formula(trt_form),
                           data = trt_data_in,
                           family = stats::binomial())
      } else {
        trtMod <- glm.trt
      }
    }

    for (regime_ind in 1:n_regime){
      trt_data_pred <- trt_data_in
      if (t.ind != 0) {
        for(regime_time in 1:length(trtofTime[which(trtofTime < t.ind)])){
          t <- trtofTime[regime_time]
          trt_data_pred[,paste0("trt_t",t)] <- regime[regime_time,regime_ind]
        }
      }


      suppressWarnings(
        pred <- predict(trtMod, newdata = trt_data_pred, type = "response")
      )

      if(t.ind == 0){
        prev_g <- 1
      }else{
        prev_t <- max(trtofTime[trtofTime < t.ind ])
        prev_g <- dat[id_include, paste0("g_", regime_ind, "_t", prev_t)]
      }


      if(regime[which(trtofTime == t.ind ),regime_ind] == 1){
        dat[id_include, paste0("g_", regime_ind, "_t", t.ind)] <-
          pred*prev_g
      }else{
        dat[id_include, paste0("g_", regime_ind, "_t", t.ind)] <-
          (1-pred)*prev_g
      }


  } ##### need to add the case when trt is not binary


  # truncate propensities
  for(a in 1:n_regime){
    eval(parse(text = paste0("dat$g_", a, "_t",t.ind, "[dat$g_", a, "_t", t.ind,
                             "< gtol]<- gtol")))
  }

  suppressWarnings(
    ###
    pred_obsz <- predict(trtMod, newdata = trt_data_in, type = "response")
  )

  if(t.ind == 0){prev_g_obs = 1}else{
    prev_t_obs <- max(trtofTime[trtofTime < t.ind ])
    prev_g_obs <- dat[id_include, paste0("g_obsz_t", prev_t_obs)]}
    dat[id_include, paste0("g_obsz_t", t.ind)] <-
      ifelse(thisY == 1, pred_obsz, 1 - pred_obsz)*prev_g_obs

}

}


out <- list()
out$dat <- dat
out$trtMod <- NULL


if(returnModels) out$trtMod <- trtMod
return(out)

}
