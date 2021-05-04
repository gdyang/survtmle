estimateCensoringT <- function(dat, adjustVars,
                               glm.trt = NULL,
                               trt, t0,
                               SL.trt = NULL,
                               SL.ctime = NULL,
                               returnModels = FALSE,
                               verbose = FALSE, glm.family,
                               gtol = 1e-3,
                               trtofTime = NULL,
                               glm.ctime = NULL,
                               trtOfInterest = NULL, ...) {

  ##### Time varing confounding censoring

  regimen <- trtOfInterest[,-1]
  trtofTime <- trtOfInterest[,1]
  n_regimen <- ncol(regimen)


  if (!is.null(glm.family)) {
    glm_family <- parse(text = paste0("stats::", glm.family, "()"))
  }

  for (t.ind in c(0, seq_len(max(t0)))){

      ## create data to fit the model
      if (t.ind != 0){
        id_include <- (dat$ftime >= t.ind)
        ftype_include <- dat[id_include, "ftype"]
        ftime_include <- dat[id_include, "ftime"]


        trt_include <-c("id", paste0("trt_t", trtofTime[trtofTime < t.ind]))
        trt_past <- trt[id_include, trt_include]


        var_include <- NULL
        for (var.ind in trtofTime[trtofTime < t.ind]){
          var_include <- c(var_include,colnames(adjustVars)[grepl(paste0("t",var.ind), colnames(adjustVars))])
        }

        var_past <- adjustVars[id_include,c("id", var_include)]


        # binarize the outcome
        #thisY <- as.numeric(cen_outcome == 0)

        thisY <- ifelse((ftime_include == t.ind & ftype_include == 0) , 1, 0)
        trt_var_past <- merge(trt_past, var_past, by = "id")[,-1]

        # fit Super Learner
        #  NEED TO FIX SUPER LEARNER PART
        if(!is.null(SL.trt)) {
          if(class(SL.trt) != "SuperLearner") {
            if(!all(thisY == 0)){
              ctimeMod  <- SuperLearner::SuperLearner(Y = thisY, X = trt_var_past,
                                                   newX = trt_var_past,
                                                   SL.library = SL.ctime,
                                                   id = id_include, verbose = verbose,
                                                   family = "binomial", method = "method.CC_nloglik", cvControl = list(V = 5))
            }else{
              trtMod <- SL.trt
              ctimeMod <- "No censoring observed"
              class(ctimeMod) <- "noCens"
            }
            }else{# if input SLlibrary.time is Super Learner object, just use that
              ctimeMod <- SL.ctime
            }

          #dat[[paste0("g_",max(dat$trt))]] <- trtMod$SL.predict
          #dat[[paste0("g_",min(dat$trt))]] <- 1 - trtMod$SL.predict

        } else if(!is.null(glm.ctime) & is.null(SL.ctime)) {
          # set up model formula and data for the treatment regression
          #trt_form <- paste("thisY", "~", glm.trt, sep = " ")
          trt_data_in <- as.data.frame(cbind(trt_var_past, thisY))

          # fit GLM if Super Learner not requested
          if(!("glm" %in% class(glm.ctime)) & !("speedglm" %in% class(glm.ctime))) {
            # fit the treatment model
            if(!all(thisY == 0)) {
              if (is.list(glm.ctime)){
                if ( length(glm.ctime)  == max(t0)){
                ctimeForm <- stats::as.formula(sprintf("%s ~ %s", "thisY", glm.ctime[[paste0("t", t.ind)]]))
                } else{
                 stop("need to specified regression form for censoring at all the time points")
                }
              } else {
                ctimeForm <- stats::as.formula(sprintf("%s ~ %s", "thisY", glm.ctime))
              }
              ctimeMod <- fast_glm(reg_form = ctimeForm,
                                   data =  trt_data_in,
                                   family = eval(glm_family))
              if (unique(class(ctimeMod) %in% c("glm", "lm"))) {
                ctimeMod <- cleanglm(ctimeMod)
              }
            }else{
              dat[id_include, paste0("G_dC_t", t.ind)]  <- 1
              ctimeMod <- "No censoring observed"
              class(ctimeMod) <- "noCens"
            }
          } else {
            ctimeMod <- glm.ctime
          }



        }
       if(all(class(ctimeMod) %in% "noCens") ){
         pred <- 0
       }else{
        suppressWarnings(
          pred <- predict(ctimeMod, newdata = trt_data_in, type = "response")
        )
       }
        dat[id_include, paste0("G_obs_dC_t", t.ind)] <- (1-pred)* dat[id_include, paste0("G_obs_dC_t", t.ind-1)]



        for (regimen_ind in 1:n_regimen){
          if(all(class(ctimeMod) != "noCens")) {
            dat[id_include, paste0("G_", regimen_ind ,"dC_t", t.ind)] <- 1
            trt_data_pred <- trt_data_in

          for(regimen_time in 1:length(trtofTime[which(trtofTime < t.ind)])){
            t <- trtofTime[regimen_time ]
            trt_data_pred[,paste0("trt_t",t)] <- regimen[regimen_time ,regimen_ind]
            }


          suppressWarnings(
            pred <- predict(ctimeMod, newdata = trt_data_pred, type = "response")
          )

            dat[id_include, paste0("G_", regimen_ind ,"dC_t", t.ind)] <- (1-pred)* dat[id_include, paste0("G_", regimen_ind ,"dC_t", t.ind-1)]
          }else {
            dat[id_include, paste0("G_", regimen_ind ,"dC_t", t.ind)] <-  dat[id_include, paste0("G_", regimen_ind ,"dC_t", t.ind-1)]
          }
        }

        # truncate propensities
        for(a in 1:n_regimen){
          eval(parse(text = paste0("dat$G_", a, "dC_t",t.ind, "[dat$G_", a, "dC_t", t.ind,
                                   "< gtol]<- gtol")))
        }



      }else{
        dat[, paste0("G_obs_dC_t", t.ind)] <- 1
        for (regimen_ind in 1:n_regimen){
          dat[, paste0("G_", regimen_ind ,"dC_t", t.ind)] <- 1
        }

      }


  # make a column of observed a

  #ind <- rep(NA, n)
  #for(j in 1:ntrt){
  #  ind[dat$trt == uniqtrt[j]] <- which(colnames(dat) == paste0("g_",uniqtrt[j]))
  #}
  #dat$g_obsz <- dat[cbind(seq_along(ind),ind)]

  #out <- list()
  #out$dat <- dat
  #out$trtMod <- NULL






  }
  out <- list(dat = dat,
              ctimeMod = if(returnModels) {
                ctimeMod
              } else {
                NULL
              }
  )
}
