#' Estimation for the Method of Iterated Means Time-varing Confounding
#'
#' This function computes an estimate of the G-computation regression at a
#' specified time \code{t} using either \code{glm} or \code{SuperLearner}. The
#' structure of the function is specific to how it is called within
#' \code{mean_tmle}. In particular, \code{wideDataList} must have a very
#' specific structure for this function to run properly. The list should consist
#' of \code{data.frame} objects. The first should have all rows set to their
#' observed value of \code{trt}. The remaining should in turn have all rows set
#' to each value of \code{trtOfInterest} in the \code{survtmle} call. Currently
#' the code requires each \code{data.frame} to have named columns for each name
#' in \code{names(adjustVars)}, as well as a column named \code{trt}. It must
#' also have a columns named \code{Nj.Y} where j corresponds with the numeric
#' values input in \code{allJ}. These are the indicators of failure due to the
#' various causes before time \code{t} and are necessary for determining who to
#' include in the regression. Similarly, each \code{data.frame} should have a
#' column call \code{C.Y} where Y is again \code{t - 1}, so that right censored
#' observations are not included in the regressions. The function will fit a
#' regression with \code{Qj.star.t+1} (also needed as a column in
#' \code{wideDataList}) on functions of \code{trt} and \code{names(adjustVars)}
#' as specified by \code{glm.ftime} or \code{SL.ftime}.
#'
#' @param wideDataList A list of \code{data.frame} objects.
#' @param t The timepoint at which to compute the iterated mean.
#' @param whichJ Numeric value indicating the cause of failure for which
#'        regression should be computed.
#' @param allJ Numeric vector indicating the labels of all causes of failure.
#' @param t0 The timepoint at which \code{survtmle} was called to evaluate.
#'        Needed only because the naming convention for the regression if
#'        \code{t == t0} is different than if \code{t != t0}.
#' @param adjustVars Object of class \code{data.frame} that contains the
#'        variables to adjust for in the regression.
#' @param SL.ftime A character vector or list specification to be passed to the
#'        \code{SL.library} argument in the call to \code{SuperLearner} for the
#'        outcome regression (either cause-specific hazards or conditional mean).
#'        See \code{?SuperLearner} for more information on how to specify valid
#'        \code{SuperLearner} libraries. It is expected that the wrappers used
#'        in the library will play nicely with the input variables, which will
#'        be called \code{"trt"} and \code{names(adjustVars)}.
#' @param glm.ftime A character specification of the right-hand side of the
#'        equation passed to the \code{formula} option of a call to \code{glm}
#'        for the outcome regression (either cause-specific hazards or
#'        conditional mean). Ignored if \code{SL.ftime != NULL}. Use \code{"trt"}
#'        to specify the treatment in this formula (see examples). The formula
#'        can additionally include any variables found in
#'        \code{names(adjustVars)}.
#' @param verbose A boolean indicating whether the function should print
#'        messages to indicate progress.
#' @param returnModels A boolean indicating whether to return the
#'        \code{SuperLearner} or \code{glm} objects used to estimate the
#'        nuisance parameters. Must be set to \code{TRUE} if the user plans to
#'        use calls to \code{timepoints} to obtain estimates at times other than
#'        \code{t0}. See \code{?timepoints} for more information.
#' @param bounds A list of bounds to be used when performing the outcome
#'        regression (Q) with the Super Learner algorithm. NOT YET IMPLEMENTED.
#' @param ... Other arguments. Not currently used.
#'
#' @importFrom SuperLearner SuperLearner SuperLearner.CV.control
#' @importFrom stats as.formula predict model.matrix optim glm binomial gaussian
#' @importFrom speedglm speedglm
#'
#' @return The function then returns a list that is exactly the same as the
#'         input \code{wideDataList}, but with a column named \code{Qj.t} added
#'         to it, which is the estimated conditional mean of \code{Qj.star.t+1}
#'         evaluated at the each of the rows of each \code{data.frame} in
#'         \code{wideDataList}.
#'

estimateIteratedMeanT <- function(wideDataList, t, whichJ, allJ, t0, adjustVars, trt,
                                 SL.ftime = NULL, glm.ftime = NULL, verbose,
                                 returnModels = FALSE, bounds = NULL,
                                 trtofTime=NULL, #### new
                                 trtOfInterest, ...){
  ## determine who to include in estimation
  include <- rep(TRUE, nrow(wideDataList[[1]]))
  n_regimen <-ncol(trtOfInterest)-1
  regimen <- trtOfInterest[,-1]

  if(t != 1) {
    for(j in allJ) {
      # exclude previously failed subjects
      include[wideDataList[[1]][[paste0("N",j,".",t-1)]]==1] <- FALSE
    }
  }
  # exclude previously censored subjects
  include[wideDataList[[1]][[paste0("C.",t)]]==1] <- FALSE


  ## determine the outcome for the regression
  outcomeName <- ifelse(t == t0, paste0("N", whichJ, ".", t0),
                        paste0("Q", whichJ, ".", t + 1, ".", t0, ".0"))

  ## create an indicator of any failure prior to t
  wideDataList <- lapply(wideDataList, function(x, t){
    if(length(allJ) > 1) {
      x[[paste0("NnotJ.",t - 1)]] <-
        rowSums(cbind(rep(0, nrow(x)), x[, paste0('N', allJ[allJ != whichJ], '.', t - 1)]))
    } else {
      x[[paste0("NnotJ.",t - 1)]] <- 0
    }
    x
  }, t = t)

  lj.t <- paste0("l", whichJ, ".", t)
  uj.t <- paste0("u", whichJ, ".", t)
  Qtildej.t <- paste0("Qtilde", whichJ, ".", t)
  Nj.tm1 <- paste0("N", whichJ, ".", t - 1)
  Qj.t <- paste0("Q", whichJ, ".", t, ".", t0, ".0")
  NnotJ.tm1 <- paste0("NnotJ.", t - 1)
  if (is.list(glm.ftime)){
    if ( length(glm.ftime)  >= max(t0)){
      Qform <- paste(outcomeName, "~", glm.ftime[[paste0("t", t)]], sep = " ")
    } else{
        stop("need to specified regression form for censoring at all the time points")
    }
  } else {
      Qform <- paste(outcomeName, "~", glm.ftime, sep = " ")
  }


  #### creating matrix for past covariates t-1 and treatment up to t-1

  trt_include <-c(paste0("trt_t", trtofTime[trtofTime < t]))
  trt_past <- trt[ , trt_include, drop = F]


  var_include <- NULL
  for (var.ind in trtofTime[trtofTime < t]){
    var_include <- c(var_include,colnames(adjustVars)[grepl(paste0("t",var.ind),
                                                            colnames(adjustVars))])
  }

  var_past <- adjustVars[ , c(var_include)]

  trt_var_past <- cbind(trt_past, var_past)
  ### check if it lines up??


  ## GLM code
  if(is.null(SL.ftime)) {
    if(is.null(bounds)) { # with no bounds


      if (t == t0) {
        glmdata <- cbind(wideDataList[[1]][include, outcomeName, drop =F],
                         trt_var_past[include, ])

      } else {
        glmdata <- vector("list", ncol(trtOfInterest) - 1)
        for (r in seq(ncol(trtOfInterest) - 1)){
          trt_data_temp <- trt_var_past
          #if (t -1 %in% trtofTime){
          #    trt_data_temp[, paste0("trt_t",trtofTime[t -1])] <-
          #      regimen[t -1,r]
          #}
          ## double check the outcome name is correct

          outcome_r <- paste0("Q", whichJ, ".", t+1, ".", t0, ".", r)
          glmdata[[r]] <- cbind(wideDataList[[1]][include, outcome_r, drop =F],
                                trt_data_temp[include, ])
          #print(paste0("t: ", t,"t0: ", t0,"J: ", whichJ, "r: ", r))
          #print(mean(wideDataList[[1]][include, outcome_r]))
          if (t <= max(trtofTime)){
            glmdata[[r]][["r.ind"]] = r
          }
          colnames(glmdata[[r]])[which( colnames(glmdata[[r]]) ==  outcome_r)] <- outcomeName
        }

        glmdata <- do.call(rbind, glmdata)
        if (t <= max(trtofTime)){
          glmdata$r.ind = factor(glmdata$r.ind)
        }
      }

      suppressWarnings({
        Qmod <- fast_glm(reg_form = stats::as.formula(Qform),
                         data = glmdata,
                         family = stats::binomial())
        if (unique(class(Qmod) %in% c("glm", "lm"))) {
          Qmod <- cleanglm(Qmod)
        }

        wideDataList <- lapply(wideDataList, function(x, whichJ, t, Qj.t, Nj.tm1, NnotJ.tm1) {
          #pred.temp <- suppressWarnings(
          ###### predict Q bar under observed treatment
          #     predict(Qmod,newdata=trt_var_past,type="response")
          #)
          ###### setting all the missing value to 1
          #pred.temp[is.na(pred.temp)] <- 1
          #x[[Qj.t]] <- x[[Nj.tm1]] + (1-x[[NnotJ.tm1]]-x[[Nj.tm1]])*pred.temp

          ###### predict Q bar under each treach regrimen
          for(regimen_ind in 1:n_regimen){
            Qj.t.r <- paste0("Q", whichJ, ".", t, ".", t0, ".", regimen_ind)
            trt_data_pred <- trt_var_past
            if ((t -1) %in% trtofTime){
              trt_data_pred[, paste0("trt_t",t -1)] <-
                regimen[which(trtofTime == t -1),regimen_ind]
            }
            if (t <= max(trtofTime)){
              trt_data_pred[["r.ind"]] = factor(regimen_ind, levels = seq_len(n_regimen))
            }

            pred.temp.r <- suppressWarnings(
              predict(Qmod,newdata=trt_data_pred,type="response")
            )
            # Maybe an issue which is the same as the :pred.temp.r[is.na(pred.temp.r)] <- 1
            pred.temp.r[x[[Nj.tm1]] == 1] <- 1
            pred.temp.r[x[[NnotJ.tm1]] == 1] <- 1
            #pred.temp.r[is.na(pred.temp.r)] <- 1

          #  past_t <- trtofTime[trtofTime < t]
          #  ind.r.t_minus1 <-  apply(x[,paste0("trt_t", past_t), drop = FALSE], 1, function(x){
          #    as.numeric(all(x == trtOfInterest[(trtofTime < t), paste0("regimen", regimen_ind)]))  } )
          #  ind.r.t_minus1[is.na(ind.r.t_minus1)] <- 1
              x[[Qj.t.r]] <- x[[Nj.tm1]] + (1-x[[NnotJ.tm1]]-x[[Nj.tm1]])*pred.temp.r
          }

          x
        }, t = t, whichJ = whichJ, Qj.t =  Qj.t, Nj.tm1 = Nj.tm1, NnotJ.tm1 =  NnotJ.tm1)
      })
    } else { # with bounds
      X <- stats::model.matrix(stats::as.formula(Qform),
                               data = cbind(wideDataList[[1]][include, outcomeName, drop =F]
                                            , trt_var_past[include, ]))
      Ytilde <- (wideDataList[[1]][include, outcomeName] -
                   wideDataList[[1]][[lj.t]][include]) /
        (wideDataList[[1]][[uj.t]][include] -
           wideDataList[[1]][[lj.t]][include])
      Qmod <- stats::optim(par = rep(0, ncol(X)), fn = LogLikelihood,
                           Y = Ytilde, X = X, method = "BFGS", gr = grad,
                           control = list(reltol = 1e-7, maxit = 50000))
      beta <- Qmod$par


      #### need to figure out what to do with the missing value in the trt_var_past
      wideDataList <- lapply(wideDataList, function(x, j, t) {
        newX <- stats::model.matrix(stats::as.formula(Qform), data = trt_var_past )
        x[[Qj.t]] <- x[[Nj.tm1]] + (1 - x[[NnotJ.tm1]] - x[[Nj.tm1]]) *
          (plogis(newX %*% beta) * (x[[uj.t]] - x[[lj.t]]) + x[[lj.t]])
        x
      }, j = whichJ, t = t)
    }
  } else if(is.null(glm.ftime)) { # Super Learner
    if(is.null(bounds)) { # with no bounds
      # some stability checks
      # number of unique outcome values
      nUniq <- length(unique(wideDataList[[1]][include, outcomeName]))
      cvControl <- SuperLearner::SuperLearner.CV.control()
      if(t == t0) {
        # if there are less than 2 events at t0, just fit regression using only Z
        nE <- sum(wideDataList[[1]][include, outcomeName])
        ignoreSL <- nE <= 2
        if(ignoreSL) {
          suppressWarnings({
            Qform_trt <- paste(outcomeName, "~", "trt", sep = " ")
            Qmod <- fast_glm(reg_form = stats::as.formula(Qform_trt),
                             data = cbind(wideDataList[[1]][include, outcomeName, drop =F],
                                          trt_var_past[include, ]),
                             family = stats::gaussian())
            wideDataList <- lapply(wideDataList, function(x, whichJ, t, Qj.t, Nj.tm1, NnotJ.tm1) {
              suppressWarnings(
                x[[Qj.t]] <- x[[Nj.tm1]] + (1-x[[NnotJ.tm1]]- x[[Nj.tm1]]) *
                  predict(Qmod,newdata=trt_var_past)
              )


              ###### predict Q bar under each treach regrimen
              for(regimen_ind in 1:n_regimen){
                Qj.t.r <- paste0("Q", whichJ, ".", t, ".", t0, ".", regimen_ind)
                trt_data_pred <- trt_var_past
                for(regimen_time in 1:length(trtofTime[which(trtofTime < t)])){
                  trt_data_pred[,paste0("trt_t",trtofTime[regime_time])] <- regimen[regimen_time,regimen_ind]
                }

                # past_t <- trtofTime[trtofTime < t]
                # ind.r.t_minus1 <-  apply(x[,paste0("trt_t", past_t), drop = FALSE], 1, function(x){
                #   as.numeric(all(x == trtOfInterest[(trtofTime < t), paste0("regimen", regimen_ind)]))  } )
                # ind.r.t_minus1[is.na(ind.r.t_minus1)] <- 1

                suppressWarnings(
                  x[[Qj.t.r]] <- x[[Nj.tm1]] + (1-x[[NnotJ.tm1]]-x[[Nj.tm1]])*
                    predict(Qmod,newdata=trt_data_pred,type="response")
                )
              }

              x
            }, t = t, whichJ = whichJ)
          })
        } else {
          simplify <- nE <= cvControl$V
          if(simplify) cvControl <- list(V = nE - 1, stratifyCV = TRUE)
          suppressWarnings(
            Qmod <- SuperLearner::SuperLearner(Y = wideDataList[[1]][include,outcomeName],
                                               X = trt_var_past[include, ],
                                               SL.library = SL.ftime,
                                               cvControl = list(V = 5), method = "method.CC_nloglik",

                                               family = "binomial",
                                               verbose = verbose)
          )
          wideDataList <- lapply(wideDataList, function(x, whichJ, t,  Qj.t, Nj.tm1, NnotJ.tm1) {
            x[[Qj.t]] <- x[[Nj.tm1]] + (1-x[[NnotJ.tm1]]-x[[Nj.tm1]])*
              predict(Qmod, newdata = trt_var_past, onlySL = TRUE)$pred
            x
          }, t = t, whichJ = whichJ, Qj.t =  Qj.t, Nj.tm1 = Nj.tm1, NnotJ.tm1 =  NnotJ.tm1)
          }
      } else {
        suppressWarnings(
          Qmod <- SuperLearner::SuperLearner(Y = wideDataList[[1]][include, outcomeName],
                                             X = cbind(trt_var_past[include, ]),
                                             SL.library = SL.ftime,
                                             cvControl = list(V = 5),
					     method = "method.CC_LS",
                                             family = "binomial",
                                             verbose = verbose)
        )
        wideDataList <- lapply(wideDataList, function(x, whichJ, t, Qj.t, Nj.tm1, NnotJ.tm1) {
          suppressWarnings(
            x[[Qj.t]] <- x[[Nj.tm1]] + (1 - x[[Nj.tm1]] - x[[NnotJ.tm1]]) *
              predict(Qmod, newdata = trt_var_past, onlySL = TRUE)$pred
          )

          ###### predict Q bar under each treach regrimen
          for(regimen_ind in 1:n_regimen){
            Qj.t.r <- paste0("Q", whichJ, ".", t, ".", t0, ".", regimen_ind)
            trt_data_pred <- trt_var_past
            for(regimen_time in 1:length(trtofTime[which(trtofTime <= t)])){
              trt_data_pred[,paste0("trt_t",trtofTime[regime_time])] <- regimen[regimen_time,regimen_ind]
            }
            suppressWarnings(
              x[[Qj.t.r]] <- x[[Nj.tm1]] + (1-x[[NnotJ.tm1]]-x[[Nj.tm1]])*
                predict(Qmod,newdata=trt_data_pred,type="response")
            )
          }


          x
        }, t = t, whichJ = whichJ, Qj.t =  Qj.t, Nj.tm1 = Nj.tm1, NnotJ.tm1 =  NnotJ.tm1)
        }
    } else {
      stop("Super Learner code with bounds not written yet")
    }
  }
  out <- list(wideDataList = wideDataList,
              ftimeMod = if(returnModels == TRUE) {
                Qmod
              } else {
                NULL
              }
  )
  return(out)
}
