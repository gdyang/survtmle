### create TVC treatment and covariates

trtDataframe<-  function(id = id, trt = trt, trtofTime = trtofTime,...){
  
  trtData <- data.frame(id = id)
  
  trt_hist.temp <- lapply(trt[paste0("t", trtofTime)], function(trtsub){
    trt.temp <- merge(data.frame(id = id), trtsub, by = "id", all.x = T)
    return(trt.temp)
  })
  
  for(list.ind in 1: length(trt_hist.temp)){
    trt_t <- names(trt_hist.temp[list.ind])
    temp <- trt_hist.temp[[list.ind]]
    colnames(temp)[which(colnames(temp) != "id")] <- paste0("trt_", trt_t)
    trtData <- merge(trtData, temp, by = "id")
  } 


return(trtData)
}


VarDataframe <-  function(id = id, adjustVars = adjustVars, trtofTime = trtofTime,...){

  
  var <- data.frame(id = id)
  
  

    var_past.temp <- lapply(adjustVars[paste0("t", trtofTime)], function(varsub){
      var.temp <- merge(data.frame(id = id), varsub, by = "id", all.x = T)
      return(var.temp)
    })
    
    for(list.ind in 1: length(var_past.temp)){
      var_t <- names(var_past.temp[list.ind])
      temp <- var_past.temp[[list.ind]]
      colnames(temp)[which(colnames(temp) != "id")] <- paste0(colnames(temp[-1]), "_", var_t)
      var <- merge( var, temp, by = "id")
    }  
  
  
  
  return(var)
}





wideDataT <- function(dat = dat, allJ = allJ, ...){

J <- allJ

n <- nrow(dat)
rankftime <- match(dat$ftime, sort(unique(dat$ftime)))

# first element used for estimation
out <- dat[rep(1:nrow(dat), rankftime), ]

out <- dat[rep(1:nrow(dat), rankftime), ]

for(j in J) {
  out[[paste0("N",j)]] <- 0
  out[[paste0("N",j)]][cumsum(rankftime)] <- as.numeric(dat$ftype==j)
}
out$C <- 0
out$C[cumsum(rankftime)] <- as.numeric(dat$ftype == 0)

n.row.ii <- nrow(out)
uniqftime <- unique(dat$ftime)
orduniqftime <- uniqftime[order(uniqftime)]
row.names(out)[row.names(out) %in%
                           paste(row.names(dat))] <- paste0(row.names(dat),".0")
out$t <- orduniqftime[as.numeric(paste(unlist(strsplit(row.names(out), ".",
                                                                 fixed = TRUE))[seq(2, n.row.ii * 2, 2)])) + 1]



if(!is.null(bounds)){
  boundFormat <- data.frame(t = bounds$t)
  for(j in J){
    if(paste("l", j, sep = "") %in% colnames(bounds)) {
      boundFormat[[paste0("l",j)]] <- bounds[,paste0("l",j)]
    } else {
      boundFormat[[paste0("l",j)]] <- 0
    }
    if(paste("u", j, sep = "") %in% names(bounds)) {
      boundFormat[[paste0("u",j)]] <- bounds[,paste0("u",j)]
    } else {
      boundFormat[[paste0("u",j)]] <- 1
    }
  }
  suppressMessages(
    out <- plyr::join(x = out, y = boundFormat,
                                type = "left")
  )
  # if any bounds are missing, add in 0 and 1
  for(j in J) {
    tmp <- is.na(out[, paste0("l", j)])
    out[tmp, paste0("l", j)] <- 0
    tmp <- is.na(out[, paste0("u", j)])
    out[tmp,paste0("u", j)] <- 1
  }
} else {
  for(j in J){
    out[[paste0("l",j)]] <- 0
    out[[paste0("u",j)]] <- 1
  }
}

dropVars <- c("ftime","ftype")
out <- data.frame(stats::reshape(out[, -which(colnames(out) %in% dropVars) ],
                                               direction = "wide",
                                               timevar = "t", idvar = "id"))

return(out)
}

