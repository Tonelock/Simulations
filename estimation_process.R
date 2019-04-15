bs_one_iter <- function(n, data, max_id, AJE, TE, times=x, h, tra){
  
  multinom <- sample(1:n, replace = TRUE)
  data_bs <- subset(data, id==data$id[multinom[1]])
  for(l in 2:n){
    data_temp <- subset(data, id==data$id[multinom[l]])
    data_temp$id <- data_temp$id + max_id * l
    data_bs <- rbind(data_bs, data_temp)
  }
  etmAJE <- etm(data_bs, state.names = 0:2, tra=tra, cens.name=3, s=s, t=t, covariance=FALSE, delta.na=FALSE)

  PijAJE <- c(0,etmAJE$est[i + 1,j + 1,])

  timesAJE <- c(s,etmAJE$time)

  good_ind <- which(data$entry < s & data$exit >s & data$from==i)

  rel_ids <- unique(data$id[good_ind])

  npeople <- length(rel_ids)
  data_s <- data[data$id %in% rel_ids,]
  #Remove all intervals ending before s
  data_s <- data_s[data_s$exit > s,]
  #Truncate all intervals starting before s
  data_s$entry[data_s$entry <s]<-s
  
  etm_s <- etm(data_s, state.names = 0:2, tra=tra, cens.name=3, s=s, t=t, covariance=FALSE, delta.na=FALSE)
  F0 <- c(0,etm_s$est[i + 1,j + 1,])
  timesTE <- c(s,etm_s$time)
  
  TE_bs <- approxfun(timesTE, F0, method="constant")
  AJE_bs <- approxfun(timesAJE, PijAJE, method="constant")
# !
  TE_bs <- ifelse(is.na(TE_bs(times)), max(TE_bs(times), na.rm=TRUE), TE_bs(times))

  AJE_bs <- ifelse(is.na(AJE_bs(times)), max(AJE_bs(times), na.rm=TRUE), AJE_bs(times))

  return(n*sum((AJE_bs-AJE-TE_bs+TE)^2 * h))
}


estimation <- function(data,s,t,i,j,BSiter,tra){
  # estimate transition probabilities
  etmAJE <- etm(data, state.names = 0:2, tra=tra, cens.name=3, s=s, t=t, covariance=FALSE, delta.na=FALSE)
  PijAJE <- c(0,etmAJE$est[i + 1,j + 1,])
  timesAJE <- c(s,etmAJE$time)
  
  good_ind <- which(data$entry < s & data$exit > s & data$from == i)
  rel_ids <- unique(data$id[good_ind])
  npeople <- length(rel_ids)
  # "good" data
  data_s <- data[data$id %in% rel_ids,]
  # Remove all intervals ending before s
  data_s <- data_s[data_s$exit > s,]
  #Truncate all intervals starting before s
  data_s$entry[data_s$entry < s] <- s
  
  etm_s <- etm(data_s, state.names = 0:2, tra=tra, cens.name=3, s=s, t=t, covariance=FALSE, delta.na=FALSE)
  F0 <- c(0,etm_s$est[i + 1,j + 1,])
  timesTE <- c(s,etm_s$time)
  
  h <- 0.5
  x <- seq(s,t,h)
  TE <- approxfun(timesTE, F0, method="constant")  
  AJE <- approxfun(timesAJE, PijAJE, method="constant")
  plot(x, ifelse(is.na(TE(x)), max(TE(x), na.rm=TRUE), TE(x)), type="l", ylim=c(0,1), xlab="time (days)", ylab="prob critical to dead")
  lines(x, ifelse(is.na(AJE(x)), max(AJE(x), na.rm=TRUE), AJE(x)), type="l", col=2)
  legend("bottomright", legend=c("Titman estimator", "Aalen-Johansen estimator"), col=1:2, lty=1)
  
  TE_pred <- ifelse(is.na(TE(x)), max(TE(x), na.rm=TRUE), TE(x))
  AJE_pred <- ifelse(is.na(AJE(x)), max(AJE(x), na.rm=TRUE), AJE(x))
  
  teststat_pred <- n*sum((AJE_pred-TE_pred)^2 * h)
  
  max_id <- max(data$id)
  
  BSstat_pred <- replicate(BSiter, bs_one_iter(n, data, max_id, AJE_pred, TE_pred, times=x, h, tra))
  
  return(sum(BSstat_pred >= teststat_pred)/BSiter)
  
}
