rm(list = ls())

library(etm)
library(mstate)

#markovian case

alpha <- 0.05
BSiter <- 300

conf_level <- rep(0,1000)
test_accepted <- rep(0,1000)

for (k in 1:1000){

  # simulation
  n <- 800
  rate_0 <- 1/100
  T_01 <- rexp(n,rate_0) # transition times from healthy to sick state
  rate_1 <- 1/250
  T_12 <- T_01 + rexp(n,rate_1) # transition times from sick to death state
  cens <- rexp(n,1/500) # censoring times
  T_02 <- rexp(n,1/250) # transition times from healthy to death state
  # we will consider the first transition (minimum of random numbers T_01,T_02,cens) as the transition that actually happens, same for the data, if we are in state 1 (minimum of T_12, cens)
  # results in time-continuous markovian multi-state process'
  id <- 1:n
  start <- 20
  
  #state 3 = censored
  
  # next we transform the generated data into a data matrix that shows transition times and states per patient
  entry_state <-  (cens > start) * {0*(start < pmin(T_01,T_02)) + 1*((T_01 < pmin(T_02,start)) & (start < T_12)) + 2*((T_02 < start) | (T_12 < start))} + 
            3*(cens < start)
  first_data <- pmin(cens,T_01,T_02)
  all_datas <- matrix(rep(0,5*4*n),nrow = 4*n)
  count <- 1
  for (i in id){
    if (entry_state[i] == 0){
      
      if (T_01[i] < cens[i]){
        all_datas[count,] <- c(i,0,1,0,T_01[i])    
        count <- count + 1
        
        if(T_12[i] < cens[i]){
          all_datas[count,] <- c(i,1,2,T_01[i],T_12[i])
        }
        else{
          all_datas[count,] <- c(i,1,3,T_01[i],cens[i])
        }
        
      }
      else{
        all_datas[count,] <- c(i,0,3,0,cens[i])
      }
    }
    else if(entry_state[i] == 1){
      
      if(T_12[i] < cens[i]){
        all_datas[count,] <- c(i,1,2,0,T_12[i])
      }
      else{
        all_datas[count + 1,] <- c(i,1,3,0,cens[i])
      }
      
    }
    
    count <- count + 1
  }
  
  datas <- data.frame(all_datas[-which(all_datas[,1] == 0),])
  names(datas) <- c("id","from","to","entry","exit")
  
  tra <- matrix(c(F,T,T,F,F,T,F,F,F),3,3,byrow=T)
  
  # from illness to death
  i <- 1
  j <- 2
  # time stamp
  s <- 100
  t <- 1200
  
  # estimate transition probabilities
  etmAJE_markov <- etm(datas, state.names = 0:2, tra=tra, cens.name=3, s=s, t=t, covariance=FALSE, delta.na=FALSE)
  PijAJE_markov <- c(0,etmAJE_markov$est[i + 1,j + 1,])
  timesAJE_markov <- c(s,etmAJE_markov$time)
  
  good_ind_markov <- which(datas$entry <= s & datas$exit >s & datas$from == i)
  rel_ids_markov <- unique(datas$id[good_ind_markov])
  npeople_markov <- length(rel_ids_markov)
  # "good" data
  datas_s <- datas[datas$id %in% rel_ids_markov,]
  # Remove all intervals ending before s
  datas_s <- datas_s[datas_s$exit >= s,]
  #Truncate all intervals starting before s
  datas_s$entry[datas_s$entry < s] <- s
  
  etm_s_markov <- etm(datas_s, state.names = 0:2, tra=tra, cens.name=3, s=s, t=t, covariance=FALSE, delta.na=FALSE)
  F0_markov <- c(0,etm_s_markov$est[i + 1,j + 1,])
  timesTE_markov <- c(s,etm_s_markov$time)
  
  h <- 0.5
  x_markov <- seq(s,t,h)
  TE_markov <- approxfun(timesTE_markov, F0_markov, method="constant")
  plot(x_markov, ifelse(is.na(TE_markov(x_markov)), max(TE_markov(x_markov), na.rm=T), TE_markov(x_markov)), type="l", ylim=c(0,1), xlab="time (days)", ylab="prob critical to dead")
  AJE_markov <- approxfun(timesAJE_markov, PijAJE_markov, method="constant")
  lines(x_markov, ifelse(is.na(AJE_markov(x_markov)), max(AJE_markov(x_markov), na.rm=T), AJE_markov(x_markov)), type="l", col=2)
  legend("bottomright", legend=c("Titman estimator", "Aalen-Johansen estimator"), col=1:2, lty=1)
  
  TE_pred_markov <- ifelse(is.na(TE_markov(x_markov)), max(TE_markov(x_markov), na.rm=T), TE_markov(x_markov))
  AJE_pred_markov <- ifelse(is.na(AJE_markov(x_markov)), max(AJE_markov(x_markov), na.rm=T), AJE_markov(x_markov))
  
  teststat_pred_markov <- n*sum((AJE_pred_markov-TE_pred_markov)^2 * h)
  
  bs_one_iter <- function(n, data, max_id, AJE, TE, times=x, h, tra){
    
    multinom <- sample(1:n, replace=T)
    
    data_bs <- subset(data, id==data$id[multinom[1]])
    for(k in 2:n){
      data_temp <- subset(data, id==data$id[multinom[k]])
      data_temp$id <- data_temp$id + max_id * k
      data_bs <- rbind(data_bs, data_temp)
    }
    
    etmAJE <- etm(data_bs, state.names = 0:2, tra=tra, cens.name=3, s=s, t=t, covariance=FALSE, delta.na=FALSE)
    PijAJE <- c(0,etmAJE$est[i + 1,j + 1,])
    timesAJE <- c(s,etmAJE$time)
    
    good_ind <- which(data$entry <= s & data$exit >s & data$from==i)
    rel_ids <- unique(data$id[good_ind])
    npeople <- length(rel_ids)
    data_s <- data[data$id %in% rel_ids,]
    #Remove all intervals ending before s
    data_s <- data_s[data_s$exit >= s,]
    #Truncate all intervals starting before s
    data_s$entry[data_s$entry <s]<-s
    
    etm_s <- etm(data_s, state.names = 0:2, tra=tra, cens.name=3, s=s, t=t, covariance=FALSE, delta.na=FALSE)
    F0 <- c(0,etm_s$est[i + 1,j + 1,])
    timesTE <- c(s,etm_s$time)
    
    TE_bs <- approxfun(timesTE, F0, method="constant")
    AJE_bs <- approxfun(timesAJE, PijAJE, method="constant")
    
    TE_bs <- ifelse(is.na(TE_bs(times)), max(TE_bs(times), na.rm=T), TE_bs(times))
    AJE_bs <- ifelse(is.na(AJE_bs(times)), max(AJE_bs(times), na.rm=T), AJE_bs(times))
    return(n*sum((AJE_bs-AJE-TE_bs+TE)^2 * h))
  }
  
  
  max_id_markov <- max(datas$id)
  #set.seed(1234)
  BSstat_pred <- replicate(BSiter, bs_one_iter(n, datas, max_id_markov, AJE_pred_markov, TE_pred_markov, times=x_markov, h, tra))
  sum(BSstat_pred >= teststat_pred_markov)/BSiter
  # [1] 0.967
  
  BSstat_pred_sort <- sort(BSstat_pred)
  Cn <- BSstat_pred_sort[BSiter*(1-alpha)]
  
  conf_level[k] <- sum(BSstat_pred >= teststat_pred_markov)/BSiter
  # [1] 0.967
  
  BSstat_pred_sort <- sort(BSstat_pred)
  Cn <- BSstat_pred_sort[BSiter*(1-alpha)]
  test_accepted[k] <- (teststat_pred_markov <= Cn) #accept H0 with conf level alpha?
  
  save.image()
  
  
}
