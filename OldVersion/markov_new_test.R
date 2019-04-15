rm(list = ls())

library(etm)
library(mstate)

#markovian case

BSiter <- 500
set.seed(1000)
d <- 1
bs_one_iter <- function(n, data, max_id, times=x, h, tra){
  
  multinom <- sample(1:n, replace=TRUE)
  
  data_bs <- subset(data, id==data$id[multinom[1]])
  for(k in 2:n){
    data_temp <- subset(data, id==data$id[multinom[k]])
    data_temp$id <- data_temp$id + max_id * k
    data_bs <- rbind(data_bs, data_temp)
  }
  save.image()
  etmAJE <- etm(data_bs, state.names = 1:3, tra=tra, cens.name=4, s=s, t=t, covariance=FALSE, delta.na=FALSE)
  PijAJE <- c(0,etmAJE$est[i ,j ,])
  timesAJE <- c(s,etmAJE$time)
  
  good_ind <- which(data$entry <= s & data$exit >s & data$from==i)
  rel_ids <- unique(data$id[good_ind])
  npeople <- length(rel_ids)
  data_s <- data[data$id %in% rel_ids,]
  #Remove all intervals ending before s
  data_s <- data_s[data_s$exit >= s,]
  #Truncate all intervals starting before s
  data_s$entry[data_s$entry <s]<-s
  message(paste0("+++ Laufindex innerhalb der Fkt c jetzt: ", c))
  save.image()
  etm_s <- etm(data_s, state.names = 1:3, tra=tra, cens.name=4, s=s, t=t, covariance=FALSE, delta.na=FALSE)
  c <- c+1
  return(data_bs)
}

for (k in 1:1000){

  # simulation
  n <- 100
  rate_0 <- 1/30
  start_state <- rep(1,n)
  T_0 <- rexp(n,rate_0) # time in state 0
  p <- 0.7
  r <- runif(n,0,1)
  state_1 <- 2*(r < p) + 3*(p <= r & r < p+(1-p)/2) + 4*(1-(1-p)/2 <= r)
  rate_1 <- 1/60
  T_1 <- rexp(n,rate_1) # time in state 1
  q <- 0.8
  state_2 <- rbinom(n,1,1-q) + 3
  # results in time-continuous markovian multi-state process'

  # state 3 = censored
  
  # next we transform the generated data into a data matrix that shows transition times and states per patient
  
  all_datas <- matrix(rep(0,5*2*n),nrow = 2*n)
  count <- 1
  for (i in 1:n){
    if(state_1[i] == 4 | state_1[i] == 3){
      all_datas[count,] <- c(i,1,state_1[i],0,T_0[i])
    }
    else{
      all_datas[count,] <- c(i,1,2,0,T_0[i])
      count <- count + 1
      all_datas[count,] <- c(i,2,state_2[i],T_0[i],T_0[i] + T_1[i])
    }
    count <- count + 1
  }
  
  datas <- data.frame(all_datas[-which(all_datas[,1] == 0),])
  names(datas) <- c("id","from","to","entry","exit")
  
  tra <- matrix(c(FALSE,TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE),3,3,byrow=TRUE)
  
  # from illness to death
  i <- 2
  j <- 3
  # time stamp
  s <- 10
  t <- 1200
  h <- 0.5
  x_markov <- seq(s,t,h)
  
  # estimate transition probabilities
  message(paste0("+++ Laufindex außerhalb der Fkt d jetzt: ", d))  
  save.image()

  etmAJE_markov <- etm(datas, state.names = 1:3, tra=tra, cens.name=4, s=s, t=t, covariance=FALSE, delta.na=FALSE)
  PijAJE_markov <- c(0,etmAJE_markov$est[i ,j ,])
  timesAJE_markov <- c(s,etmAJE_markov$time)
  
  good_ind_markov <- which(datas$entry <= s & datas$exit >s & datas$from == i)
  rel_ids_markov <- unique(datas$id[good_ind_markov])
  npeople_markov <- length(rel_ids_markov)
  # "good" data
  datas_s <- data.frame(datas[datas$id %in% rel_ids_markov,])
  # Remove all intervals ending before s
  datas_s <- datas_s[datas_s$exit >= s,]
  #Truncate all intervals starting before s
  datas_s$entry[datas_s$entry < s] <- s
  message(paste0("+++ Laufindex außerhalb der Fkt d jetzt: ", d))
  etm_s_markov <- etm(datas_s, state.names = 1:3, tra=tra, cens.name=4, s=s, t=t, covariance=FALSE, delta.na=FALSE)
  F0_markov <- c(0,etm_s_markov$est[i ,j ,])
  timesTE_markov <- c(s,etm_s_markov$time)
  
  # Computing the Aalen Johansen and Titman estimators


  # Computing the test statistic
  c <- d
  
  max_id_markov <- max(datas$id)
  # BSstat_pred <- replicate(BSiter, bs_one_iter(n, datas, max_id_markov, times=x_markov, h, tra))
  d <- d + 1
  
}

