rm(list = ls())

library(etm)
library(mstate)
source("estimation_process.R")
# continuous case: First non-markovian, then markovian
set.seed(1234)

alpha <- 0.05
BSiter <- 20 # maybe more ~2000 for final simus 
number_iteration <- 50 # more for final simus ~5000

rate_0 <- 10
rate_1 <- 10
rate_2 <- 5
rate_cens <- 18

rate_21 <- 10
rate_22 <- 5

conf_level <- rep(0,number_iteration)
test_rejected <- rep(0,number_iteration)

tra <- matrix(c(FALSE,TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE),3,3,byrow=TRUE)

# from illness to death
i <- 1
j <- 2
# time stamp
s <- 12
t <- 30

for (n in c(1000)){
  
  for (k in 1:number_iteration){
    
    # simulation
    
    start_state <- rep(0,n)
    T_01 <- rpois(n,rate_0) + 1 # time health to ill state
    T_02 <- rpois(n,rate_1) + 1 # time health to death state 
    T_12 <- rpois(n,rate_21)*(pmin(T_01,T_02) < 17) + rpois(n,rate_22)*(pmin(T_01,T_02) >= 17) # time from ill to death state, dependend on the state of each individual at time 25
    cens <- rpois(n,rate_cens) # time to censoring
    state_1 <- 1*(T_01 < pmin(T_02,cens)) + 2*(T_02 < pmin(T_01,cens)) + 3*(cens <= pmin(T_01,T_02))
    state_2 <- 2*(cens >= T_01+T_12) + 3*(T_01 < cens & cens < T_01+T_12)
    
    # results in time-discrete non-markovian multi-state process
    # state 3 = censored
    
    # next we transform the generated data into a data matrix that shows transition times and states per patient
    
    all_data <- matrix(rep(0,5*2*n),nrow = 2*n)
    count <- 1
    for (m in 1:n){
      if(state_1[m] == 2){
        all_data[count,] <- c(m,0,state_1[m],0,T_02[m])
      }
      else if(state_1[m] == 3){
        all_data[count,] <- c(m,0,state_1[m],0,cens[m])
      }
      else if(state_1[m] == 1 & state_2[m] == 2){
        all_data[count,] <- c(m,0,1,0,T_01[m])
        count <- count + 1
        all_data[count,] <- c(m,1,2,T_01[m],T_01[m] + T_12[m])
      }
      else if(state_1[m] == 1 & state_2[m] == 3){
        all_data[count,] <- c(m,0,1,0,T_01[m])
        count <- count + 1
        all_data[count,] <- c(m,1,3,T_01[m],cens[m])
      }
      count <- count + 1
    }
    if (all_data[2*n,1] == 0){data <- data.frame(all_data[-which(all_data[,1] == 0),])} else {data == all_data}
    names(data) <- c("id","from","to","entry","exit")
    
    conf_level[k] <- estimation(data,s,t,i,j,BSiter,tra)
    test_rejected[k] <- (conf_level[k] <= alpha) 
    
    # BSstat_pred_sort <- sort(BSstat_pred)
    # Cn <- BSstat_pred_sort[BSiter*(1-alpha)]
    # test_accepted[k] <- (teststat_pred_markov <= Cn) #accept H0 with conf level alpha?
  }
  string <- paste("nonmarkovdisc",number_iteration,"n",n,".RData",sep = "")
  save.image(string)
}

# next: Markovian case

for (n in c(1000)){
  
  for (k in 1:number_iteration){
  
    # simulation
  
    start_state <- rep(0,n)
    T_01 <- rpois(n,rate_0) + 1 # time health to ill state
    T_02 <- rpois(n,rate_1) + 1 # time health to death state 
    T_12 <- rpois(n,rate_2) + 1 # time from ill to death state
    cens <- rpois(n,rate_cens) + 1 # time to censoring
    state_1 <- 1*(T_01 < pmin(T_02,cens)) + 2*(T_02 < pmin(T_01,cens)) + 3*(cens <= pmin(T_01,T_02))
    state_2 <- 2*(cens >= T_01+T_12) + 3*(T_01 < cens & cens < T_01+T_12)
    
    # results in time-discrete markovian multi-state process'
    # state 3 = censored
    
    # next we transform the generated data into a data matrix that shows transition times and states per patient
    
    all_data <- matrix(rep(0,5*2*n),nrow = 2*n)
    count <- 1
    for (m in 1:n){
      if(state_1[m] == 2){
        all_data[count,] <- c(m,0,state_1[m],0,T_02[m])
      }
      else if(state_1[m] == 3){
        all_data[count,] <- c(m,0,state_1[m],0,cens[m])
      }
      else if(state_1[m] == 1 & state_2[m] == 2){
        all_data[count,] <- c(m,0,1,0,T_01[m])
        count <- count + 1
        all_data[count,] <- c(m,1,2,T_01[m],T_01[m] + T_12[m])
      }
      else if(state_1[m] == 1 & state_2[m] == 3){
        all_data[count,] <- c(m,0,1,0,T_01[m])
        count <- count + 1
        all_data[count,] <- c(m,1,3,T_01[m],cens[m])
      }
      count <- count + 1
    }
    
    if (all_data[2*n,1] == 0){data <- data.frame(all_data[-which(all_data[,1] == 0),])} else {data == all_data}
    names(data) <- c("id","from","to","entry","exit")
    
    conf_level[k] <- estimation(data,s,t,i,j,BSiter,tra)
    test_rejected[k] <- (conf_level[k] <= alpha) 
    # BSstat_pred_sort <- sort(BSstat_pred)
    # Cn <- BSstat_pred_sort[BSiter*(1-alpha)]
    # test_accepted[k] <- (teststat_pred_markov <= Cn) #accept H0 with conf level alpha?
    }
  string <- paste("markovdisc",number_iteration,"n",n,".RData",sep = "")
  save.image(string)
}