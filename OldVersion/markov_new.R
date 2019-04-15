rm(list = ls())

library(etm)
library(mstate)
source("estimation_process.R")
#markovian case

set.seed(1234)

alpha <- 0.05
BSiter <- 1000 # maybe more ~2000 for final simus 
number_iteration <- 1000 # more for final simus ~5000

rate_0 <- 1/50
rate_1 <- 1/50
rate_2 <- 1/100
rate_cens <- 1/120

conf_level <- rep(0,1000)
# test_accepted <- rep(0,1000)
test_rejected <- rep(0,1000)

tra <- matrix(c(FALSE,TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE),3,3,byrow=TRUE)

# from illness to death
i <- 1
j <- 2
# time stamp
s <- 20
t <- 180

for (n in c(50,100,300,500,700,1000)){
  
  for (k in 1:number_iteration){
  
    # simulation
  
    start_state <- rep(0,n)
    T_01 <- rexp(n,rate_0) # time health to ill state
    T_02 <- rexp(n,rate_1) # time health to death state 
    T_12 <- rexp(n,rate_2) # time from ill to death state
    cens <- rexp(n,rate_cens) # time to censoring
    state_1 <- 1*(T_01 < pmin(T_02,cens)) + 2*(T_02 < pmin(T_01,cens)) + 3*(cens <= pmin(T_01,T_02))
    state_2 <- 2*(cens > T_01+T_12) + 3*(T_01 <= cens & cens <= T_01+T_12)
    
    # results in time-continuous markovian multi-state process'
    # state 3 = censored
    
    # next we transform the generated data into a data matrix that shows transition times and states per patient
    
    all_datas <- matrix(rep(0,5*2*n),nrow = 2*n)
    count <- 1
    for (m in 1:n){
      if(state_1[m] == 2){
        all_datas[count,] <- c(m,0,state_1[m],0,T_02[m])
      }
      else if(state_1[m] == 3){
        all_datas[count,] <- c(m,0,state_1[m],0,cens[m])
      }
      else if(state_1[m] == 1 & state_2[m] == 2){
        all_datas[count,] <- c(m,0,1,0,T_01[m])
        count <- count + 1
        all_datas[count,] <- c(m,1,2,T_01[m],T_01[m] + T_12[m])
      }
      else if(state_1[m] == 1 & state_2[m] == 3){
        all_datas[count,] <- c(m,0,1,0,T_01[m])
        count <- count + 1
        all_datas[count,] <- c(m,1,3,T_01[m],cens[m])
      }
      count <- count + 1
    }
    
    datas <- data.frame(all_datas[-which(all_datas[,1] == 0),])
    names(datas) <- c("id","from","to","entry","exit")
    
    conf_level[k] <- estimation(datas,s,t,i,j,BSiter,tra)
    test_rejected[k] <- (conf_level[k] <= alpha) 
    # BSstat_pred_sort <- sort(BSstat_pred)
    # Cn <- BSstat_pred_sort[BSiter*(1-alpha)]
    # test_accepted[k] <- (teststat_pred_markov <= Cn) #accept H0 with conf level alpha?
    }
  string <- paste("nonmarkov",number_iteration,"n",n,".RData",sep = "")
  save.image(string)
}