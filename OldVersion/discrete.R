rm(list = ls())

library(etm)
library(mstate)
source("estimation_process.R")
# discrete case: First non-markovian, then markovian
set.seed(1234)


alpha <- 0.05
BSiter <- 50
number_iteration <- 50 # more for final simus ~5000


lambda_01 <- 15
lambda_02 <- 20
lambda_12 <- 15
rate_cens <- 1/300

p <- 0.8  

conf_level <- rep(0,number_iteration)
test_rejected <- rep(0,number_iteration)

tra <- matrix(c(FALSE,TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE),3,3,byrow=TRUE)

# from illness to death
i <- 1
j <- 2
# time stamp
s <- 30
t <- 100

size <- 50
coinbox <- rexp(size,1/2)
# box groesse 50, replizieren (aneinanderhaengen) um Markov-Fall zu simulieren, oder einfach mit zuruecklegen ziehen (um Markov Fall zu generieren)
# fuer stetigen Fall wie in Titman Paper
# evtl. spaeter noch fuer andere Zensierungsraten

# Markovian case can be simulated, if we draw the transfer time from state_1 to state_2 from a different box
for (n in c(500,700,1000)){
  for (k in 1:number_iteration){
  
    nc01 <- pmin(rpois(n,lambda_01)+1, size/2) #number of coins to draw in state 0 to state 1
    nc02 <- pmin(rpois(n,lambda_02)+1, size/2) #number of coins to draw in state 0 to state 2
    nc12 <- pmin(rpois(n,lambda_12)+1, size/2)  #number of coins to draw in state 1 to state 2
    r <- runif(n,0,1)
    cens <- rexp(n,rate_cens)
    
    # 'draw' coins from the box with replacement and note the sum to get the time of transfer between states
    # results in a markov multi state process
    
    T_01 <- rep(0,n)
    T_02 <- rep(0,n)
    T_12 <- rep(0,n)
    
    for (m in 1:n){
      zero_ind <- sample(size,nc01[m])
      zero_samp <- coinbox[zero_ind]
      T_01[m] <- sum(zero_samp)
      remain <- coinbox[-zero_ind]
      zero_ind <- sample(size,nc02[m])
      zero_samp <- coinbox[zero_ind]
      T_02[m] <- sum(zero_samp)
      one_ind <- sample(size - nc01[m],nc12[m])
      one_samp <- remain[one_ind]
      T_12[m] <- sum(one_samp)
    }
    
    state_1 <- 1*(T_01 < pmin(T_02,cens)) + 2*(T_02 < pmin(T_01,cens)) + 3*(cens <= pmin(T_01,T_02))
    state_2 <- 2*(cens > T_01+T_12) + 3*(T_01 <= cens & cens <= T_01+T_12)
    
    #state 3 = censored
    
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
    
    data <- data.frame(all_data[-which(all_data[,1] == 0),])
    names(data) <- c("id","from","to","entry","exit")
    
    conf_level[k] <- estimation(data,s,t,i,j,BSiter,tra)
    test_rejected[k] <- (conf_level[k] <= alpha) 
   
  }
  string <- paste("nonmarkovdisc",number_iteration,"n",n,".RData",sep = "")
  save.image(string)
}
coinbox <- rep(coinbox,50)
# next: Markovian case
# Markovian case can be simulated, if we draw the transfer time from state_1 to state_2 from a different box
for (n in c(500,700,1000)){
  for (k in 1:number_iteration){
    nc01 <- rpois(n,lambda_01)+1 #number of coins to draw in state 0 to state 1
    nc02 <- rpois(n,lambda_02)+1 #number of coins to draw in state 0 to state 2
    nc12 <- rpois(n,lambda_12)+1 #number of coins to draw in state 1 to state 2
    r <- runif(n,0,1)
    cens <- rexp(n,rate_cens)
    
    # 'draw' coins from the box with replacement and note the sum to get the time of transfer between states
    # results in a markov multi state process
    
    T_01 <- rep(0,n)
    T_02 <- rep(0,n)
    T_12 <- rep(0,n)
    
    for (m in 1:n){
      zero_ind <- sample(size,nc01[m], replace = TRUE)
      zero_samp <- coinbox[zero_ind]
      T_01[m] <- sum(zero_samp)
      zero_ind <- sample(size,nc02[m], replace = TRUE)
      zero_samp <- coinbox[zero_ind]
      T_02[m] <- sum(zero_samp)
      one_ind <- sample(size,nc12[m], replace = TRUE)
      one_samp <- coinbox[one_ind]
      T_12[m] <- sum(one_samp)
    }

    state_1 <- 1*(T_01 < pmin(T_02,cens)) + 2*(T_02 < pmin(T_01,cens)) + 3*(cens <= pmin(T_01,T_02))
    state_2 <- 2*(cens > T_01+T_12) + 3*(T_01 <= cens & cens <= T_01+T_12)
    
    #state 3 = censored
    
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
    
    data <- data.frame(all_data[-which(all_data[,1] == 0),])
    names(data) <- c("id","from","to","entry","exit")
    
    conf_level[k] <- estimation(data,s,t,i,j,BSiter,tra)
    test_rejected[k] <- (conf_level[k] <= alpha) 
    
  }
  string <- paste("markovdisc",number_iteration,"n",n,".RData",sep = "")
  save.image(string)
}


