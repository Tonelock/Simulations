rm(list = ls())

library(etm)
library(mstate)
source("estimation_process.R")
# simulation
set.seed(1234)


alpha <- 0.05
BSiter <- 500
number_iteration <- 50 # more for final simus ~5000


lambda_0 <- 10    
lambda_1 <- 15
rate_cens <- 1/120

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
box <- rexp(size,1/2)
# box groesse 50, replizieren (aneinanderhaengen) um Markov-Fall zu simulieren, oder einfach mit zurücklegen ziehen (um Markov Fall zu generieren)
# fuer stetigen Fall wie in Titman Paper
# evtl. spaeter noch fuer andere Zensierungsraten

# Markovian case can be simulated, if we draw the transfer time from state_1 to state_2 from a different box
for (n in c(1000)){
  for (k in 1:number_iteration){
  
    nc0 <- rpois(n,lambda_0)+1 #number of coins to draw in state 0
    nc1 <- rpois(n,lambda_1)+1  #number of coins to draw in state 1  
    r <- runif(n,0,1)
    cens <- rexp(n,rate_cens)
    
    # 'draw' coins from the box without replacement and note the sum to get the time of transfer between states
    # results in a non-markov multi state process
    
    T_0 <- rep(0,n)
    T_1 <- rep(0,n)
    
    for (m in 1:n){
      zero_ind <- sample(size,min(nc0[m],size - 2))
      zero_samp <- box[zero_ind]
      T_0[m] <- sum(zero_samp)
      remain <- box[-zero_samp]
      if(nc0[m] > (size-2)){one_ind <- 1} else {one_ind <- sample(size - nc0[m],nc1[m])}
      one_samp <- remain[one_ind]
      T_1[m] <- sum(one_samp)
    }
    state_1 <- (1*(r < p) + 2*(p <= r))*(cens > T_0) + 3*(cens <= T_0)
    state_2 <- 2*(cens > T_0+T_1) + 3*(T_0 <= cens & cens <= T_0+T_1) + 0*(cens <= T_0)
    
    #state 3 = censored
    
    all_data <- matrix(rep(0,5*2*n),nrow = 2*n)
    count <- 1
    for (m in 1:n){
      if(state_1[m] == 2){
        all_data[count,] <- c(m,0,state_1[m],0,T_0[m])
      }
      else if(state_1[m] == 3){
        all_data[count,] <- c(m,0,state_1[m],0,cens[m])
      }
      else if(state_1[m] == 1 & state_2[m] == 2){
        all_data[count,] <- c(m,0,1,0,T_0[m])
        count <- count + 1
        all_data[count,] <- c(m,1,2,T_0[m],T_0[m] + T_1[m])
      }
      else if(state_1[m] == 1 & state_2[m] == 3){
        all_data[count,] <- c(m,0,1,0,T_0[m])
        count <- count + 1
        all_data[count,] <- c(m,1,3,T_0[m],cens[m])
      }
      count <- count + 1
    }
    
    data <- data.frame(all_data[-which(all_data[,1] == 0),])
    names(data) <- c("id","from","to","entry","exit")
    
    conf_level[k] <- estimation(data,s,t,i,j,BSiter,tra)
    test_rejected[k] <- (conf_level[k] <= alpha) 
   
  }
  string <- paste("nonmarkov",number_iteration,"n",n,".RData",sep = "")
  save.image(string)
}

