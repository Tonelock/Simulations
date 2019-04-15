rm(list = ls())

library(etm)
library(mstate)

# simulation

n <- 3000
id <- 1:n
size <- 500
box <- rexp(size,1/2)
sick <- rpois(n,40) #number of coins to draw until sick
death <- rpois(n,50)  #number of coins to draw until dead
earldeath <- rpois(n,50)  #number of coins to draw until early death (transfer from healthy to dead state)
censored <- rpois(n,100) #number of coins to draw until censored
T_01 <- rep(0,n)
T_02 <- rep(0,n)
T_12 <- rep(0,n)
cens <- rep(0,n)

# 'draw' coins from the box without replacement and note the sum to get the time of transfer between states
# results in a time-continuous non-markov multi state process
# again, the first transition that happens to an individual is considered to be the one that actually happens

for (i in 1:n){

  sick_ind <- sample(size,sick[i])
  sick_samp <- box[sick_ind]
  T_01[i] <- sum(sick_samp)
  remain <- box[-sick_samp]
  
  death_ind <- sample(size - sick[i],death[i])
  death_samp <- remain[death_ind]
  T_12[i] <- T_01[i] + sum(death_samp)
  
  earldeath_ind <- sample(size,earldeath[i])
  earldeath_samp <- box[earldeath_ind]
  T_02[i] <- sum(earldeath_samp)
  
  cens_ind <- sample(size,censored[i])
  cens_samp <- box[cens_ind]
  cens[i] <- sum(cens_samp)
}

start <- 10

#state 3 = censored

entry_state <-  (cens > start) * {0*(start < pmin(T_01,T_02)) + 1*((T_01 < pmin(T_02,start)) & (start < T_12)) + 2*((T_02 < start) | (T_12 < start))} + 
  3*(cens < start)
first_transition <- pmin(cens,T_01,T_02)
all_data <- matrix(rep(0,5*4*n),nrow = 4*n)
count <- 1
for (i in id){
  if (entry_state[i] == 0){
    
    if (T_01[i] < cens[i]){
      all_data[count,] <- c(i,0,1,0,T_01[i])    
      count <- count + 1
      
      if(T_12[i] < cens[i]){
        all_data[count,] <- c(i,1,2,T_01[i],T_12[i])
      }
      else{
        all_data[count,] <- c(i,1,3,T_01[i],cens[i])
      }
      
    }
    else{
      all_data[count,] <- c(i,0,3,0,cens[i])
    }
  }
  else if(entry_state[i] == 1){
    
    if(T_12[i] < cens[i]){
      all_data[count,] <- c(i,1,2,0,T_12[i])
    }
    else{
      all_data[count + 1,] <- c(i,1,3,0,cens[i])
    }
    
  }
  
  count <- count + 1
}

data <- data.frame(all_data[-which(all_data[,1] == 0),])
names(data) <- c("id","from","to","entry","exit")

tra <- matrix(c(F,T,T,F,F,T,F,F,F),3,3,byrow=T)

# from illness to death
i <- 1
j <- 2
# time stamp
s <- 70
t <- 250

# estimate transition probabilities
etmAJE_nonmarkov <- etm(data, state.names = 0:2, tra=tra, cens.name=3, s=s, t=t, covariance=FALSE, delta.na=FALSE)
PijAJE_nonmarkov <- c(0,etmAJE_nonmarkov$est[i + 1,j + 1,])
timesAJE_nonmarkov <- c(s,etmAJE_nonmarkov$time)

good_ind_nonmarkov <- which(data$entry <= s & data$exit >s & data$from == i)
rel_ids_nonmarkov <- unique(data$id[good_ind_nonmarkov])
npeople_nonmarkov <- length(rel_ids_nonmarkov)
# "good" data
data_s <- data[data$id %in% rel_ids_nonmarkov,]
# Remove all intervals ending before s
data_s <- data_s[data_s$exit >= s,]
# Truncate all intervals starting before s
data_s$entry[data_s$entry < s] <- s

etm_s_nonmarkov <- etm(data_s, state.names = 0:2, tra=tra, cens.name=3, s=s, t=t, covariance=FALSE, delta.na=FALSE)
F0_nonmarkov <- c(0,etm_s_nonmarkov$est[i + 1,j + 1,])
timesTE_nonmarkov <- c(s,etm_s_nonmarkov$time)

h <- 0.5
x_nonmarkov <- seq(s,t,h)
TE_nonmarkov <- approxfun(timesTE_nonmarkov, F0_nonmarkov, method="constant")
plot(x_nonmarkov, ifelse(is.na(TE_nonmarkov(x_nonmarkov)), max(TE_nonmarkov(x_nonmarkov), na.rm=T), TE_nonmarkov(x_nonmarkov)), type="l", ylim=c(0,1), xlab="time (days)", ylab="prob critical to dead")
AJE_nonmarkov <- approxfun(timesAJE_nonmarkov, PijAJE_nonmarkov, method="constant")
lines(x_nonmarkov, ifelse(is.na(AJE_nonmarkov(x_nonmarkov)), max(AJE_nonmarkov(x_nonmarkov), na.rm=T), AJE_nonmarkov(x_nonmarkov)), type="l", col=2)
legend("bottomright", legend=c("Titman estimator", "Aalen-Johansen estimator"), col=1:2, lty=1)

# the plot shows that Titman estimator and Aalen-Johansen estimator converge to different distributions for the transition probability

TE_pred_nonmarkov <- ifelse(is.na(TE_nonmarkov(x_nonmarkov)), max(TE_nonmarkov(x_nonmarkov), na.rm=T), TE_nonmarkov(x_nonmarkov))
AJE_pred_nonmarkov <- ifelse(is.na(AJE_nonmarkov(x_nonmarkov)), max(AJE_nonmarkov(x_nonmarkov), na.rm=T), AJE_nonmarkov(x_nonmarkov))

teststat_pred_nonmarkov <- n*sum((AJE_pred_nonmarkov-TE_pred_nonmarkov)^2 * h)

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


 max_id_nonmarkov <- max(data$id)
 alpha <- 0.95
 BSiter <- 1000
 set.seed(1234)
 BSstat_pred <- replicate(BSiter, bs_one_iter(n, data, max_id_nonmarkov, AJE_pred_nonmarkov, TE_pred_nonmarkov, times=x_nonmarkov, h, tra))
 sum(BSstat_pred >= teststat_pred_nonmarkov)/BSiter

 BSstat_pred_sort <- sort(BSstat_pred)
 Cn <- BSstat_pred_sort[BSiter*(1-alpha)]
 (teststat_pred_nonmarkov <= Cn)