library(etm)
library(mstate)
data("prothr")
head(prothr)
prothr[150:200,]


data <- prothr
data$from <- as.integer(data$from)
data$to <- as.integer(data$to)
# if Tstop before Tstart, increase Tstop by 0.5
data$Tstop[which(data$Tstart >= data$Tstop)] <- data$Tstop[which(data$Tstart >= data$Tstop)] + 0.5
n <- length(unique(data$id))

# last transitions per id
last_trans <- cumsum(unname(table(data$id)))
# censoring 
data$to[last_trans[which(data$status[last_trans] == data$status[last_trans-1])]] <- "cens"
data$status[last_trans[which(data$status[last_trans] == data$status[last_trans-1])]] <- "cens"
data <- subset(data, status!=0)
# first transitions per id
first_trans <- cumsum(c(1,unname(table(data$id))))
first_trans <- first_trans[-length(first_trans)]
first_trans <- first_trans[which(!(data$from[first_trans]== 2 & data$to[first_trans]== 1))]
last_trans <- cumsum(unname(table(data$id)))
# sort by transition times
data <- data[unique(sort(c(first_trans,last_trans))),]
data <- data[,c(1,2,3,5,6,8)]
names(data)[4:5] <- c("entry", "exit")
first_trans <- cumsum(c(1,unname(table(data$id))))
first_trans <- first_trans[-length(first_trans)]
data$entry[first_trans] <- 0

data_pred <- subset(data, treat=="Prednisone")
head(data_pred)

data_plac <- subset(data, treat=="Placebo")
head(data_plac)


tra <- matrix(c(F,T,T,F,F,T,F,F,F),3,3,byrow=T)
# from illness to death
i <- 2
j <- 3
# time stamp
s <- 1000


# first for Prednisone dataset
data <- data_pred

# estimate transition probabilities
etmAJE <- etm(data, state.names = 1:3, tra=tra, cens.name="cens", s=1000, t=3000, covariance=FALSE, delta.na=FALSE)
PijAJE <- c(0,etmAJE$est[i,j,])
timesAJE <- c(s,etmAJE$time)

first_trans <- cumsum(c(1,unname(table(data$id))))
first_trans <- first_trans[-length(first_trans)]

good_ind <- which(data$entry <= s & data$exit >s & data$from==i)
rel_ids <- unique(data$id[good_ind])
npeople <- length(rel_ids)
# "good" data
data_s <- data[data$id %in% rel_ids,]
#Remove all intervals ending before s
data_s <- data_s[data_s$exit >= s,]
#Truncate all intervals starting before s
data_s$entry[data_s$entry <s]<-s

etm_s <- etm(data_s, state.names = 1:3, tra=tra, cens.name="cens", s=1000, t=3000, covariance=FALSE, delta.na=FALSE)
F0 <- c(0,etm_s$est[i,j,])
timesTE <- c(s,etm_s$time)

#cbind(etm_s$time, etm_s$n.risk)

h <- 0.5
x <- seq(1000,3000,h)
TE <- approxfun(timesTE, F0, method="constant")
plot(x, ifelse(is.na(TE(x)), max(TE(x), na.rm=T), TE(x)), type="l", ylim=c(0,1), xlab="time (days)", ylab="prob critical to dead", main="prednisone group")
AJE <- approxfun(timesAJE, PijAJE, method="constant")
lines(x, ifelse(is.na(AJE(x)), max(AJE(x), na.rm=T), AJE(x)), type="l", col=2)
legend("bottomright", legend=c("Titman estimator", "Aalen-Johansen estimator"), col=1:2, lty=1)

TE_pred <- ifelse(is.na(TE(x)), max(TE(x), na.rm=T), TE(x))
AJE_pred <- ifelse(is.na(AJE(x)), max(AJE(x), na.rm=T), AJE(x))

teststat_pred <- n*sum((AJE_pred-TE_pred)^2 * h)



# next for placebo dataset
data <- data_plac

etmAJE <- etm(data, state.names = 1:3, tra=tra, cens.name="cens", s=1000, t=3000, covariance=FALSE, delta.na=FALSE)
PijAJE <- c(0,etmAJE$est[i,j,])
timesAJE <- c(s,etmAJE$time)

first_trans <- cumsum(c(1,unname(table(data$id))))
first_trans <- first_trans[-length(first_trans)]

good_ind <- which(data$entry <= s & data$exit >s & data$from==i)
rel_ids <- unique(data$id[good_ind])
npeople <- length(rel_ids)
data_s <- data[data$id %in% rel_ids,]
#Remove all intervals ending before s
data_s <- data_s[data_s$exit >= s,]
#Truncate all intervals starting before s
data_s$entry[data_s$entry <s]<-s

etm_s <- etm(data_s, state.names = 1:3, tra=tra, cens.name="cens", s=1000, t=3000, covariance=FALSE, delta.na=FALSE)
F0 <- c(0,etm_s$est[i,j,])
timesTE <- c(s,etm_s$time)

#cbind(etm_s$time, etm_s$n.risk)

h <- 0.5
x <- seq(1000,3000,h)
TE <- approxfun(timesTE, F0, method="constant")
plot(x, ifelse(is.na(TE(x)), max(TE(x), na.rm=T), TE(x)), type="l", ylim=c(0,1), xlab="time (days)", ylab="prob critical to dead", main="placebo group")
AJE <- approxfun(timesAJE, PijAJE, method="constant")
lines(x, ifelse(is.na(AJE(x)), max(AJE(x), na.rm=T), AJE(x)), type="l", col=2)
legend("bottomright", legend=c("Titman estimator", "Aalen-Johansen estimator"), col=1:2, lty=1)

TE_plac <- ifelse(is.na(TE(x)), max(TE(x), na.rm=T), TE(x))
AJE_plac <- ifelse(is.na(AJE(x)), max(AJE(x), na.rm=T), AJE(x))

teststat_plac <- n*sum((AJE_plac-TE_plac)^2 * h)


bs_one_iter <- function(n, data, max_id, AJE, TE, times=x, h, tra){
  
  multinom <- sample(1:n, replace=T)
  
  data_bs <- subset(data, id==data$id[multinom[1]])
  for(k in 2:n){
    data_temp <- subset(data, id==data$id[multinom[k]])
    data_temp$id <- data_temp$id + max_id * k
    data_bs <- rbind(data_bs, data_temp)
  }
  
  etmAJE <- etm(data_bs, state.names = 1:3, tra=tra, cens.name="cens", s=1000, t=3000, covariance=FALSE, delta.na=FALSE)
  PijAJE <- c(0,etmAJE$est[i,j,])
  timesAJE <- c(s,etmAJE$time)
  
  good_ind <- which(data$entry <= s & data$exit >s & data$from==i)
  rel_ids <- unique(data$id[good_ind])
  npeople <- length(rel_ids)
  data_s <- data[data$id %in% rel_ids,]
  #Remove all intervals ending before s
  data_s <- data_s[data_s$exit >= s,]
  #Truncate all intervals starting before s
  data_s$entry[data_s$entry <s]<-s
    
  etm_s <- etm(data_s, state.names = 1:3, tra=tra, cens.name="cens", s=1000, t=3000, covariance=FALSE, delta.na=FALSE)
  F0 <- c(0,etm_s$est[i,j,])
  timesTE <- c(s,etm_s$time)
  
  TE_bs <- approxfun(timesTE, F0, method="constant")
  AJE_bs <- approxfun(timesAJE, PijAJE, method="constant")
  
  TE_bs <- ifelse(is.na(TE_bs(times)), max(TE_bs(times), na.rm=T), TE_bs(times))
  AJE_bs <- ifelse(is.na(AJE_bs(times)), max(AJE_bs(times), na.rm=T), AJE_bs(times))
  return(n*sum((AJE_bs-AJE-TE_bs+TE)^2 * h))
}


max_id <- max(data$id)
alpha <- 0.05
BSiter <- 1000
set.seed(1234)
BSstat_pred <- replicate(BSiter, bs_one_iter(n, data_pred, max_id, AJE_pred, TE_pred, times=x, h, tra))
sum(BSstat_pred >= teststat_pred)/BSiter

BSstat_pred_sort <- sort(BSstat_pred)
Cn <- BSstat_pred_sort[BSiter*(1-alpha)]
(teststat_pred <= Cn)
#0.989
BSstat_plac <- replicate(BSiter, bs_one_iter(n, data_plac, max_id, AJE_plac, TE_plac, times=x, h, tra))
sum(BSstat_plac >= teststat_plac)/BSiter

BSstat_plac_sort <- sort(BSstat_plac)
Cn <- BSstat_plac_sort[BSiter*(1-alpha)]
(teststat_plac <= Cn)
#0.98