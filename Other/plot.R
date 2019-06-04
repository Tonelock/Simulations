samplesizes <- c(50,100,300,500,700,1000,2000)
markovvalues <- c(0.126,0.128,0.102,0.074,0.086,0.068,0.056)
nonmarkovvalues <- c(0.204,0.256,0.472,0.628,0.806,0.892,0.996)

plot(samplesizes,markovvalues,type = 'l')
plot(samplesizes,markovvalues, type="l", ylim=c(0,1), xlab="samplesize", ylab="rejection rate")
lines(samplesizes,nonmarkovvalues, type="l", col=2)
lines(c(-10,2500),c(0.05,0.05), lty = 'dotted')
legend("topleft", legend=c("Markov Case", "Non-Markov case"), col=1:2, lty=1)