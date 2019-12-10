source('MCMC_BP_univariateRE.R')

#---------------------------------------------------------------------------------------------------#
# Model setup
m = 30; n = 10; model = "PH";distr=3;BP=1
crcoef_s = c(-1,0.2);pcr = 2;FTcoef_s = c(0.5,-0.5);pFT=2;phi_s = 0.3
#    n---number of individuals under study.
#    m---number of units within each individual.
#    d_r---dimension of the random effect vector u_i.
#    model---the survival model for the latent distribution.
#    distribution---the assumed parametric centering distribution for the latent distribution.
#    crcoef_s---cure rate coefficients with dimesion pcr.
#    FTcoef_s--- coefficients in the survival model of the latent distribution; dimension of it is pFT. 
#    phi_s---the coefficient linking latent distribution random effects and cure rate random effects.FTw = phi*cru. 

#----------------------------------------------------------------------------------------------------#
#set.seed(12348)

#------------------------------------------------------------------------------------------------------------#
#Simulate observed failure times, censoring indicators, and covariates 
crx = cbind(rep(1,n*m),rnorm(n*m,0,0.5))
FTx = cbind(rbinom(n*m,1,0.5),rnorm(n*m,0,0.5))

sigmau = 0.25

#    adj.g---adjacency matrix
#    adj.g.wish---the upper left half of adj.g
#    Dn---matrix with diagnoal elements being the number of neighbors for each unit
#    K---A sample of precision matrix from G-Wishart distribution with graph Dn-0.9*adj.g and 4 degrees of freedom
#    sigmau--covariance matrix of the random effects


cru_s = rep(0,n);FTw_s = rep(0,n)
for(i in 1:m){
  cru_s[i] = rnorm(1,0,sqrt(sigmau))
  FTw_s[i] = phi_s*cru_s[i]
}

#    cru---cure rate random effects, sampled from multivariate normal with mean zero and covariance sigmau.
#    FTw---latent distribution random effects
nt = sampleTs(model,n,m,crx, rep(cru_s,each=n),crcoef_s, FTx,FTcoef_s, rep(FTw_s,each=n))
type=rep(0,n*m);ci = rep(seq(1,m,1),each = n)

#t1--lower bound of the observed interval for failure time or observed; default is zero
#t2--upper bound of the observed interval for failure time or observed; default is Inf

ind = 1;t1=rep(0,n*m);t2=rep(Inf,n*m);censortu=NULL;censortl=NULL
for(i in 1:(n*m)){
  censortu[i] = 20
  censortl[i] = rexp(1,1)
   if(nt[i]<censortl[i]){type[i]=2;t2[i]=censortl[i]}
    if(nt[i]>censortu[i]){type[i]=3;t1[i]=censortu[i]}
    if(nt[i]<censortu[i] && nt[i]>censortl[i]){
      if(runif(1)<0.6){type[i] = 1;t1[i]=nt[i];t2[i] = nt[i]}
      else{type[i] = 4;t1[i]=censortl[i];t2[i]=censortu[i]}
      }

}

data = list(t1 = t1, t2 = t2,type = type, FTx = FTx, crx = crx, ci = ci)
mcmc.setup = list(nrun = 40000, nburn = 1000, nskip = 5)
BP.setup = list(Jw = 20, a.alpha = 1, b.alpha=1)
th.initial=c(-2,0.5)

mcmc.init (model="PH", distr=3, SR=1,data, BP.setup$Jw, th.initial)

simulresult<-mcmc(model = "PH",BP=1,SR=1,distr=3,data,mcmc.setup,BP.setup,th.initial)


par(mfrow=c(3,3))
plot(simulresult$theta[,1],type="l",ylab=expression(theta[1]),xlab="iteration")
plot(simulresult$theta[,2],type="l",ylab=expression(theta[2]),xlab="iteration")
plot(simulresult$crcoef[,1],type="l",ylab=expression(gamma[1]),xlab="iteration")
plot(simulresult$crcoef[,2],type="l",ylab=expression(gamma[2]),xlab="iteration")
plot(simulresult$FTcoef[,1],type="l",ylab=expression(beta[1]),xlab="iteration")
plot(simulresult$FTcoef[,2],type="l",ylab=expression(beta[2]),xlab="iteration")
plot(simulresult$weight[,1],type="l",ylab=expression(w[1]),xlab="iteration")
plot(simulresult$weight[,2],type="l",ylab=expression(w[1]),xlab="iteration")
plot(simulresult$alpha,type="l",ylab=expression(alpha),xlab="iteration")



cbind(crcoef_s,apply(simulresult$crcoef,2,mean))
cbind(FTcoef_s,apply(simulresult$FTcoef,2,mean))

round(mean(simulresult$sigma),2)
round(sigmau,2)
mean(simulresult$phi)

par(mfrow=c(3,1))
plot(simulresult$sigma,type="l",ylab=expression(Sigma),xlab="iteration")
plot(simulresult$random_effect[,3],type="l",ylab=expression(u[1]),xlab="iteration")
plot(simulresult$phi,type="l",ylab=expression(phi),xlab="iteration")

t=seq(0.1,30,0.1)
densit0 = NULL
survt0=NULL
densit = rep(0,length(t))
survt = rep(0,length(t))
seq = 1:500

for(k in 1:length(seq)){
  for(i in 1:length(t)){
    th1 = simulresult$theta[seq[k],1]
    th2 = simulresult$theta[seq[k],2]
    w = simulresult$weight[seq[k],]
    densit[i] = densit[i]+exp(logf0BP(t[i], th1, th2, w, 1, 1))/length(seq)
    survt[i] = survt[i]+S0BP(t[i], th1, th2, w, 1, 1)/length(seq)
    
  }
}


for(i in 1:length(t)){
  th1 = mean(simulresult$theta[seq,1])
  th2 = mean(simulresult$theta[seq,2])
  w = apply(simulresult$weight[seq,],2,mean)
  densit0[i] = densi0(t[i])
  survt0[i] = 1-F0(t[i])
 
}

par(mfrow=c(2,1))
plot(t,densit,type="l",col="red",ylim=c(0,0.2))
lines(t,densit0,col="blue")

plot(t,survt,type="l",col="red",ylim=c(0,1))
lines(t,survt0,col="blue")
