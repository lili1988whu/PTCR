setwd("~/Dropbox/research/Active/Dipankar/Analysis/Data/BP/Ht")
source('MCMC_BP_multivariateRE.R')

#---------------------------------------------------------------------------------------------------#
# Model setup
n = 10;m = 20;d_r = 2; model = "AFT";distr=1;BP=1
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

Randomx = matrix(0,ncol=d_r,nrow=n*m)
Randomx[,1] = rnorm(n*m,0,1)
Randomx[,2] = rnorm(n*m,0,1)


#    crx---covariates in the model for cure rate
#    FTx---covariates in the survival model for the latent distribution

sigmau = cbind(c(0.5,0),c(0,0.5))


cru_s = rep(0,n*d_r);FTw_s = rep(0,n*d_r)
for(i in 1:n){
  seqi = ((i-1)*d_r+1):(i*d_r)
  if(d_r>1){  cru_s[seqi] = rmnorm(rep(0,d_r), sigmau)}
  else{  cru_s[seqi] = rnorm(1,0,sigmau)}
  FTw_s[seqi] = phi_s*cru_s[seqi]
}
#    cru---cure rate random effects, sampled from multivariate normal with mean zero and covariance sigmau.
#    FTw---latent distribution random effects
nt = sampleTs(model,n,m,d_r,crx, cru_s,crcoef_s, FTx,FTcoef_s, FTw_s,Randomx)
type=rep(0,n*m)

ind = 1;t1=rep(0,n*m);t2=rep(0,n*m);censortu=NULL;censortl=NULL
for(i in 1:(n*m)){
  censortu[i] = 20
  censortl[i] = rexp(1,1)
    if(nt[i]<censortl[i]){type[i]=2;t1[i] = censortl[i]; t2[i]=censortl[i]}
    if(nt[i]>censortu[i]){type[i]=3;t1[i]=censortu[i]; t2[i] = censortu[i]}
    if(nt[i]<censortu[i] && nt[i]>censortl[i]){
      if(runif(1)<0.6){type[i] = 1;t1[i]=nt[i]}
      else{type[i] = 4;t1[i]=censortl[i];t2[i]=censortu[i]}
      }

}

#Setup for MCMC
nrun =50000 ; nskip = 5;nskip_r=50; nburn=1000;nburn_r = 100;Jw=20;a_alpha = 1; b_alpha=1;maxc = 15;BP=1;distr=1;SR=1
#    FTJ--number of level for the tail-free process
#    nrun--- total number of iterations

th_initial=c(-2,0.5)

simulresult<-mcmc(model,distr,maxc,t1,t2,type,m,d_r,BP,SR,crx,FTx,1,Randomx,nrun,nskip,nskip_r,nburn,nburn_r,Jw,a_alpha,b_alpha,th_initial)

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

sigmaumean=matrix(0,ncol=d_r,nrow=d_r)
for(i in 1:d_r){
  for(j in 1:d_r){
    sigmaumean[i,j] = mean(simulresult$sigma[,i,j])
  }
}

cbind(crcoef_s,apply(simulresult$crcoef,2,mean))
cbind(FTcoef_s,apply(simulresult$FTcoef,2,mean))

round(sigmaumean,2)
round(sigmau,2)
mean(simulresult$phi)
mean(simulresult$GPphi)
plot(simulresult$sigma[,1,1],type="l",ylab=expression(Sigma[11]),xlab="iteration")
plot(simulresult$random_effect[,1,4],type="l",ylab=expression(u[11]),xlab="iteration")
plot(simulresult$phi,type="l",ylab=expression(phi),xlab="iteration")
plot(simulresult$GPtau,type="l",ylab=expression(tau),xlab="iteration")
plot(simulresult$GPphi,type="l",ylab=expression(tau),xlab="iteration")

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
