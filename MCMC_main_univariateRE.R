source('MCMC_BP_univariateRE.R')

#---------------------------------------------------------------------------------------------------#
# Model setup
m = 30; n = 30; model = "PH";distr=2;BP=1;SR=3
crcoef_s = c(1,0.2);pcr = 2;FTcoef_s = c(0.5,-0.5);pFT=2
nsimul = 1

t=seq(0.1,30,0.1)
densit0 = NULL
survt0=NULL
densit = matrix(0,ncol = length(t),nrow= nsimul)
survt = matrix(0,ncol = length(t),nrow= nsimul)
densit.p = matrix(0,ncol = length(t),nrow= nsimul)
survt.p = matrix(0,ncol = length(t),nrow= nsimul)
crcoefresult = matrix(0,ncol=4,nrow=nsimul)
FTcoefresult = matrix(0,ncol=4,nrow=nsimul)
if(SR==1|SR==2|SR==3){sigmaresult = matrix(0,ncol=2,nrow=nsimul) }
if(SR==4){sigmaresult = matrix(0,ncol=4,nrow=nsimul) }
if(SR==1){phiresult = matrix(0,ncol=2,nrow=nsimul)}
#    n---number of individuals under study.
#    m---number of units within each individual.
#    d_r---dimension of the random effect vector u_i.
#    model---the survival model for the latent distribution.
#    distr--1:Logistic distribution; 2:Normal distribution; 3: Weibull distribution
#    crcoef_s---cure rate coefficients with dimesion pcr.
#    FTcoef_s--- coefficients in the survival model of the latent distribution; dimension of it is pFT. 
#    phi_s---the coefficient linking latent distribution random effects and cure rate random effects.FTv = phi*crv. 

#----------------------------------------------------------------------------------------------------#
#set.seed(12348)
for(ns in 1:nsimul){
#------------------------------------------------------------------------------------------------------------#
#Simulate observed failure times, censoring indicators, and covariates 
crx = cbind(rep(1,n*m),rnorm(n*m,0,1))
FTx = cbind(rbinom(n*m,1,0.5),rnorm(n*m,0,1))


#    adj.g---adjacency matrix
#    adj.g.wish---the upper left half of adj.g
#    Dn---matrix with diagnoal elements being the number of neighbors for each unit
#    K---A sample of precision matrix from G-Wishart distribution with graph Dn-0.9*adj.g and 4 degrees of freedom
#    sigmau--covariance matrix of the random effects

if(SR==0){crv_s = rep(0,m);FTv_s = rep(0,m)}
if(SR==1){ sigmau = 0.25;phi_s = 0.3;crv_s = rnorm(m,0,sqrt(sigmau));FTv_s = phi_s*crv_s}
if(SR==2){sigmau = 0.25; crv_s = rnorm(m,0,sqrt(sigmau));FTv_s = rep(0,m)}
if(SR==3){sigmau = 0.25;FTv_s = rnorm(m,0,sqrt(sigmau));crv_s = rep(0,m)}
if(SR==4){sigmau = c(0.25,0.25); crv_s = rnorm(m,0,sqrt(sigmau));FTv_s = rnorm(m,0,sqrt(sigmau))}

#    crv---cure rate random effects, sampled from multivariate normal with mean zero and covariance sigmau.
#    FTv---latent distribution random effects
nt = sampleTs(model,n,m,crx, rep(crv_s,each=n),crcoef_s, FTx,FTcoef_s, rep(FTv_s,each=n))
ci = rep(seq(1,m,1),each = n)

#t1--lower bound of the observed interval for failure time or observed; default is zero
#t2--upper bound of the observed interval for failure time or observed; default is Inf

type=rep(0,n*m);ind = 1;t1=rep(0,n*m);t2=rep(Inf,n*m);censortu=NULL;censortl=NULL
for(i in 1:(n*m)){
  censortu[i] = 30
  censortl[i] = rexp(1,1)
   if(nt[i]<censortl[i]){type[i]=2;t2[i]=censortl[i]}
    if(nt[i]>censortu[i]){type[i]=3;t1[i]=censortu[i]}
    if(nt[i]<censortu[i] && nt[i]>censortl[i]){
      if(runif(1)<0.9){type[i] = 1;t1[i]=nt[i];t2[i] = nt[i]}
      else{type[i] = 4;t1[i]=censortl[i];t2[i]=censortu[i]}
      }

}

i = 5;type[i];t1[i];t2[i]
data = list(t1 = t1, t2 = t2,type = type, FTx = FTx, crx = crx, ci = ci)
mcmc.setup = list(nrun = 20000, nburn = 1000, nskip = 5)
BP.setup = list(Jw = 15, a.alpha = 0.1, b.alpha=0.1)
th.initial=c(-2,0.5)

mcmc.int.o <- mcmc.init (model, distr, SR,data, BP.setup$Jw, th.initial)

simulresult<-mcmc(model,BP,SR,distr,data,mcmc.setup,BP.setup,th.initial)



crcoefresult[ns,] = cbind(crcoef_s,apply(simulresult$crcoef,2,mean))
FTcoefresult[ns,] = cbind(FTcoef_s,apply(simulresult$FTcoef,2,mean))
if(SR==1|SR==2|SR==3){sigmaresult[ns,] = c(sigmau,mean(simulresult$sigma))}
if(SR==4){sigmaresult[ns,] = c(sigmau,mean(simulresult$sigma))}
if(SR==1){phiresult[ns,] = c(phi_s,mean(simulresult$phi))}

seq = 1:length(simulresult$theta[,1])
for(k in 1:length(seq)){
  for(i in 1:length(t)){
    th1 = simulresult$theta[seq[k],1]
    th2 = simulresult$theta[seq[k],2]
    w = simulresult$weight[seq[k],]
    densit[ns,i] = densit[ns,i]+exp(logf0BP(t[i], th1, th2, w, BP, distr))/length(seq)
    survt[ns,i] = survt[ns,i]+S0BP(t[i], th1, th2, w, BP, distr)/length(seq)
  }
}

 for(i in 1:length(t)){
   w = apply(simulresult$weight[seq,],2,mean)
   densit.p[ns,i] = exp(logf0BP(t[i], mcmc.int.o$theta[1], mcmc.int.o$theta[2], w, 0, distr))
   survt.p[ns,i] = S0BP(t[i], mcmc.int.o$theta[1], mcmc.int.o$theta[2], w, 0, distr)
  }

}

for(i in 1:length(t)){
  densit0[i] = densi0(t[i])
  survt0[i] = 1-F0(t[i])
}


par(mfrow=c(3,3))
plot(simulresult$theta[,1],type="l",ylab=expression(theta[1]),xlab="iteration")
plot(simulresult$theta[,2],type="l",ylab=expression(theta[2]),xlab="iteration")
plot(simulresult$crcoef[,1],type="l",ylab=expression(gamma[1]),xlab="iteration")
plot(simulresult$crcoef[,2],type="l",ylab=expression(gamma[2]),xlab="iteration")
plot(simulresult$FTcoef[,1],type="l",ylab=expression(beta[1]),xlab="iteration")
plot(simulresult$FTcoef[,2],type="l",ylab=expression(beta[2]),xlab="iteration")
plot(simulresult$weight[,3],type="l",ylab=expression(w[1]),xlab="iteration")
plot(simulresult$weight[,4],type="l",ylab=expression(w[1]),xlab="iteration")
plot(simulresult$alpha,type="l",ylab=expression(alpha),xlab="iteration")



par(mfrow=c(2,2))
plot(simulresult$sigma,type="l",ylab=expression(Sigma),xlab="iteration")
plot(simulresult$crv[,3],type="l",ylab=expression(u[1]),xlab="iteration")
plot(simulresult$phi,type="l",ylab=expression(phi),xlab="iteration")

plot(t,apply(densit,2,mean),type="l",col="red",ylim=c(0,0.2))
lines(t,densit0,col="blue")
lines(t,apply(densit.p,2,mean),col="black")
legend("topright",legend=c("NP","T","P"),col=c("red","blue","black"),lty=1)

plot(t,apply(survt,2,mean),type="l",col="red",ylim=c(0,1))
lines(t,survt0,col="blue")
lines(t,apply(survt.p,2,mean),col="black")
legend("topright",legend=c("NP","T","P"),col=c("red","blue","black"),lty=1)

apply(crcoefresult,2,mean)
apply(FTcoefresult,2,mean)
