#Source the functions into the R invironment
setwd("~/Dropbox/research/Active/Dipankar/Analysis/Data/BP/Ht")

source('~/Dropbox/research/Active/Dipankar/Analysis/Data/BP/Ht/Datastep.R')
source('~/Dropbox/research/Active/Dipankar/Analysis/Data/BP/Ht/MCMC_BP_univariateRE.R')
nrun =400000 ; nskip = 10 ; nskip_r = 40; nburn=20000; nburn_r=5000;Jw=16;a_alpha = 1; b_alpha=1;maxc = 50;m=28;
d_r = 2; BP=1; a_alpha = 1; b_alpha = 1;SR=1
th_init_L = c(-4,1)
dataresult_AFT_LogL<-mcmc("AFT",3,maxc,t1,t2,type_d,m,BP,SR,crx1,FTx1,1,
                          nrun,nskip,nskip_r,nburn,nburn_r,Jw,a_alpha,b_alpha,th_init_L)
save.image("~/Dropbox/research/Active/Dipankar/Analysis/Data/BP/Ht/NPA_univariate.RData")

