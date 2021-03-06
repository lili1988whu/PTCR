Rcpp::sourceCpp('MCMC_BP_univariateRE.cpp')
produ<-function(x,beta){
  if(length(beta)>1){return(x%*%beta)}
  else{return(x*beta)}
}

likelihood.wrapper<-function(model,BP,distr,data,parameters){
  #model---"AFT" or "PH" or "PO"
  #BP--0: no Bernstein Polynomial modeling of the baseline hazard function; 1: Otherwise
  #distr--1:Logistic distribution; 2:Normal distribution; 3: Weibull distribution
  #data includes:
               #t1--lower bound of the observed interval for failure time or observed; default is zero
               #t2--upper bound of the observed interval for failure time or observed; default is Inf
               #ci -- cluster indicator
               #type--1: failure time is observed; 2: failure time is left censored; 3: failure time is
                      #right censored; 4: failure time is interval-censored
               #FTx--covariates for the latent failure time distribution; could be NULL
               #crx--covariates for the cure rate 
  #parameter includes  
               #z--transformed weight parameters
               #th1, th2--parameters of the centering distribution
               #FTcoef--coefficients for FTx; could be NULL
               #crcoef--coefficients for crx
               #crv--random effects for cure rate; could be NULL
               #FTv--random effects for latent failure time distribution; could be NULL
  #browser()
  t1 = data$t1
  t2 = data$t2
  type = data$type
  FTx = data$FTx
  crx = data$crx
  ci = data$ci
  n.per.cluster = as.vector(table(ci))
  
  bz = parameters$bz
  th1 = parameters$theta[1]; th2 = parameters$theta[2]
  FTcoef = parameters$FTcoef; crcoef = parameters$crcoef
  
  Fvc = parameters$Fv[ci]
  crvc = parameters$crv[ci]
  
  dimFTcoef = if(is.null(dim(FTx))) NULL else dim(FTx)[2]
  FXbeta = if(is.null(dimFTcoef))  rep(0,length(ci)) else produ(FTx,FTcoef)
  dimcrcoef = if(is.null(dim(crx))) NULL else dim(crx)[2]
  crXbeta = if (is.null(dimcrcoef)) rep(0,length(ci)) else produ(crx,crcoef)
  weight = Ys_to_weight(bz)
  
  
  likelihoodall = likelihoodv(model,BP, t1,  t2, type, th1, th2, weight, distr, FXbeta, Fvc, crXbeta, crvc)
  tempdata = data.frame(likelihoodall,ci)
  likelihoodi = as.vector(tapply(tempdata$likelihoodall,INDEX = tempdata$ci, FUN = sum))
  likelihoodsum = sum( likelihoodall)  
  return(list(likelihoodi = likelihoodi,likelihoodsum = likelihoodsum))
}

mcmc.init <- function (model = "PH", distr=3, SR=1, data,Jw, th.initial){
  #model---"AFT" or "PH" or "PO"
  #distr--1:Logistic distribution; 2: log-Normal distribution; 3: Weibull distribution
  #SR random effects indicators:
     #0--no random effects
     #1--Fv = rho*crv
     #2--crv, no Fv
     #3--Fv, no crv
     #4--both crv and Fv
  #data includes:
       #t1--lower bound of the observed interval for failure time or observed; default is zero
       #t2--upper bound of the observed interval for failure time or observed; default is Inf
       #ci -- cluster indicator
       #type--1: failure time is observed; 2: failure time is left censored; 3: failure time is
           #right censored; 4: failure time is interval-censored
       #FTx--covariates for the latent failure time distribution; could be NULL
       #crx--covariates for the cure rate 
  
  t1 = data$t1
  t2 = data$t2
  type = data$type
  FTx = data$FTx
  crx = data$crx
  ci = data$ci
  m = length(unique(ci))
  pcr = if(is.null(crx)) NULL else dim(crx)[2]
  pFT = if(is.null(FTx)) NULL else dim(FTx)[2]
  #Parametric fit to obtain priors for theta
  likelihoodoptim<-function(pa.input){
    #pa.input includes  
    #th1, th2--parameters of the centering distribution
    #FTcoef--coefficients for FTx; could be NULL
    #crcoef--coefficients for crx
    parameters = list(theta =pa.input[1:2])
    parameters$crcoef = if(is.null(pcr)) NULL else pa.input[3:(pcr+2)]
    templength = if(is.null(pcr)) 2 else pcr+2
    parameters$FTcoef = if(is.null(pFT)) NULL else pa.input[(templength+1):(templength+pFT)]
    parameters$bz = c(0,0)
    parameters$Fv = rep(0,length(t1))
    parameters$crv = rep(0,length(t1))
    p = likelihood.wrapper(model,BP = 0,distr,data,parameters)$likelihoodsum
    return(-p);
  }
  
  parastart = if(is.null(pFT)) c(th.initial,rep(0,pcr)) else c(th.initial,rep(0,pcr),rep(0,pFT))
  fit = optim(par=parastart,likelihoodoptim,hessian="TRUE")
  
  #--------------------------------------------------------------------#
  #initial values
  #-------------------------------------------------------------------------------------------------------------#
  pmean.th = fit$par[1:2]; pcov.th = 0.01*fit$hessian[1:2,1:2]
  crcoef  = fit$par[3:(pcr+2)]
  FTcoef  = if (is.null(pFT)) NULL else fit$par[(pcr+3):(pcr+2+pFT)]
  bz  = rep(0,Jw-1)
  if(SR==1|SR==2|SR==3){sigma  = 0.2}
  if(SR==4){sigma = c(0.2,0.2)}
  if(SR==0){sigma = NULL;phi  = NULL}
  if(SR==1){phi = 0.3}else{phi=NULL}
  if(BP>0){alpha  = 5}else{alpha=1}
  
  #SR random effects indicators:
  #0--no random effects
  #1--Fv = rho*crv
  #2--crv, no Fv
  #3--Fv, no crv
  #4--both crv and Fv
  
  
  Fv  = if (SR==1|SR==3|SR==4) rnorm(m,0,0.01) else rep(0,m)
  crv  = if (SR==1|SR==2|SR==4) rnorm(m,0,0.01) else rep(0,m)
  
  return(list(theta = pmean.th, theta.mean = pmean.th, theta.cov = pcov.th, crcoef = crcoef, 
              FTcoef = FTcoef, bz = bz, sigma = sigma, phi = phi, 
              logalpha = log(alpha),Fv= Fv, crv = crv))
}


adaptive.v.mean.cov <-function(d,samplesize,proposal.old,new.x,smalln, adaptive.c){
  new.mean = recursivemean_vector(proposal.old$mean, new.x, samplesize)
  new.cov = recursivecov_vector(d, samplesize,proposal.old$cov,proposal.old$mean,new.x,smalln, adaptive.c/d)
  return(list(mean = new.mean, cov = new.cov))
}

adaptive.mean.var <-function(samplesize,proposal.old,new.x,smalln, adaptive.c){
    new.mean = recursivemean(proposal.old$mean, new.x, samplesize)
    new.var = recursivecov(samplesize,proposal.old$var,proposal.old$mean,new.x,smalln, adaptive.c)
    return(list(mean = new.mean, var = new.var))
  }
  
update.wrapper<-function(likelihood.n,likelihood.o,prior.o,prior.n,para.new,para.current){
  if(log(runif(1))<(sum(likelihood.n)+prior.n-sum(likelihood.o)-prior.o)){
    para.updated= para.new; likelihood.updated = likelihood.n; accept = TRUE
  }else{
    para.updated = para.current; likelihood.updated = likelihood.o;accept = FALSE
  }
  return(list(likelihood.updated =likelihood.updated ,para.updated = para.updated,accept = accept))
}
mcmc<-function(model = "PH",BP=1,SR=1,distr=3,data, mcmc.setup,BP.setup,th.initial){
  
  #model---"AFT" or "PH" or "PO"
  #BP--0: no Bernstein Polynomial modeling of the baseline hazard function; 1: otherwise
  #SR random effects indicators:
        #0--no random effects
        #1--Fv = rho*crv
        #2--crv, no Fv
        #3--Fv, no crv
        #4--both crv and Fv
  #distr--1:Logistic distribution; 2: log-Normal distribution; 3: Weibull distribution
 
  #data includes:
           #t1--lower bound of the observed interval for failure time or observed; default is zero
           #t2--upper bound of the observed interval for failure time or observed; default is Inf
           #ci -- cluster indicator
           #type--1: failure time is observed; 2: failure time is left censored; 3: failure time is
                     #right censored; 4: failure time is interval-censored
           #FTx--covariates for the latent failure time distribution; could be NULL
           #crx--covariates for the cure rate 
  
  #mcmc.setup includes:
           #nrun--total number of interations
           #nskip--number of skipped interations
           #nburn--number of burn-in interations
  #BP setup when a Bernstein Polynomial prior (BP) is assumed 
           #Jw--number of basis functions
           #a.alpha and b.alpha--hyperparameters for alpha
  
  #th.initial--initial parameters for the centering distribution
  
  t1 = data$t1
  t2 = data$t2
  type = data$type
  FTx = data$FTx
  crx = data$crx
  ci = data$ci
  m = length(unique(ci)) # number of clusters
  
  
  nrun = mcmc.setup$nrun
  nburn = mcmc.setup$nburn
  nskip = mcmc.setup$nskip
  
  Jw = BP.setup$Jw
  a.alpha = BP.setup$a.alpha
  b.alpha = BP.setup$b.alpha
  
  pcr = if (is.null(crx)) NULL else dim(crx)[2]
  pFT = if(is.null(FTx)) NULL else dim(FTx)[2]
  
  #-------------------------------------------------------------------------------------------------------------#
  # Adaptive MCMC
  
  para.current = mcmc.init (model, distr, SR,data,Jw, th.initial)
  
  
  n.I = 100; smalln = 1.0e-6; c.n = 1000
  theta.prop = list(mean =rep(0,2), cov = c.n*smalln*diag(2)); thchain.I = array(0,dim=c(n.I,2)); 
  if(!is.null(pcr)){crcoef.prop = list(mean = rep(0,pcr), cov = c.n*smalln*diag(pcr)); crcoefchain.I = array(0,dim=c(n.I,pcr))}
  if(!is.null(pFT)){FTcoef.prop = list(mean = rep(0,pFT), cov = c.n*smalln*diag(pFT)); FTcoefchain.I = array(0,dim=c(n.I,pFT));}
  if(BP>0){bz.prop = list(mean = rep(0,Jw-1), cov = c.n*smalln*diag(Jw-1));  bzchain.I = array(0,dim=c(n.I,Jw-1)); 
          logalpha.prop=list(mean = 0, var = c.n*smalln); alphachain.I = rep(0,n.I)
           phi.prop=list(mean = 0, var = c.n*smalln); phichain.I = rep(0,n.I)
  }
  

  if(SR==1|SR==2|SR==4|SR==0){crv.prop=list(mean = rep(0,m), var = c.n*rep(smalln,m)); crvchain.I = array(0,dim=c(n.I,m));}
  if(SR==3|SR==4){Fv.prop=list(mean = rep(0,m), var = c.n*rep(smalln,m)); FTvchain.I = array(0,dim=c(n.I,m));}
  #0--no random effects
  #1--Fv = rho*crv
  #2--crv, no Fv
  #3--Fv, no crv
  #4--both crv and Fv
   
  
  # Things to save
  indsave = 0; 
  nsave = nrun/nskip - nburn
  if(!is.null(pcr)){crcoefchain = array(0,dim=c(nsave,pcr));} 
  if(!is.null(pFT)){FTcoefchain = array(0,dim=c(nsave,pFT));} 
  thetachain = array(0,dim=c(nsave,2)); 
  if(BP>0){weightchain = array(0,dim=c(nsave,Jw));
           bzchain = array(0,dim=c(nsave,Jw-1)); 
           alphachain = rep(0,nsave)
  }else{weightchain = NULL; bzchain = NULL; alphachain = NULL}
  if(SR==1|SR==2|SR==4){
           crvchain = array(0,dim=c(nsave,m)); 
  }else{crvchain = NULL}
  if(SR==3|SR==4){
    FTvchain = array(0,dim=c(nsave,m)); 
  }else{FTvchain = NULL}
  
  if(SR==1){phichain = rep(0,nsave)}else{phichain = NULL}
  if(SR==0){sigmauchain = NULL}
  if(SR==1|SR==2|SR==3){sigmauchain = rep(0,nsave)}
  if(SR==4){sigmauchain = array(0,dim=c(nsave,2))}
  
  acceptance = list(theta = 0, FTcoef = 0, crcoef = 0, bz = 0, alpha = 0, phi = 0, sigmau = 0, crv = rep(0,m), Fv = rep(0,m))
  likelihoodchain =array(0,dim=c(nsave,m));
  
  
  #MCMC iterations
  for(iscan in 2:nrun){
    adaptive.c = 2.4^2
    
    if(iscan==n.I){
      n.eff = iscan -2
      theta.prop$mean = apply(thchain.I[1:(n.eff),], 2,mean)
      theta.prop$cov = adaptive.c/2*cov(thchain.I[1:(n.eff),])+smalln*diag(2)
      if(!is.null(pcr)){crcoef.prop$mean = apply(crcoefchain.I[1:(n.eff),], 2,mean); crcoef.prop$cov = adaptive.c/pcr*cov(crcoefchain.I[1:(n.eff),])+smalln*diag(pcr)}
      if(!is.null(pFT)){FTcoef.prop$mean = apply(FTcoefchain.I[1:(n.eff),], 2,mean); FTcoef.prop$cov = adaptive.c/pFT*cov(FTcoefchain.I[1:(n.eff),])+smalln*diag(pFT)}
      if(BP>0){
          bz.prop$mean = apply(bzchain.I[1:(n.eff),], 2,mean)
          bz.prop$cov = adaptive.c/(Jw-1)*cov(bzchain.I[1:(n.eff),])+smalln*diag(Jw-1)
          logalpha.prop$mean = mean(log(alphachain.I[1:(n.eff)]))
          logalpha.prop$var= adaptive.c*var(log(alphachain.I[1:(n.eff)]))+smalln
      }
      if(SR==1|SR==2|SR==4){
          if(SR==1){phi.prop$mean = mean(phichain.I[1:(n.eff)])
          phi.prop$var = adaptive.c*var(phichain.I[1:(n.eff)])+smalln}
          crv.prop$mean = apply(crvchain.I[1:n.eff,],2,mean)
          crv.prop$var = adaptive.c*apply(crvchain.I[1:n.eff,],2,var)+rep(smalln,m)}
    
    if(SR==3|SR==4){
      Fv.prop$mean = apply(FTvchain.I[1:n.eff,],2,mean)
      Fv.prop$var = adaptive.c*apply(FTvchain.I[1:n.eff,],2,var)+rep(smalln,m)
      }
    }
    
    if(iscan>n.I){
   
      n.eff = iscan -2
      theta.prop = adaptive.v.mean.cov(2,n.eff,theta.prop,para.current$theta,smalln,adaptive.c)
      if(is.null(pcr)){crcoef.prop = adaptive.v.mean.cov(pcr,n.eff,crcoef.prop,para.current$crcoef,smalln,adaptive.c) }
      if(is.null(pFT)){FTcoef.prop = adaptive.v.mean.cov(pFT,n.eff,FTcoef.prop,para.current$FTcoef,smalln,adaptive.c) }
      if(BP>0){
       bz.prop = adaptive.v.mean.cov(Jw-1,n.eff,bz.prop,para.current$bz,smalln,adaptive.c) 
       logalpha.prop = adaptive.mean.var(n.eff,logalpha.prop,para.current$logalpha,smalln,adaptive.c) 
      }
      if(SR==1|SR==2|SR==4){
      if(SR==1){phi.prop = adaptive.mean.var(n.eff,phi.prop,para.current$phi,smalln,adaptive.c)} 
      for(i in 1:m){
        crvi.prop = list(mean = crv.prop$mean[i],var =  crv.prop$var[i])
        crvi.prop = adaptive.mean.var(n.eff,crvi.prop,para.current$crv[i],smalln,adaptive.c) 
        crv.prop$mean[i]=crvi.prop$mean
        crv.prop$var[i]=crvi.prop$var
      }}
      if(SR==3|SR==4){
        Fvi.prop = list(mean = Fv.prop$mean[i],var =  Fv.prop$var[i])
        Fvi.prop = adaptive.mean.var(n.eff,Fvi.prop,para.current$Fv[i],smalln,adaptive.c) 
        Fv.prop$mean[i]=Fvi.prop$mean
        Fv.prop$var[i]=Fvi.prop$var
      }
      }
    
   
    likelihood.c = likelihood.wrapper(model,BP,distr,data,para.current)$likelihoodsum
                              
    #update theta
    #Maybe updating theta with bz?
    para.new = para.current
    prior.c = log_dnorm (para.current$theta, para.current$theta.mean, para.current$theta.cov,2)
    para.new$theta = rmnorm(para.current$theta,theta.prop$cov)
    likelihood.n = likelihood.wrapper(model,BP,distr,data,para.new)$likelihoodsum
    prior.n = log_dnorm (para.new$theta, para.current$theta.mean, para.current$theta.cov,2)
    theta.result = update.wrapper(likelihood.n,likelihood.c,prior.c,prior.n,para.new,para.current)

    likelihood.c = theta.result$likelihood.updated
    para.current = theta.result$para.updated
    acceptance$theta = 1*isTRUE(theta.result$accept)+acceptance$theta
    
   
    if(BP>0){
      #update weight parameters bz
      para.new = para.current
      prior.c = exp(para.current$logalpha)*sum(log(Ys_to_weight(para.current$bz)))
      para.new$bz = rmnorm(para.current$bz,bz.prop$cov)
      likelihood.n = likelihood.wrapper(model,BP,distr,data,para.new)$likelihoodsum
      prior.n = exp(para.current$logalpha)*sum(log(Ys_to_weight(para.new$bz)))
      bz.result = update.wrapper(likelihood.n,likelihood.c,prior.c,prior.n,para.new,para.current)
      likelihood.c = bz.result$likelihood.updated
      para.current = bz.result$para.updated
      acceptance$bz = 1*isTRUE(bz.result$accept)+acceptance$bz
      #update alpha
      #browser()
      weight.c = Ys_to_weight(para.current$bz)
      alpha.o = exp(para.current$logalpha)
      prior.alpha.o = lgamma(alpha.o*Jw)-Jw*lgamma(alpha.o)+sum((alpha.o-1)*log(weight.c))+(a.alpha-1)*log(alpha.o)-b.alpha*alpha.o
      alpha.n = exp(rnorm(1,para.current$logalpha,sqrt( logalpha.prop$var)))
      prior.alpha.n = lgamma(alpha.n*Jw)-Jw*lgamma(alpha.n)+sum((alpha.n-1)*log(weight.c))+(a.alpha-1)*log(alpha.n)-b.alpha*alpha.n
      
      if(log(runif(1))<(prior.alpha.n-prior.alpha.o)){
        para.current$logalpha = log(alpha.n); acceptance$alpha = acceptance$alpha+1
      }
      
    }
   
    #update crcoef parameters 
    if(!is.null(pcr)){
    para.new = para.current
    prior.c = -t(para.current$crcoef)%*%(para.current$crcoef)/5^2/2
    para.new$crcoef = if(pcr==1) rnorm(1,para.current$crcoef,sqrt(crcoef.prop$cov)) else rmnorm(para.current$crcoef,crcoef.prop$cov)
    likelihood.n = likelihood.wrapper(model,BP,distr,data,para.new)$likelihoodsum
    prior.n = -t(para.new$crcoef)%*%(para.new$crcoef)/5^2/2
    crcoef.result = update.wrapper(likelihood.n,likelihood.c,prior.c,prior.n,para.new,para.current)
    likelihood.c = crcoef.result$likelihood.updated
    para.current = crcoef.result$para.updated
    acceptance$crcoef = 1*isTRUE(crcoef.result$accept)+acceptance$crcoef 
    }
    
    
    #update FTcoef parameters 
    if(!is.null(pFT)){
      para.new = para.current
      prior.c = -t(para.current$FTcoef)%*%(para.current$FTcoef)/5^2/2
      para.new$FTcoef = if(pFT==1) rnorm(1,para.current$FTcoef,sqrt(FTcoef.prop$cov)) else rmnorm(para.current$FTcoef,FTcoef.prop$cov)
      likelihood.n = likelihood.wrapper(model,BP,distr,data,para.new)$likelihoodsum
      prior.n = -t(para.new$FTcoef)%*%(para.new$FTcoef)/5^2/2
      FTcoef.result = update.wrapper(likelihood.n,likelihood.c,prior.c,prior.n,para.new,para.current)
      likelihood.c = FTcoef.result$likelihood.updated
      para.current = FTcoef.result$para.updated
      acceptance$FTcoef = isTRUE(FTcoef.result$accept)+acceptance$FTcoef
    }
    
    if(SR==1){
    #update random effects
    likelihood.cv = likelihood.wrapper(model,BP,distr,data,para.current)
    para.new = para.current
    para.new$crv = samplerancrv(para.current$crv, crv.prop$var, m)
    para.new$FTv = para.new$phi* para.new$crv
    likelihood.nv = likelihood.wrapper(model,BP,distr,data,para.new)
    prior.c = -para.current$crv^2/para.current$sigma^2/2
    prior.n = -para.new$crv^2/para.new$sigma^2/2
    eval.v = (prior.n+likelihood.nv$likelihoodi)>(prior.c+likelihood.cv$likelihoodi)
    para.current$crv = 1*(eval.v)*para.new$crv + (1-1*(eval.v))*para.current$crv
    acceptance$crv = 1*(eval.v)+acceptance$crv
    
    
    #update sigma.u
    Lambda= sampleLambda(a = 2, zetasq = 10^3, para.current$sigma)
    para.current$sigma = samplesigma(n, a=2, para.current$crv, Lambda)
    
    #update phi
    para.new = para.current
    para.new$phi = rnorm(1,para.current$phi,sqrt(phi.prop$var))
    para.new$FTv = para.new$phi*para.new$crv;
    prior.c = - para.current$phi^2/25/2
    prior.n = - para.new$phi^2/25/2
    likelihood.c = likelihood.wrapper(model,BP,distr,data,para.current)$likelihoodsum 
    likelihood.n = likelihood.wrapper(model,BP,distr,data,para.new)$likelihoodsum 
    phi.result = update.wrapper(likelihood.n,likelihood.c,prior.c,prior.n,para.new,para.current)
    likelihood.c = phi.result$likelihood.updated
    para.current = phi.result$para.updated
    para.current$FTv = para.current$phi*para.current$crv
    acceptance$phi= isTRUE(phi.result$accept)+acceptance$phi
    }
    
    if(SR==2){
      #update random effects
      likelihood.cv = likelihood.wrapper(model,BP,distr,data,para.current)
      para.new = para.current
      para.new$crv = samplerancrv(para.current$crv, crv.prop$var, m)
      para.new$FTv = rep(0,m)
      likelihood.nv = likelihood.wrapper(model,BP,distr,data,para.new)
      prior.c = -para.current$crv^2/para.current$sigma^2/2
      prior.n = -para.new$crv^2/para.new$sigma^2/2
      eval.v = (prior.n+likelihood.nv$likelihoodi)>(prior.c+likelihood.cv$likelihoodi)
      para.current$crv = 1*(eval.v)*para.new$crv + (1-1*(eval.v))*para.current$crv
      acceptance$crv = 1*(eval.v)+acceptance$crv
      
      #update sigma.u
      Lambda= sampleLambda(a = 2, zetasq = 10^3, para.current$sigma)
      para.current$sigma = samplesigma(n, a=2, para.current$crv, Lambda)
    }
    
    if(SR==3){
      #update random effects
      likelihood.cv = likelihood.wrapper(model,BP,distr,data,para.current)
      para.new = para.current
      para.new$Fv = samplerancrv(para.current$Fv, Fv.prop$var, m)
      para.new$crv = rep(0,m)
      likelihood.nv = likelihood.wrapper(model,BP,distr,data,para.new)
      prior.c = -para.current$Fv^2/para.current$sigma^2/2
      prior.n = -para.new$Fv^2/para.new$sigma^2/2
      eval.v = (prior.n+likelihood.nv$likelihoodi)>(prior.c+likelihood.cv$likelihoodi)
      para.current$Fv = 1*(eval.v)*para.new$Fv + (1-1*(eval.v))*para.current$Fv
      acceptance$Fv = 1*(eval.v)+acceptance$Fv
      
      #update sigma.u
      Lambda= sampleLambda(a = 2, zetasq = 10^3, para.current$sigma)
      para.current$sigma = samplesigma(n, a=2, para.current$Fv, Lambda)
    }
    
    if(SR==4){
      #update random effects
      
      likelihood.cv = likelihood.wrapper(model,BP,distr,data,para.current)
      para.new = para.current
      para.new$crv = samplerancrv(para.current$crv, crv.prop$var, m)
      para.new$FTv = rep(0,m)
      likelihood.nv = likelihood.wrapper(model,BP,distr,data,para.new)
      prior.c = -para.current$crv^2/para.current$sigma[1]^2/2
      prior.n = -para.new$crv^2/para.new$sigma[1]^2/2
      eval.v = (prior.n+likelihood.nv$likelihoodi)>(prior.c+likelihood.cv$likelihoodi)
      para.current$crv = 1*(eval.v)*para.new$crv + (1-1*(eval.v))*para.current$crv
      acceptance$crv = 1*(eval.v)+acceptance$crv
      
      #update sigma.u
      Lambdacrv= sampleLambda(a = 2, zetasq = 10^3, para.current$sigma[2])
      para.current$sigma[1] = samplesigma(n, a=2, para.current$crv, Lambdacrv)
      
      
      likelihood.cv = likelihood.wrapper(model,BP,distr,data,para.current)
      para.new = para.current
      para.new$Fv = samplerancrv(para.current$Fv, Fv.prop$var, m)
      para.new$crv = rep(0,m)
      likelihood.nv = likelihood.wrapper(model,BP,distr,data,para.new)
      prior.c = -para.current$Fv^2/para.current$sigma[2]^2/2
      prior.n = -para.new$Fv^2/para.new$sigma[2]^2/2
      eval.v = (prior.n+likelihood.nv$likelihoodi)>(prior.c+likelihood.cv$likelihoodi)
      para.current$Fv = 1*(eval.v)*para.new$Fv + (1-1*(eval.v))*para.current$Fv
      acceptance$Fv = 1*(eval.v)+acceptance$Fv
      
      #update sigma.u
      LambdaFv= sampleLambda(a = 2, zetasq = 10^3, para.current$sigma[2])
      para.current$sigma[2] = samplesigma(n, a=2, para.current$Fv, LambdaFv)
    }
    
    
    
    #print(c(para.current$FTcoef,para.current$crcoef,para.current$logalpha,para.current$phi))
    if(iscan<=n.I){
      crcoefchain.I[iscan-1,] = para.current$crcoef; 
      FTcoefchain.I[iscan-1,] = para.current$FTcoef; 
      thchain.I[iscan-1,] = para.current$theta; 
      if(BP>0){bzchain.I[iscan-1,] = para.current$bz; 
      alphachain.I[iscan-1] = exp(para.current$logalpha)}
      if(SR==1){phichain.I[iscan-1] = para.current$phi} 
      if(SR==1|SR==2|SR==4){crvchain.I[iscan-1,] = para.current$crv;} 
      if(SR==3|SR==4){FTvchain.I[iscan-1,] = para.current$Fv;} 
      
    }
    
    if(iscan%%nskip==0 && iscan>nskip*nburn){
      indsave = indsave +1
      #print(cbind(iscan,indsave,sum(likelihood.s),alpha.c,phi.c,acceptance[1]))
      crcoefchain[indsave,] = para.current$crcoef; 
      FTcoefchain[indsave,] = para.current$FTcoef; 
      thetachain[indsave,] = para.current$theta; 
      if(BP>0){weightchain[indsave,] = Ys_to_weight(para.current$bz); 
      bzchain[indsave,] = para.current$bz
      alphachain[indsave] = exp(para.current$logalpha)}
      if(SR==1){phichain[indsave] = para.current$phi} 
      if(SR==1|SR==2|SR==3){sigmauchain[indsave] = para.current$sigma}
      if(SR==4){sigmauchain[indsave,] = para.current$sigma}
      if(SR==1|SR==2|SR==4){crvchain[indsave,] = para.current$crv}
      if(SR==3|SR==4){FTvchain[indsave,] = para.current$Fv}
      
      likelihoodchain[indsave,] = likelihood.c
    }
    
  }
  #Compute WAIC
  lppd = sum(log(apply(exp(likelihoodchain),2,mean)))
  WAIC1 = lppd-(2*sum(log(apply(exp(likelihoodchain),2,mean))-apply(likelihoodchain,2,mean)))
  WAIC2 = lppd-sum(apply(likelihoodchain,2,var))#preferred
  LPML = -sum(log(apply(exp(-likelihoodchain),2,mean)))
  #browser()
  
  result<-list(theta = thetachain,
               crcoef = crcoefchain,
               FTcoef = FTcoefchain,
               weight = weightchain,
               bz = bzchain,
               alpha = alphachain,
               crv = crvchain,
               FTv = FTvchain,
               likelihood = likelihoodchain,
               sigma = sigmauchain,
               phi = phichain,
               acceptance = acceptance,
               lppd = lppd,
               WAIC1 = WAIC1,
               WAIC2 = WAIC2,
               LPML = LPML);
  return(result);
}



