#Gwishart random effect cov
#dimFTcoef--indicator whether there are covariates for the survival part

Rcpp::sourceCpp("~/Dropbox/research/Active/Dipankar/Analysis/Data/BP/Ht/MCMC_BP_multivariateRE.cpp")

produ<-function(x,beta){
  if(length(beta)>1){return(x%*%beta)}
  else{return(x*beta)}
}
mcmc<-function(model,distr,maxc,t1,t2,type,m,d_r,BP,SR,crx,FTx,dimFTcoef,Randomx,nrun,nskip,nskip_r,nburn,nburn_r,Jw,a_alpha,b_alpha,th_initial){
  if(is.vector(crx)){pcr = 1}
  else{pcr = length(crx[1,])}
  if(is.vector(FTx)){pFT = 1}
  else{  pFT = length(FTx[1,])}
 #Parametric fit to obtain priors for theta
  likelihoodoptim<-function(parameters){
    th1_p = parameters[1] 
    th2_p = parameters[2]
    crcoef_p=parameters[3:(pcr+2)] 
    if(dimFTcoef>0){FTcoef_p=parameters[(pcr+3):(pcr+2+pFT)]
    FXbeta = produ(FTx,FTcoef_p)}
    else{FXbeta = rep(0,length(t1))}
    w = rep(1/Jw,Jw)
    crXbeta = produ(crx,crcoef_p)
    Fv = rep(0,length(t1))
    crv = rep(0,length(t1))
    p=likelihoodv_para(model, 0,m, t1, t2, type, th1_p, th2_p, w, distr, FXbeta, 
                       Fv, crXbeta, crv) 
    return(-sum(p));
  }
  if(dimFTcoef>0){parastart=c(th_initial,rep(0,pcr),rep(0,pFT))}
  else{parastart=c(th_initial,rep(0,pcr))}
  fit = optim(par=parastart,likelihoodoptim,hessian="TRUE")
  pmean_th = fit$par[1:2]
  if(BP>0){pcov_th = fit$hessian[1:2,1:2]}
  else{pcov_th = 0.001*fit$hessian[1:2,1:2]}
  
  #--------------------------------------------------------------------#
  #initial values
  #-------------------------------------------------------------------------------------------------------------#
  th_c = fit$par[1:2]
  crcoef_c = fit$par[3:(pcr+2)]
  if(dimFTcoef>0){FTcoef_c = fit$par[(pcr+3):(pcr+2+pFT)]}
  else{FTcoef_c = 0}
  bz_c = rep(0,Jw-1)
  Sigma_c = diag(0.2,d_r)
  phi_c = 0.3
  alpha_c = 0.5
  if(SR>0){randomeffect_c=matrix(rnorm(n*d_r,0,0.01),ncol=n,byrow=TRUE)}
  else{randomeffect_c=matrix(rep(0,n*d_r),ncol=n,byrow=TRUE)}
  
  #-------------------------------------------------------------------------------------------------------------#
  # Adaptive MCMC
  indsave = 0; indsave_r=0
  mean_th =rep(0,2); cov_th = matrix(0,ncol=2,nrow=2)
  if(pcr>1){mean_crcoef = rep(0,pcr); cov_crcoef = matrix(0,ncol=pcr,nrow=pcr)}
  else{  mean_crcoef = 0; cov_crcoef = 0}
  if(pFT>1){mean_FTcoef = rep(0,pFT); cov_FTcoef = matrix(0,ncol=pFT,nrow=pFT)}
  else{mean_FTcoef = 0; cov_FTcoef = 0}
  mean_bz = rep(0,Jw-1); cov_bz = matrix(0,ncol=Jw-1,nrow=Jw-1)
  mean_logalpha = 0; var_logalpha = 0
  mean_phi = 0; var_phi = 0
  mean_randomeffect = rep(0,n*d_r);cov_randomeffect = array(0,dim=c(d_r,d_r,n))
  
  # Things to save
  nsave = nrun/nskip - nburn
  n = length(t1)/m
  crcoefchain = array(0,dim=c(nsave,pcr)); 
  FTcoefchain = array(0,dim=c(nsave,pFT)); 
  thchain = array(0,dim=c(nsave,2)); 
  weightchain = array(0,dim=c(nsave,Jw));
  bzchain = array(0,dim=c(nsave,Jw-1)); 
  alphachain = rep(0,nsave)
  phichain = rep(0,nsave)  
  likelihoodchain =array(0,dim=c(nsave,n)); 
  Sigmauchain = array(0,dim=c(nsave,d_r,d_r))
  randomeffectchain = array(0,dim=c(nrun/nskip_r-nburn_r,d_r,n)); 
  acceptance = rep(0,n+6)
  
  # For computing initial covariance array
  n_I = 100
  crcoefchain_I = array(0,dim=c(n_I,pcr)); 
  FTcoefchain_I = array(0,dim=c(n_I,pFT)); 
  thchain_I = array(0,dim=c(n_I,2)); 
  bzchain_I = array(0,dim=c(n_I,Jw-1)); 
  alphachain_I = rep(0,n_I)
  phichain_I = rep(0,n_I) 
  randomeffectchain_I = array(0,dim=c(n_I,d_r,n));
  
  #MCMC iterations
  for(iscan in 2:nrun){
    
    th_o = th_c
    crcoef_o = crcoef_c
    FTcoef_o = FTcoef_c
    bz_o = bz_c
    alpha_o = alpha_c
    phi_o = phi_c
    cru_o = c(randomeffect_c)
    FTw_o = phi_c*cru_o
    Sigma_o = Sigma_c
    
    smalln = 1e-6
    if(iscan<n_I){
      cov_th = 1000*smalln*diag(2)
      cov_crcoef = 1000*smalln*diag(pcr)
      cov_FTcoef = 1000*smalln*diag(pFT)
      cov_bz = 1000*smalln*diag(Jw-1)
      var_logalpha= 0.01
      var_phi = 0.001
      for(l in 1:n){
         cov_randomeffect[,,l] = 1000*smalln*diag(d_r)
      }
    }
    if(iscan==n_I){
      n_eff = iscan -2
      mean_th = apply(thchain_I[1:(n_eff),], 2,mean)
      cov_th = 2.4^2/2*cov(thchain_I[1:(n_eff),])+smalln*diag(2)
      if(pcr>1){mean_crcoef = apply(crcoefchain_I[1:(n_eff),], 2,mean)
      cov_crcoef = 2.4^2/pcr*cov(crcoefchain_I[1:(n_eff),])+smalln*diag(pcr)}
      else{
        mean_crcoef = mean(crcoefchain_I[1:(n_eff)])
        cov_crcoef = 2.4^2*var(crcoefchain_I[1:(n_eff)])+smalln
      }
      if(pFT>1){mean_FTcoef = apply(FTcoefchain_I[1:(n_eff),], 2,mean)
      cov_FTcoef = 2.4^2/pFT*cov(FTcoefchain_I[1:(n_eff),])+smalln*diag(pFT)}
      else{
        mean_FTcoef = mean(FTcoefchain_I[1:(n_eff)])
        cov_FTcoef = 2.4^2*var(FTcoefchain_I[1:(n_eff)])+smalln
      }
      mean_bz = apply(bzchain_I[1:(n_eff),], 2,mean)
      cov_bz = 2.4^2/(Jw-1)*cov(bzchain_I[1:(n_eff),])+smalln*diag(Jw-1)
      mean_logalpha = mean(log(alphachain_I[1:(n_eff)]))
      var_logalpha= 2.4^2*var(log(alphachain_I[1:(n_eff)]))+smalln
      
      if(SR>0){
      mean_phi = mean(phichain_I[1:(n_eff)])
      var_phi = 2.4^2*var(phichain_I[1:(n_eff)])+smalln
      
      for(l in 1:n){
        seq = ((l-1)*d_r+1):(l*d_r)
        if(d_r>1){mean_randomeffect[seq] = apply(randomeffectchain_I[1:n_eff,,l],2,mean)
        cov_randomeffect[,,l] = 2.4^2/d_r*cov(randomeffectchain_I[1:n_eff,,l])+smalln*diag(d_r)}
    }
    }}
    
    if(iscan>n_I){
      n_eff = iscan -2
      cov_th = recursivecov_vector(2,n_eff  ,cov_th,mean_th,th_c,smalln,2.4^2/2)
      mean_th = recursivemean_vector(mean_th,th_c,n_eff ) 
      if(pcr>1){
      cov_crcoef = recursivecov_vector(pcr,n_eff ,cov_crcoef,mean_crcoef,crcoef_c,smalln,2.4^2/pcr)
      mean_crcoef = recursivemean_vector(mean_crcoef,crcoef_c,n_eff)
      }
      else{
      cov_crcoef = recursivecov(n_eff ,cov_crcoef,mean_crcoef,crcoef_c,smalln,2.4^2/pcr)
      mean_crcoef = recursivemean(mean_crcoef,crcoef_c,n_eff)
      }
      if(pFT>1){cov_FTcoef = recursivecov_vector(pFT,n_eff ,cov_FTcoef,mean_FTcoef,FTcoef_c,smalln,2.4^2/pFT)
      mean_FTcoef = recursivemean_vector(mean_FTcoef,FTcoef_c,n_eff)
      }
      else{
        cov_FTcoef = recursivecov(n_eff ,cov_FTcoef,mean_FTcoef,FTcoef_c,smalln,2.4^2/pFT)
        mean_FTcoef = recursivemean(mean_FTcoef,FTcoef_c,n_eff)
      }
      cov_bz = recursivecov_vector(Jw-1,n_eff ,cov_bz,mean_bz,bz_c,smalln,2.4^2/(Jw-1))
      mean_bz = recursivemean_vector(mean_bz,bz_c,n_eff)
      var_logalpha = recursivecov(n_eff ,var_logalpha,mean_logalpha,log(alpha_c),smalln,2.4^2)
      mean_logalpha = recursivemean(mean_logalpha,log(alpha_c),n_eff)
      if(SR>0){
      var_phi = recursivecov(n_eff ,var_phi,mean_phi,phi_c,smalln,2.4^2)
      mean_phi = recursivemean(mean_phi,phi_c,n_eff)
      cov_randomeffect = recursivecov_randomeffect(n,d_r,n_eff ,cov_randomeffect,mean_randomeffect,randomeffect_c,smalln,2.4^2/d_r)
      mean_randomeffect = recursivemean_vector(mean_randomeffect,randomeffect_c,n_eff)
    }}
    
    if(dimFTcoef>0){FXbeta_o = produ(FTx,FTcoef_o)}
    else{FXbeta_o=rep(0,n*m)}
    crXbeta_o = produ(crx,crcoef_o)
    z_o = bz_o
    weight_o = Ys_to_weight(z_o)
    
    #cru_o = cru_o-mean(cru_o) #mean constraint
    crv_o = matrixproduct( Randomx, cru_o, n, m, d_r);Fv_o = phi_o*crv_o;
    missingM =  missing_impute_t( model, BP, distr,maxc,t1, t2,type, th_o[1], th_o[2], weight_o, FXbeta_o,Fv_o,crXbeta_o,crv_o)
    
    #update theta
    likelihood_o = likelihoodv(model, BP,m, t1, t2, type, th_o[1], th_o[2], weight_o, distr, FXbeta_o, Fv_o, crXbeta_o, crv_o, missingM) 
    prior_th_o = log_dnorm (th_o, pmean_th, pcov_th,2)
    th_n = rmnorm(th_o,cov_th)
    likelihood_n = likelihoodv(model, BP,m, t1, t2, type, th_n[1], th_n[2], weight_o, distr, FXbeta_o, Fv_o, crXbeta_o, crv_o, missingM) 
    prior_th_n = log_dnorm (th_n, pmean_th, pcov_th,2)
    
    if(log(runif(1))<(sum(likelihood_n)+prior_th_n-sum(likelihood_o)-prior_th_o)){
      th_c = th_n; likelihood_s = likelihood_n; acceptance[1] = acceptance[1]+1
    }
    else{
      th_c = th_o; likelihood_s = likelihood_o;
    }
    
    if(BP>0){
      #update weight parameters bz
      likelihood_o = likelihood_s
      prior_weight_o = alpha_o*sum(log(weight_o))
      bz_n = rmnorm(bz_o,cov_bz)
      z_n = bz_n;weight_n = Ys_to_weight(z_n)
      likelihood_n = likelihoodv(model, BP,m, t1, t2, type, th_c[1], th_c[2], weight_n, distr, FXbeta_o, Fv_o, crXbeta_o, crv_o, missingM) 
      prior_weight_n = alpha_o*sum(log(weight_n))
      
      if(log(runif(1))<(sum(likelihood_n)+prior_weight_n-sum(likelihood_o)-prior_weight_o)){
        bz_c = bz_n; weight_c = weight_n; likelihood_s = likelihood_n; acceptance[2] = acceptance[2]+1
      }
      else{
        bz_c = bz_o; weight_c = weight_o; likelihood_s = likelihood_o;
      }
      
      #update alpha
      z_c = bz_c
      weight_c = Ys_to_weight(z_c)
      prior_alpha_o = lgamma(alpha_o*Jw)-Jw*lgamma(alpha_o)+sum((alpha_o-1)*log(weight_c))+(a_alpha-1)*log(alpha_o)-b_alpha*alpha_o
      alpha_n = exp(rnorm(1,log(alpha_o),sqrt(var_logalpha)))
      prior_alpha_n = lgamma(alpha_n*Jw)-Jw*lgamma(alpha_n)+sum((alpha_n-1)*log(weight_c))+(a_alpha-1)*log(alpha_n)-b_alpha*alpha_n
      if(log(runif(1))<(prior_alpha_n-prior_alpha_o)){
        alpha_c = alpha_n; acceptance[3] = acceptance[3]+1
      }
      else{
        alpha_c = alpha_o; 
      }
    }
    else{
      weight_c = rep(1/Jw,Jw);
    }

    #update crcoef parameters 
    likelihood_o = likelihood_s
    prior_crcoef_o = -t(crcoef_o)%*%(crcoef_o)/3^2/2
    if(pcr>1){crcoef_n = rmnorm(crcoef_o,cov_crcoef)}
    else{crcoef_n = rnorm(1,crcoef_o,sqrt(cov_crcoef))}
    crXbeta_n = produ(crx,crcoef_n)
    likelihood_n = likelihoodv(model, BP,m, t1, t2, type, th_c[1], th_c[2], weight_c, distr, FXbeta_o, Fv_o, crXbeta_n, crv_o, missingM) 
    prior_crcoef_n = -t(crcoef_n)%*%(crcoef_n)/3^2/2
    
    if(log(runif(1))<(sum(likelihood_n)+prior_crcoef_n[1]-sum(likelihood_o)-prior_crcoef_o[1])){
      crcoef_c = crcoef_n; crXbeta_c = crXbeta_n; likelihood_s = likelihood_n; acceptance[4] = acceptance[4]+1
    }
    else{
      crcoef_c = crcoef_o; crXbeta_c = crXbeta_o; likelihood_s = likelihood_o;
    }
   
    #update FTcoef parameters 
    if(dimFTcoef>0){
      likelihood_o = likelihood_s
      prior_FTcoef_o = -t(FTcoef_o)%*%(FTcoef_o)/10^2/2
      if(pFT>1){FTcoef_n = rmnorm(FTcoef_o,cov_FTcoef)}
      else{FTcoef_n = rnorm(1,FTcoef_o,sqrt(cov_FTcoef))}
      FXbeta_n = produ(FTx,FTcoef_n)
      likelihood_n = likelihoodv(model, BP,m, t1, t2, type, th_c[1], th_c[2], weight_c, distr, FXbeta_n, Fv_o, crXbeta_c, crv_o, missingM) 
      prior_FTcoef_n = -t(FTcoef_n)%*%(FTcoef_n)/10^2/2
      
      if(log(runif(1))<(sum(likelihood_n)+prior_FTcoef_n-sum(likelihood_o)-prior_FTcoef_o)){
        FTcoef_c = FTcoef_n; FXbeta_c = FXbeta_n;likelihood_s = likelihood_n; acceptance[5] = acceptance[5]+1
      }
      else{
        FTcoef_c = FTcoef_o; FXbeta_c = FXbeta_o;likelihood_s = likelihood_o;
      }
    }
    else{
      FTcoef_c = 0
      FXbeta_c = rep(0,n*m)
    }
  
    
    if(SR>0){
    #update random effects
    likelihood_o = likelihood_s
    invsigmastar = solve(Sigma_o)
    temp = samplerandomeffect(FTw_o, cru_o, cov_randomeffect,phi_o,d_r, n)
    
    FTw_n = temp[,1]; cru_n = temp[,2]
    crv_n = matrixproduct( Randomx, cru_n, n, m, d_r);Fv_n = phi_o*crv_n;
    likelihood_n = likelihoodv(model, BP,m, t1, t2, type, th_c[1], th_c[2], weight_c, distr, FXbeta_c, Fv_n, crXbeta_c, crv_n, missingM) 
    
    tempmatrix=randomeffect (cru_o,cru_n,likelihood_o,likelihood_n,invsigmastar,n,m,d_r)
    randomeffect_c=tempmatrix[1:d_r,]
    
    crummean = apply(randomeffect_c,1,mean)
    for(j in 1:d_r){
      randomeffect_c[j,] = randomeffect_c[j,] - crummean[j] #recenter
    }
    
    acceptance[5:(n+4)] = tempmatrix[d_r+1,]+acceptance[6:(n+5)]
    cru_c = c(randomeffect_c)
    crv_c = matrixproduct( Randomx, cru_c, n, m, d_r);
    Fv_oc = phi_o*crv_c;
    
    #update Sigma_u
    
    
    a_u=2
    zetasq_u=rep(100,d_r)
    Lambdau= sampleLambda(a_u, d_r, zetasq_u, Sigma_o)
    Sigma_c = samplesigma(n, a_u, d_r, t(randomeffect_c), Lambdau)
    
    
    
    #update phi
    phi_n = rnorm(1,phi_o,sqrt(var_phi))
    Fv_c = phi_n*crv_c;
    prior_phi_o = -phi_o^2/2/9
    prior_phi_n = -phi_n^2/2/9
    likelihood_o = likelihoodv(model, BP,m, t1, t2, type, th_c[1], th_c[2], weight_c, distr, FXbeta_c, Fv_oc, crXbeta_c, crv_c, missingM) 
    likelihood_n = likelihoodv(model, BP,m, t1, t2, type, th_c[1], th_c[2], weight_c, distr, FXbeta_c, Fv_c, crXbeta_c, crv_c, missingM)
    if(log(runif(1))<(sum(likelihood_n)+prior_phi_n-sum(likelihood_o)-prior_phi_o)){
      phi_c = phi_n;acceptance[n+6] = acceptance[n+6]+1;likelihood_s=likelihood_n
    }
    else{phi_c=phi_o;likelihood_s=likelihood_o}
    }
    else{
      randomeffect_c = matrix(rep(0,n*d_r),ncol=n,byrow=TRUE)
    }
    
    #print(cbind(iscan,indsave,sum(likelihood_s),alpha_c,acceptance[1]))
    if(iscan<=n_I){
      crcoefchain_I[iscan-1,] = crcoef_c; 
      FTcoefchain_I[iscan-1,] = FTcoef_c; 
      thchain_I[iscan-1,] = th_c; 
      bzchain_I[iscan-1,] = bz_c; 
      alphachain_I[iscan-1] = alpha_c
      phichain_I[iscan-1] = phi_c 
      randomeffectchain_I[iscan-1,,] = randomeffect_c; 
      
    }
    
    if(iscan%%nskip==0 && iscan>nskip*nburn){
      indsave = indsave +1
      #print(cbind(iscan,indsave,sum(likelihood_s),alpha_c,phi_c,acceptance[1]))
      crcoefchain[indsave,] = crcoef_c; 
      FTcoefchain[indsave,] = FTcoef_c; 
      thchain[indsave,] = th_c; 
      weightchain[indsave,] = weight_c; 
      bzchain[indsave,] = bz_c
      alphachain[indsave] = alpha_c
      phichain[indsave] = phi_c 
      Sigmauchain[indsave,,] = Sigma_c
      likelihoodchain[indsave,] = likelihood_s
    }
    
    if(iscan%%nskip_r==0 && iscan>nskip_r*nburn_r){
      indsave_r = indsave_r +1
      randomeffectchain[indsave_r,,] = randomeffect_c; 
      print(cbind(iscan,indsave,sum(likelihood_s)))
    }
    
  }
  #Compute WAIC
  lppd = sum(log(apply(exp(likelihoodchain),2,mean)))
  WAIC1 = lppd-(2*sum(log(apply(exp(likelihoodchain),2,mean))-apply(likelihoodchain,2,mean)))
  WAIC2 = lppd-sum(apply(likelihoodchain,2,var))#preferred
  LPML = -sum(log(apply(exp(-likelihoodchain),2,mean)))
  
  result<-list(theta = thchain,
               crcoef = crcoefchain,
               FTcoef = FTcoefchain,
               weight = weightchain,
               bz = bzchain,
               alpha = alphachain,
               random_effect = randomeffectchain,
               likelihood = likelihoodchain,
               sigma = Sigmauchain,
               phi = phichain,
               acceptance = acceptance,
               lppd = lppd,
               WAIC1 = WAIC1,
               WAIC2 = WAIC2,
               LPML = LPML);
  return(result);
}



