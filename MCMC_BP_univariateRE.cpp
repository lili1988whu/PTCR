#include <RcppArmadillo.h> 
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

#define ESMALL 1e-10  /* small number */
#define ELARGE 1e+10 /* large number */
#define SYSMIN 1e-30  /* small number */
#define SYSMAX 1e+30 /* large number */
#define LOGSYSMAX 69.07755 /* large number */
#define LOGSYSMIN -69.07755 /* large number */

// [[Rcpp::export]]
double F0(double u){
  double pu;
  pu=0.5*Rf_pnorm5(std::log(u),1.5,0.8,true, false) + 0.5*Rf_pnorm5(std::log(u),2.5,0.3,true, false);
  return(pu);
}
// [[Rcpp::export]]
double densi0(double u){
  double pu;
  pu=0.5*Rf_dnorm4(std::log(u),1.5,0.8,false)/u + 0.5*Rf_dnorm4(std::log(u),2.5,0.3,false)/u;
  return(pu);
}

// [[Rcpp::export]]
double invF0(double u){
  double a = 0.0; double b = 100.0;
  int maxrun = 10000;double tol = 1.0e-10;
  double pa=F0(a);
  double pb=F0(b);
  double newa, newb, ahalf;
  int krun;  
  double sample;
  krun = 1;
  if(u<=pa){sample=a;}if(u>=pb){sample=b;}
  if(u<pb && u>pa){
    newa = a;  newb = b;
    ahalf = (newa+newb)/2;
    while(std::abs(F0(ahalf)-u)>tol && krun<=maxrun){
      if(u<=F0(ahalf)){
        newb = ahalf;
      }
      else{
        newa = ahalf;
      }
      ahalf=(newa+newb)/2;
      krun=krun+1;         
    }
    sample = ahalf;
  }
  return(sample) ;
}


// [[Rcpp::export]]
double sampleT(std::vector<std::string> model, double Fu,arma::rowvec FTx,arma::vec FTcoef,double FTwij){
  
  arma::mat tempmat(1,1);
  double scalef;
  double ratio;
  double FT,sample;
  tempmat = FTx*FTcoef;
  scalef = std::exp(tempmat(0,0)+FTwij);
  
  if(model[0]=="PH"){
    FT = 1.0-std::pow((1.0-Fu),1.0/scalef);
    sample = invF0(FT);
  }
  if(model[0]=="AFT"){
    FT = Fu;
    sample = invF0(FT)/scalef;
  }
  if(model[0]=="PO"){
    ratio = Fu/(1.0-Fu);
    FT = ratio/(ratio+scalef);
    sample = invF0(FT);
  }                      
  return(sample);
}

// [[Rcpp::export]]
arma::vec sampleTs(std::vector<std::string> model,int n,int m,arma::mat crx, arma::vec cru,arma::vec crcoef, 
                   arma::mat FTx,arma::vec FTcoef, arma::vec FTw){
  double u;
  double Fu;
  int ind=1;
  arma::vec t(n*m);
  arma::mat cr(n*m,1);cr.zeros();
  cr = crx*crcoef;
  cr.col(0) = cr.col(0)+cru;
  cr.col(0) = arma::exp(cr.col(0))/(1.0+arma::exp(cr.col(0)));
  arma::mat FT(n*m,1); FT.col(0) = FTw;
    for(int i=1;i<=n;i++){
        for(int j=1;j<=m;j++){
            u = Rf_runif(0,1);
            if (u<= cr(ind-1,0)){t(ind-1) = 10000.00;}
            else {
                Fu = std::log(u)/log(cr(ind-1,0));
                t(ind-1) = sampleT(model, Fu, FTx.row(ind-1),FTcoef,FT(ind-1,0));}
            ind = ind+1;
        }}
  return(t);
}
/* from conditional Ys to cumulative probs*/
// [[Rcpp::export]]
arma::vec Ys_to_weight(arma::vec Ys){
  int nYs = Ys.size();
  arma::vec weight(nYs+1);weight.zeros();
  arma::vec newYs(nYs+1);newYs.ones();
  for(int j=0; j<nYs; j++) newYs[j] = exp(Ys[j]);
  double den = arma::sum(newYs);
  for(int j=0; j<=nYs; j++) weight[j] = newYs[j]/den; 
  return(weight);
}


/////////////////////////////////////////////////////////////////////////
/////////////////// baseline suvival functions///////////////////////////
/////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double S0BP(double t, double th1, double th2, arma::vec w, int BP, int distr){
  if(t<SYSMIN) return(1.0);
  double tmp1 = (std::log(t)+th1)*std::exp(th2);
  int J = w.size();
  double surv=0, logitSt=0, logtemp=0, Ixprev=0, Ft=0;
  if(BP>0){
    if(distr==1){
      Ft = std::exp(tmp1)/(1.0+std::exp(tmp1));
    }else if (distr==2){
      Ft = Rf_pnorm5(tmp1, 0, 1, true, false);
    }else{
      Ft = 1.0-std::exp(-std::exp(tmp1));
    }
    if(Ft<SYSMIN) {Ft=SYSMIN;}
    else{ if(Ft>(1.0-SYSMIN)){Ft = 1.0-SYSMIN;}}
    logitSt = std::log(1.0-Ft)-std::log(Ft); //log(St/(1-St));
    if(logitSt<LOGSYSMIN) return(SYSMIN);
    logtemp = J*std::log(Ft); // J*log(1-St);
    Ixprev = 1.0-std::exp(logtemp); 
    surv = w(0)*Ixprev;
    for(int j=1; j<J; j++){
      logtemp += logitSt + std::log((J-j+1.0)/(j+0.0));
      Ixprev -= std::exp(logtemp);
      surv += w[j]*Ixprev;
    }
  }else{
    if(distr==1){
      surv = 1.0/(1.0+std::exp(tmp1));
    }else if (distr==2){
      surv = Rf_pnorm5(tmp1, 0, 1, false, false);
    }else{
      surv = std::exp(-std::exp(tmp1));
    }
  }

  if(surv>SYSMIN){
    return(surv);
  }else{
    return(SYSMIN);
  }
}
// [[Rcpp::export]]
double F0BP(double t, double th1, double th2, arma::vec w, int BP, int distr){
  double F = 1.0-S0BP(t, th1, th2, w, BP, distr);
  return(F);
}
// [[Rcpp::export]]
double logf0BP(double t, double th1, double th2,arma::vec w, int BP, int distr){
  if(t<SYSMIN) return(LOGSYSMIN);
  double tmp1 = (std::log(t)+th1)*std::exp(th2);
  if(tmp1>LOGSYSMAX) return(LOGSYSMIN);
  int J = w.size();
  double ll=0, logft=0, logitSt=0, logtemp=0, Ft=0;
  if(BP>0){
    if(distr==1){
      Ft = std::exp(tmp1)/(1.0+std::exp(tmp1));
      logft = th2 + tmp1*(1.0-std::exp(-th2)) + th1 - 2.0*std::log(1.0+std::exp(tmp1));
    }else if (distr==2){
      Ft = Rf_pnorm5(tmp1, 0, 1, true, false);
      logft = Rf_dlnorm(t, -th1, std::exp(-th2), true);
    }else{
      Ft = 1.0-std::exp(-std::exp(tmp1));
      logft = th2 + tmp1*(1.0-std::exp(-th2)) + th1 - std::exp(tmp1);
    }
    if(Ft<SYSMIN) return(LOGSYSMIN);
    logitSt = std::log(1.0-Ft)-std::log(Ft); //log(St/(1-St));
    if(logitSt<LOGSYSMIN) return(LOGSYSMIN);
    logtemp = std::log(J) + (J-1.0)*std::log(Ft); // log(J)+(J-1)*log(1-St);
    ll = w(0)*std::exp(logtemp+logft);
    for(int j=1; j<J; j++){
      logtemp += logitSt + std::log((J-j+0.0)/(j+0.0));
      ll += w(j)*std::exp(logtemp+logft);
    }
    ll = std::log(ll);
  }else{
    if(distr==1){
      ll = th2 + tmp1*(1.0-std::exp(-th2)) + th1 - 2.0*std::log(1.0+std::exp(tmp1));
    }else if (distr==2){
      ll = Rf_dlnorm(t, -th1, std::exp(-th2), true);
    }else{
      ll = th2 + tmp1*(1.0-std::exp(-th2)) + th1 - std::exp(tmp1);
    }
  }
  if(ll>LOGSYSMIN) {
    return(ll);
  } else {
    return(LOGSYSMIN);
  }
}

/////////////////////////////////////////////////////////////////////////
//////////////////////////// AFT model //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
// [[Rcpp::export]]
double AFT_BP_logpdf(double t, double th1, double th2, arma::vec w,
                     int BP, int distr, double xibeta){
  double ll = xibeta + logf0BP(exp(xibeta)*t, th1, th2, w, BP, distr);
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log survival function of t given xi
// [[Rcpp::export]]
double AFT_BP_logsurv(double t, double th1, double th2, arma::vec w, 
                      int BP, int distr, double xibeta){
  double temp = S0BP(exp(xibeta)*t, th1, th2, w, BP, distr);
  double ll;
  if(temp<SYSMIN){ll = LOGSYSMIN;}
  else{ll = std::log(temp);}
  return(ll);
}
// log cdf of t given xi
// [[Rcpp::export]]
double AFT_BP_logcdf(double t, double th1, double th2, arma::vec w, 
                     int BP, int distr, double xibeta){
  double temp = 1.0-S0BP(exp(xibeta)*t, th1, th2, w, BP, distr);
  double ll;
  if(temp<SYSMIN){ll = LOGSYSMIN;}
  else{ll = std::log(temp);}
  return(ll);
}
// log (S(t1|xi)-S(t2|xi))
// [[Rcpp::export]]
double AFT_BP_logsurvdiff(double t1, double t2, double th1, double th2, arma::vec w, 
                          int BP, int distr, double xibeta){
  double St1 = S0BP(exp(xibeta)*t1, th1, th2, w, BP, distr);
  double St2 = S0BP(exp(xibeta)*t2, th1, th2, w, BP, distr);
  double temp = std::abs(St1 - St2);
  double ll;
  if(temp<SYSMIN){ll = LOGSYSMIN;}
  else{ll = std::log(temp);}
  return(ll);
}

// Calculate loglikelihood for each obervation i
// [[Rcpp::export]]
arma::vec AFT_BP_logliki(arma::vec t1, arma::vec t2, 
                         Rcpp::IntegerVector type, double th1, double th2, arma::vec w,
                         int BP, int distr, arma::vec FXbeta,arma::vec Fv,
                         arma::vec crXbeta,arma::vec crv,arma::mat missing){
  arma::vec res(type.size());res.zeros();
  double temp,temp2;
  double logl;
  for(int i=0; i<type.size(); i++){
    temp2 = exp(crXbeta(i)+crv(i))/(1.0+exp(crXbeta(i)+crv(i)));
    if(temp2<SYSMIN){temp = -LOGSYSMIN;}
    else{temp =-log(temp2);}
    //temp = exp(crXbeta(i)+crv(i));
    if(type[i]==1){
      res(i) = AFT_BP_logpdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i))+log(temp)
      -exp(AFT_BP_logcdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp; /*Observed*/
    }else if(type[i]==2){
      logl = 1.0-exp(-exp(AFT_BP_logcdf(t2(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp);
      if (logl<SYSMIN){res(i) = LOGSYSMIN;}
      else{ res(i) = log(logl);}/*Left censored*/
    }else if(type[i]==3){
      res(i) = -exp(AFT_BP_logcdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp;/*Right censored*/
    }else if (type[i]==4){
      logl = exp(-exp(AFT_BP_logcdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp)
      -exp(-exp(AFT_BP_logcdf(t2(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp);
      if (logl<SYSMIN){res(i) = LOGSYSMIN;}
      else{res(i) = log(logl);}
      /*Interval-observed*/
/*Interval-observed*/
    }else{
      if(missing(i,0)>0.5){res(i) = AFT_BP_logpdf(missing(i,1), th1, th2, w, BP, distr, FXbeta(i)+Fv(i))+log(temp)
        -exp(AFT_BP_logcdf(missing(i,1), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp; /*Observed*/
      }
      else{res(i) = -exp(AFT_BP_logcdf(missing(i,1), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp;}
    }
  } 
  return(res);
}

// [[Rcpp::export]]
arma::vec AFT_BP_logliki_para(arma::vec t1, arma::vec t2, 
                         Rcpp::IntegerVector type, double th1, double th2, arma::vec w,
                         int BP, int distr, arma::vec FXbeta,arma::vec Fv,
                         arma::vec crXbeta,arma::vec crv){
  arma::vec res(type.size());res.zeros();
  double temp,temp2=0.0;
  for(int i=0; i<type.size(); i++){
    temp2 = exp(crXbeta(i)+crv(i))/(1.0+exp(crXbeta(i)+crv(i)));
    if(temp2<SYSMIN){temp = -LOGSYSMIN;}
    else{temp =-log(temp2);}
    //temp = exp(crXbeta(i)+crv(i));
    if(type[i]==1){
      res(i) = AFT_BP_logpdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i))+log(temp)
      -exp(AFT_BP_logcdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp; /*Observed*/
    }else if(type[i]==2){
      res(i) = log(1.0-exp(-exp(AFT_BP_logcdf(t2(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp ) ) ;/*Left censored*/
    }else if(type[i]==3){
      res(i) = -exp(AFT_BP_logcdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp;/*Right censored*/
    }else if (type[i]==4){
      res(i) = log(exp(-exp(AFT_BP_logcdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp)
                     -exp(-exp(AFT_BP_logcdf(t2(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp));/*Interval-observed*/
    }
  } 
  return(res);
}


/////////////////////////////////////////////////////////////////////////
//////////////////////////// PH model //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
// [[Rcpp::export]]
double PH_BP_logpdf(double t, double th1, double th2, arma::vec w,
                    int BP, int distr, double xibeta){
  double ll = xibeta + logf0BP(t, th1, th2, w, BP, distr);
  ll += (exp(xibeta)-1.0)*log(S0BP(t, th1, th2, w, BP, distr));
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log survival function of t given xi
// [[Rcpp::export]]
double PH_BP_logsurv(double t, double th1, double th2, arma::vec w, 
                     int BP, int distr, double xibeta){
  double ll = exp(xibeta)*log(S0BP(t, th1, th2, w, BP, distr));
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log cdf of t given xi
// [[Rcpp::export]]
double PH_BP_logcdf(double t, double th1, double th2, arma::vec w, 
                    int BP, int distr, double xibeta){
  double temp = 1.0-exp( exp(xibeta)*log(S0BP(t, th1, th2, w, BP, distr)) );
  double ll;
  if(temp<SYSMIN){ll = LOGSYSMIN;}
  else{ll = std::log(temp);}
  return(ll);
}
// log (S(t1|xi)-S(t2|xi))
// [[Rcpp::export]]
double PH_BP_logsurvdiff(double t1, double t2, double th1, double th2, arma::vec w, 
                         int BP, int distr, double xibeta){
  double St1 = exp( exp(xibeta)*log(S0BP(t1, th1, th2, w, BP, distr)) );
  double St2 = exp( exp(xibeta)*log(S0BP(t2, th1, th2, w, BP, distr)) );
  double temp = std::abs(St1 - St2);
  double ll;
  if(temp<SYSMIN){ll = LOGSYSMIN;}
  else{ll = std::log(temp);}
  return(ll);
}

// Calculate loglikelihood for each obervation i
// [[Rcpp::export]]
arma::vec PH_BP_logliki(arma::vec t1, arma::vec t2, 
                         Rcpp::IntegerVector type, double th1, double th2, arma::vec w,
                         int BP, int distr, arma::vec FXbeta,arma::vec Fv,
                         arma::vec crXbeta,arma::vec crv,arma::mat missing){
  arma::vec res(type.size());res.zeros();
  double temp,temp2;
  double logl;
  for(int i=0; i<type.size(); i++){
    /*temp = exp(crXbeta(i)+crv(i));*/
    temp2 = exp(crXbeta(i)+crv(i))/(1.0+exp(crXbeta(i)+crv(i)));
    if(temp2<SYSMIN){temp = -LOGSYSMIN;}
    else{temp =-log(temp2);}
    if(type[i]==1){
      res(i) = PH_BP_logpdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i))+log(temp)
      -exp(PH_BP_logcdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp; /*Observed*/
    }else if(type[i]==2){
      logl = 1.0-exp(-exp(PH_BP_logcdf(t2(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp);
      if (logl<SYSMIN){res(i) = LOGSYSMIN;}
      else{ res(i) = log(logl);}/*Left censored*/
    }else if(type[i]==3){
      res(i) = -exp(PH_BP_logcdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp;/*Right censored*/
    }else if (type[i]==4){
      logl = exp(-exp(PH_BP_logcdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp)
      -exp(-exp(PH_BP_logcdf(t2(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp);
      if (logl<SYSMIN){res(i) = LOGSYSMIN;}
      else{res(i) = log(logl);}
      /*Interval-observed*/
    }else{
      if(missing(i,0)>0.5){res(i) = PH_BP_logpdf(missing(i,1), th1, th2, w, BP, distr, FXbeta(i)+Fv(i))+log(temp)
        -exp(PH_BP_logcdf(missing(i,1), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp; /*Observed*/
      }
      else{res(i) = -exp(PH_BP_logcdf(missing(i,1), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp;}
    }
  } 
  return(res);
}
// [[Rcpp::export]]
arma::vec PH_BP_logliki_para(arma::vec t1, arma::vec t2, 
                              Rcpp::IntegerVector type, double th1, double th2, arma::vec w,
                              int BP, int distr, arma::vec FXbeta,arma::vec Fv,
                              arma::vec crXbeta,arma::vec crv){
  arma::vec res(type.size());res.zeros();
  double temp,temp2;
  for(int i=0; i<type.size(); i++){
    temp2 = exp(crXbeta(i)+crv(i))/(1.0+exp(crXbeta(i)+crv(i)));
    if(temp2<SYSMIN){temp = -LOGSYSMIN;}
    else{temp =-log(temp2);}
    if(type[i]==1){
      res(i) = PH_BP_logpdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i))+log(temp)
      -exp(PH_BP_logcdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp; /*Observed*/
    }else if(type[i]==2){
      res(i) = log(1.0-exp(-exp(PH_BP_logcdf(t2(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp ) ) ;/*Left censored*/
    }else if(type[i]==3){
      res(i) = -exp(PH_BP_logcdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp;/*Right censored*/
    }else if (type[i]==4){
      res(i) = log(exp(-exp(PH_BP_logcdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp)
                     -exp(-exp(PH_BP_logcdf(t2(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp));/*Interval-observed*/
    }
  } 
  return(res);
}
/////////////////////////////////////////////////////////////////////////
//////////////////////////// PO model //////////////////////////////////
/////////////////////////////////////////////////////////////////////////
// log densitt of t given xi
// [[Rcpp::export]]
double PO_BP_logpdf(double t, double th1, double th2, arma::vec w,
                    int BP, int distr, double xibeta){
  double ll = logf0BP(t, th1, th2, w, BP, distr)-xibeta;
  double temp = 1.0+(exp(-xibeta)-1.0)*S0BP(t, th1, th2, w, BP, distr);
  if(temp<SYSMIN){ll +=-2.0*LOGSYSMIN;}
  else{ll +=-2.0*std::log(temp);}
  if(ll<LOGSYSMIN) ll=LOGSYSMIN;
  return(ll);
}
// log survival function of t given xi
// [[Rcpp::export]]
double PO_BP_logsurv(double t, double th1, double th2, arma::vec w, 
                     int BP, int distr, double xibeta){
  double S0t = S0BP(t, th1, th2, w, BP, distr);
  double temp = 1.0+(exp(-xibeta)-1.0)*S0t;
  double ll;
  if(temp <SYSMIN){ll = log(S0t)-xibeta-LOGSYSMIN; }
  else{ll = log(S0t)-xibeta-std::log(temp); }
  return(ll);
}
// log cdf of t given xi
// [[Rcpp::export]]
double PO_BP_logcdf(double t, double th1, double th2, arma::vec w, 
                    int BP, int distr, double xibeta){
  double S0t = S0BP(t, th1, th2, w, BP, distr);
  double temp1, temp2;
  double y1,y2;
  temp1 = 1.0-S0t; temp2 = 1.0+(exp(-xibeta)-1.0)*S0t;
  if(temp1<SYSMIN){y1=SYSMIN;}
  else{y1 = std::log(temp1);}
  if(temp2<SYSMIN){y2=SYSMIN;}
  else{y2 = std::log(temp2);}
  double ll = y1-y2;
  return(ll);
}
// log (S(t1|xi)-S(t2|xi))
// [[Rcpp::export]]
double PO_BP_logsurvdiff(double t1, double t2, double th1, double th2, arma::vec w, 
                         int BP, int distr, double xibeta){
  double S0t1 = S0BP(t1, th1, th2, w, BP, distr);
  double S0t2 = S0BP(t2, th1, th2, w, BP, distr);
  double exibeta = exp(-xibeta);
  double temp = std::abs(exibeta*S0t1/(1.0+(exibeta-1.0)*S0t1) - exibeta*S0t2/(1.0+(exibeta-1.0)*S0t2));
  double ll;
  if(temp<SYSMIN){ll = LOGSYSMIN;}
  else{ll = std::log(temp);}
  return(ll);
}


// Calculate loglikelihood for each obervation i
// [[Rcpp::export]]
arma::vec PO_BP_logliki(arma::vec t1, arma::vec t2, 
                        Rcpp::IntegerVector type, double th1, double th2, arma::vec w,
                        int BP, int distr, arma::vec FXbeta,arma::vec Fv,
                        arma::vec crXbeta,arma::vec crv,arma::mat missing){
  arma::vec res(type.size());res.zeros();
  double temp,temp2;
  double logl;
  for(int i=0; i<type.size(); i++){
    temp2 = exp(crXbeta(i)+crv(i))/(1.0+exp(crXbeta(i)+crv(i)));
    if(temp2<SYSMIN){temp = -LOGSYSMIN;}
    else{temp =-log(temp2);}
    if(type[i]==1){
      res(i) = PO_BP_logpdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i))+log(temp)
      -exp(PO_BP_logcdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp; /*Observed*/
    }else if(type[i]==2){
      logl = 1.0-exp(-exp(PO_BP_logcdf(t2(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp);
      if (logl<SYSMIN){res(i) = LOGSYSMIN;}
      else{ res(i) = log(logl);}/*Left censored*/
    }else if(type[i]==3){
      res(i) = -exp(PO_BP_logcdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp;/*Right censored*/
    }else if (type[i]==4){
      logl = exp(-exp(PO_BP_logcdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp)
      -exp(-exp(PO_BP_logcdf(t2(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp);
      if (logl<SYSMIN){res(i) = LOGSYSMIN;}
      else{res(i) = log(logl);}
      /*Interval-observed*/
      /*Interval-observed*/
    }else{
      if(missing(i,0)>0.5){res(i) = PO_BP_logpdf(missing(i,1), th1, th2, w, BP, distr, FXbeta(i)+Fv(i))+log(temp)
        -exp(PO_BP_logcdf(missing(i,1), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp; /*Observed*/
      }
      else{res(i) = -exp(PO_BP_logcdf(missing(i,1), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp;}
    }
  } 
  return(res);
}

// [[Rcpp::export]]
arma::vec PO_BP_logliki_para(arma::vec t1, arma::vec t2, 
                             Rcpp::IntegerVector type, double th1, double th2, arma::vec w,
                             int BP, int distr, arma::vec FXbeta,arma::vec Fv,
                             arma::vec crXbeta,arma::vec crv){
  arma::vec res(type.size());res.zeros();
  double temp,temp2;
  for(int i=0; i<type.size(); i++){
    temp2 = exp(crXbeta(i)+crv(i))/(1.0+exp(crXbeta(i)+crv(i)));
    if(temp2<SYSMIN){temp = -LOGSYSMIN;}
    else{temp =-log(temp2);}
    if(type[i]==1){
      res(i) = PO_BP_logpdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i))+log(temp)
      -exp(PO_BP_logcdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp; /*Observed*/
    }else if(type[i]==2){
      res(i) = log(1.0-exp(-exp(PO_BP_logcdf(t2(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp ) ) ;/*Left censored*/
    }else if(type[i]==3){
      res(i) = -exp(PO_BP_logcdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp;/*Right censored*/
    }else if (type[i]==4){
      res(i) = log(exp(-exp(PO_BP_logcdf(t1(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp)
                     -exp(-exp(PO_BP_logcdf(t2(i), th1, th2, w, BP, distr, FXbeta(i)+Fv(i)))*temp));/*Interval-observed*/
    }
  } 
  return(res);
}

// [[Rcpp::export]]
arma::mat matrixproduct(arma::mat Randomz, arma::vec randomeffect, int n, int m, int sizere){
  arma::mat outcome (n*m,1);outcome.zeros();
  arma::mat temp(n*sizere,1);temp.zeros();
  temp.col(0) = randomeffect;
  int ibegin;
  int iend;
  if(sizere == m){outcome = randomeffect;}
  else{
  for(int i=1;i<=n;i++){
    ibegin = (i-1)*m;
    iend = i*m-1;
    outcome.submat(ibegin,0,iend,0) = Randomz.submat(ibegin,0,iend,sizere-1)*temp.submat((i-1)*sizere,0,i*sizere-1,0);
  }}
  return (outcome);
  }



// [[Rcpp::export]]
arma::vec likelihoodv(std::vector<std::string> model,int BP,int m, arma::vec t1, arma::vec t2, 
                  Rcpp::IntegerVector type, double th1, double th2, arma::vec w,
                   int distr, arma::vec FXbeta,arma::vec Fv,
                  arma::vec crXbeta,arma::vec crv,arma::mat missing){
  arma::vec res(type.size());
  if(model[0] =="PH"){
    res = PH_BP_logliki(t1, t2, type, th1, th2, w, BP, distr, FXbeta,Fv,crXbeta,crv,missing);
  }
  if(model[0] =="AFT"){
    res = AFT_BP_logliki(t1, t2, type, th1, th2, w, BP, distr, FXbeta,Fv,crXbeta,crv,missing);
  }
  if(model[0] =="PO"){
    res = PO_BP_logliki(t1, t2, type, th1, th2, w, BP, distr, FXbeta,Fv,crXbeta,crv,missing);
  }
  int n = type.size()/m;
  arma::vec loglikeliv(n);loglikeliv.zeros();
  int ibegin, iend;
  for(int i=1;i<=n;i++){
    ibegin = (i-1)*m;
    iend = i*m-1;
    loglikeliv(i-1) = arma::sum(res.subvec(ibegin,iend));
  }
  return(loglikeliv);
}

// [[Rcpp::export]]
arma::vec likelihoodv_para(std::vector<std::string> model,int BP,int m,arma::vec t1, arma::vec t2, 
                      Rcpp::IntegerVector type, double th1, double th2, arma::vec w, int distr, arma::vec FXbeta,arma::vec Fv,
                      arma::vec crXbeta,arma::vec crv){
  arma::vec res(type.size());
  if(model[0] =="PH"){
    res = PH_BP_logliki_para(t1, t2, type, th1, th2, w, BP, distr, FXbeta,Fv,crXbeta,crv);
  }
  if(model[0] =="AFT"){
    res = AFT_BP_logliki_para(t1, t2, type, th1, th2, w, BP, distr, FXbeta,Fv,crXbeta,crv);
  }
  if(model[0] =="PO"){
    res = PO_BP_logliki_para(t1, t2, type, th1, th2, w, BP, distr, FXbeta,Fv,crXbeta,crv);
  }
  int n = type.size()/m;
  arma::vec loglikeliv(n);loglikeliv.zeros();
  int ibegin, iend;
  for(int i=1;i<=n;i++){
    ibegin = (i-1)*m;
    iend = i*m-1;
    loglikeliv(i-1) = arma::sum(res.subvec(ibegin,iend));
  }
  return(loglikeliv);
}

// [[Rcpp::export]]
double invF(double u, double th1, double th2, arma::vec w, int BP, int distr){
  arma::vec ends(2),pends(2);ends.zeros();pends.zeros();
  arma::vec newends(2),pnewends(2);newends.zeros();pnewends.zeros();
  ends(0) = 0.0; ends(1) = 500.0;
  int maxrun = 10000;double tol = 1.0e-10;
  pends(1) = F0BP(ends(1),th1, th2, w,BP,distr);
  pends(0) = F0BP(ends(0),th1, th2, w,BP,distr);
  double newa, newb, ahalf,sample;
  int krun =1;  
  if(u<=pends(0)){sample=ends(0);}
  if(u>=pends(1)){sample=ends(1);}
  if(u<pends(1) && u>pends(0)){
    newa = ends(0);  newb = ends(1);
    ahalf = (newa+newb)/2;
    newends(0)=newa;newends(1)=ahalf;
    pnewends(1) = F0BP(newends(1),th1, th2, w,BP,distr);
    pnewends(0) = F0BP(newends(0),th1, th2, w,BP,distr);
    while(std::abs(pnewends(1)-u)>tol && krun<=maxrun){
      if(u<=pnewends(1)){
        newb = ahalf;
      }
      else{
        newa = ahalf;
      }
      ahalf=(newa+newb)/2;
      krun=krun+1; 
      newends(0)=newa;newends(1)=ahalf;
      pnewends(1) = F0BP(newends(1),th1, th2, w,BP,distr);
      pnewends(0) = F0BP(newends(0),th1, th2, w,BP,distr);
    }
    sample = ahalf;
  }
  return(sample) ;
}



// [[Rcpp::export]]
arma::mat missing_impute_t(std::vector<std::string> model,int BP, int distr, double maxc,arma::vec t1, arma::vec t2, 
                           Rcpp::IntegerVector type, double th1, double th2, arma::vec w,
                           arma::vec FXbeta,arma::vec Fv,
                           arma::vec crXbeta,arma::vec crv){
    int n = type.size();
    arma::mat nt(n,2); nt.zeros(); 
    double z, u,nz,cr,FX;
    
    for(int j=1;j<=n;j++){
      nt(j-1,1) = 0.0;
      if(type(j-1)==5){
       cr = -log(exp(crXbeta(j-1)+crv(j-1))/(1.0+exp(crXbeta(j-1)+crv(j-1))));
       // cr = exp(crXbeta(j-1)+crv(j-1));
       FX = exp(FXbeta(j-1)+Fv(j-1));
       u = Rf_runif(0.0, 1.0);
      if(u<=std::exp(-cr)){nt(j-1,1) = 1000.0;} 
      else{
        /*obtain the baseline probability from different models*/
          z = -std::log(u)/cr; 
        
        if(model[0]=="PH"){
          nz = 1.0-std::pow((1.0-z),1.0/FX);
        }
        if(model[0]=="PO"){
          double ratio = z/(1.0-z);
          nz = ratio/(ratio+FX);
        }
        
        if(model[0]=="AFT"){
          nz = z;
        }
        
       /*Sample from the baseline*/
        nt(j-1,1) = invF(nz, th1,th2,w, BP, distr);
       
          if(model[0]=="AFT"){
            nt(j-1,1)=nt(j-1,1)/FX;
          }
         
        } 
      if(nt(j-1,1)>maxc){nt(j-1,1)=maxc;nt(j-1,0) = 0.0;}
      else{nt(j-1,0)=1.0;}
      
      }
      
      }
  return(nt);  
}



/*Priors*/




// [[Rcpp::export]]
double log_dnorm (arma::vec x, arma::vec meanv, arma::mat precisionm,int d) {
  arma::rowvec diff1 (d);
  arma::colvec diff2 (d);
  diff1.zeros();
  diff2.zeros();
  diff2 = x-meanv;
  diff1 = diff2.t();
  double tempp=-arma::as_scalar(diff1*precisionm*diff2)/(double)2;
  return(tempp);
}

// [[Rcpp::export]]
double pwu(int n,int m, arma::mat FTwcru, arma::mat precisionm) {
  arma::rowvec diff1 (2*m);
  arma::colvec diff2 (2*m);
  diff1.zeros();
  diff2.zeros();
  double tempp=0.0;
  for(int i=1;i<=n;i++){
    diff1 = FTwcru.row(i-1);
    diff2 = diff1.t(); 
    tempp=tempp-arma::as_scalar(diff1*precisionm*diff2)/(double)2;
  }
  return(tempp);
}

// [[Rcpp::export]]
arma::vec rmnorm(arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(1, ncols);
  arma::rowvec sample (ncols);
  sample = arma::repmat(mu, 1, 1).t() + Y * arma::chol(sigma);
  return (sample.t());
}



// [[Rcpp::export]] 
arma::mat rwishart(int nu,  arma::mat V){
  // function to draw from Wishart (nu,V) and IW
  // 
  // W ~ W(nu,V) E[W]=nuV
  // 
  // WI=W^-1 E[WI]=V^-1/(nu-m-1)
  
  RNGScope rngscope;
  
  int m = V.n_rows;
  
  arma::mat Vm(V.begin(), m, m, false);
  arma::mat T(m,m);
  arma::mat IW(m,m);
  for(int i = 0; i < m; i++) {
    T(i,i) = sqrt(R::rchisq(nu-i));
  }
  
  for(int j = 0; j < m; j++) {  
    for(int i = j+1; i < m; i++) {    
      T(i,j) = R::norm_rand();
    }}
  arma::mat C = trans(trimatl(T)) * chol(Vm);
  arma::mat CI = C.i();
  
  // C is the upper triangular root of Wishart 
  // therefore, W=C'C  this is the LU decomposition 
  // Inv(W) = CICI'    this is the UL decomp *not LU*!
  
  // W is Wishart draw, IW is W^-1
  IW = CI * trans(CI);
  return(IW);
} 
// [[Rcpp::export]] 
double samplesigma(int n, int a, arma::vec us, double Lambda){
  double df = double(n+a)/2.0;
  double tempmat; 
  double samplesig;
  tempmat=(double)(a)/Lambda; //a is the degree of freedom
  for(int i=1;i<=n;i++){
    tempmat=tempmat+us(i-1)*us(i-1)/2.0;
  }
  samplesig = tempmat /Rf_rgamma(df,1.0);
  return(samplesig);
}

// [[Rcpp::export]] 
double sampleLambda(int a, double zetasq, double Sigma){
  double sampledLambda; 
  double shapeg = (double)(a+1.0)/(double)(2);
  double rateg =0.0;
  rateg = 1.0/zetasq+(double)(a)/Sigma;
  sampledLambda=rateg/Rf_rgamma(shapeg,1.0);
  
  return(sampledLambda);
}

// [[Rcpp::export]] 
arma::mat Sigmastarf(arma::mat Sigmau, arma::mat Sigmae, double phi,int m){
  arma::mat temp (2*m,2*m); temp.zeros();
  
  temp.submat(0,0,m-1,m-1) = phi*phi*Sigmau+Sigmae;
  temp.submat(m,0,2*m-1,m-1) = phi*Sigmau;
  temp.submat(0,m,m-1,2*m-1) = phi*Sigmau;
  temp.submat(m,m,2*m-1,2*m-1) = Sigmau;
 return(temp);
}



// [[Rcpp::export]] 
double recursivemean(double pre_mean, double new_x, int samplesize){
  double new_mean;
  new_mean = pre_mean*(double(samplesize)-1.0)/(double (samplesize))+new_x/(double(samplesize));
  return(new_mean);
}


// [[Rcpp::export]] 
double recursivecov(int samplesize,double pre_cov,double pre_mean,double new_x,double epsilon,double s_d){
  double n_eff = (double) (samplesize);
  double new_cov;
  double new_mean;
  new_mean = recursivemean(pre_mean, new_x, n_eff);
  new_cov = (n_eff-2.0)/(n_eff-1.0)*pre_cov+s_d*(pre_mean*pre_mean-(n_eff)/(n_eff-1.0)*new_mean*new_mean
                                                   +new_x*new_x/(n_eff-1.0)+epsilon/(n_eff-1.0));
  
  return(new_cov);
}
// [[Rcpp::export]] 
arma::vec recursivemean_vector(arma::vec pre_mean, arma::vec new_x, int samplesize){
  int n = pre_mean.size ();
  arma::vec new_mean(n);new_mean.zeros();
  new_mean = pre_mean*(double(samplesize)-1.0)/(double (samplesize))+new_x/(double(samplesize));
  return(new_mean);
}


// [[Rcpp::export]] 
arma::mat recursivecov_vector(int d, int samplesize,arma::mat pre_cov,arma::vec pre_mean,arma::vec new_x,double epsilon,double s_d){
  double n_eff = (double) (samplesize);
  arma::mat new_cov (d,d); new_cov.zeros();
  arma::vec new_mean (d);new_mean.zeros();
  arma::mat temp_o (d,1);temp_o.zeros();
  arma::mat temp_n (d,1);temp_n.zeros();
  arma::mat temp_x (d,1);temp_x.zeros();
  arma::mat covm(d,d);covm.zeros();
  covm.eye();covm = epsilon*covm;
  new_mean = recursivemean_vector(pre_mean, new_x, n_eff);

    temp_o.col(0) = pre_mean;
    temp_n.col(0) = new_mean;
    temp_x.col(0) = new_x;
    new_cov = (n_eff-2.0)/(n_eff-1.0)*pre_cov+s_d*(temp_o*trans(temp_o)-(n_eff)/(n_eff-1.0)*temp_n*trans(temp_n)
                                                                                 +temp_x*trans(temp_x)/(n_eff-1.0)+covm/(n_eff-1.0));

  return(new_cov);
}



// [[Rcpp::export]] 
arma::mat samplerandomeffect(arma::vec FTw, arma::vec cru, arma::vec cov_randomeffect, double phi, int n)
{
  double wuiold,wuinew,covm;
  arma::mat newFTwcru(n,2);newFTwcru.zeros();
  
  for(int l=1;l<=n;l++){
    wuiold = cru(l-1);
    covm = cov_randomeffect(l-1);
    wuinew = Rf_rnorm(wuiold, sqrt(covm));
    newFTwcru(l-1,0) = phi*wuinew;
    newFTwcru(l-1,1) = wuinew;
  }
  
  return(newFTwcru);
}
// [[Rcpp::export]] 
arma::mat randomeffect (arma::vec cru_o,arma::vec cru_n,arma::vec likelihoodv_o, arma::vec likelihoodv_n, double invsigmastar,int n,int m)
{
  int ibegin1;
  int iend1;
  double wuiold;
  double wuinew;
  arma::mat randomu(2,n);randomu.zeros();
  double likelihoodi_o,likelihoodi_n;
  double priori_o,priori_n;
  for(int l=1;l<=n;l++){
    wuiold = cru_o(l-1);
    wuinew = cru_n(l-1);
    likelihoodi_o = likelihoodv_o(l-1);
    likelihoodi_n = likelihoodv_n(l-1);
    priori_o = -wuiold*invsigmastar*wuiold/2.0;
    priori_n = -wuinew*invsigmastar*wuinew/2.0;
    if(log(Rf_runif(0.0, 1.0))<(likelihoodi_n+priori_n-likelihoodi_o-priori_o)){
      randomu(0,l-1) = wuinew;randomu(1,l-1)=randomu(1,l-1)+1.0;
    }
    else{
      randomu(0,l-1) = wuiold;
    }
    /*cout <<"u"<< likelihoodi_n+priori_n <<"z"<< likelihoodi_o-priori_o<< " " <<Rf_runif(0.0, 1.0)<<".\n";*/
  }
return(randomu);
}
  
