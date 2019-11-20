load("~/Dropbox/research/Active/Dipankar/compiled5.RData")
data_I = cbind(seq(1,1586,1),X,L,U)
data= na.omit(data_I)
sub = seq(1,1584,1)#1584
L.d = data[sub,5:32]
U.d = data[sub,33:60]
n = length(data[sub,1])
m = 28

#Baseline CAL averages
ind=1
cal=rep(0,n*m)
for(i in 1:n){
  for(j in 1:m){
    cal[ind] = mean(na.omit(Y[[i]][1,]))
    ind = ind+1
  }
}
cal = (cal-mean(cal))/max(cal)
bs<-function(x,bk)
{a31=bk[3]-bk[1]; a32=bk[3]-bk[2];a42=bk[4]-bk[2];a43=bk[4]-bk[3];a41=bk[4]-bk[1];a21=bk[2]-bk[1]
t=x
if(x<=bk[1]) bsv=0
if(x<=bk[2] && x>bk[1]) bsv=(t-bk[1])**2/(a31*a21)
if(x<=bk[3] && x>bk[2]) bsv=(t-bk[1])*(bk[3]-t)/(a31*a32)+(bk[4]-t)*(t-bk[2])/(a42*a32)
if(x<=bk[4] && x>bk[3]) bsv=(t-bk[4])**2/(a42*a43)
if(x>bk[4]) bsv=0 
return(bsv)
}
N_Bspline = 6
Jbeta=N_Bspline
calcov = matrix(0,nrow =n*m,ncol=Jbeta)

min = min(cal)
max = max(cal)
N_knots = min +(max-min)/(N_Bspline-2)*seq(-2,N_Bspline,1)
Knotsmatrix = matrix(0,ncol=4,nrow=N_Bspline)
for(i in 1:N_Bspline){
  Knotsmatrix[i,]=N_knots[i:(i+3)]
}
for(i in 1:(n*m)){
  for(j in 1:(N_Bspline)){calcov[i,j] = bs(cal[i],Knotsmatrix[j,])}
}
#X has columns for gender (M/F), tobacco (Y/N) and diabetes (Y/N). 
data[,2 ] = data[,2]-1
type_d = rep(0,n*m);t1 = rep(0,n*m);t2 = rep(0,n*m)
tempx = matrix(0,nrow = n*m,ncol = 3)

ind=1
for(i in 1:n){
  for(s in 1:m){
    if(L.d[i,s]==U.d[i,s]){ t1[ind] = L.d[i,s];type_d[ind] = 1;tempx[ind,1:3] = data[i,c(2,3,4)];ind=ind+1 }#observed
    else{
      if(L.d[i,s]==0 && U.d[i,s]!=Inf){ t2[ind] = U.d[i,s];type_d[ind] = 2;tempx[ind,1:3] = data[i,c(2,3,4)];ind=ind+1  }#left censoring
      else{
        if(L.d[i,s]>0 && U.d[i,s]==Inf){ t1[ind] = L.d[i,s];type_d[ind] = 3;tempx[ind,1:3] = data[i,c(2,3,4)];ind=ind+1 }#right censoring
        else{
          if(L.d[i,s]==0 && U.d[i,s]==Inf){ type_d[ind] = 5;tempx[ind,1:3] = data[i,c(2,3,4)];ind=ind+1 }#missing
        }}
    }
  }
}
incisor<-c(6,7,8,9,20,21,22,23) #baseline
canine<-c(5,10,19,24)
premolar<-c(3,4,11,12,17,18,25,26)
molar<-c(1,2,13,14,15,16,27,28)

teethind = NULL
for(i in 1:m){
  teethind[i] = 1*is.element(i, premolar)+1*is.element(i, molar)
}
teethx=NULL
for(i in 1:n){
  seq = ((i-1)*m+1):(i*m)
  teethx[seq] = teethind
}

#setup1


FTx1 = cbind(tempx,teethx,cal)
crx1 = rep(1,n*m)




##########################################################################################################
#nonparametric fits
toothneighbor<-matrix(0,ncol=28,nrow=28)
toothneighbor[1,2] = 1
toothneighbor[2,1] = 1; toothneighbor[2,3] = 1
toothneighbor[3,2] = 1; toothneighbor[3,4] = 1
toothneighbor[4,3] = 1; toothneighbor[4,5] = 1
toothneighbor[5,4] = 1; toothneighbor[5,6] = 1
toothneighbor[6,5] = 1; toothneighbor[6,7] = 1
toothneighbor[7,6] = 1; toothneighbor[7,8] = 1
toothneighbor[8,7] = 1; toothneighbor[8,9] = 1
toothneighbor[9,8] = 1; toothneighbor[9,10] = 1
toothneighbor[10,9] = 1; toothneighbor[10,11] = 1
toothneighbor[11,10] = 1; toothneighbor[11,12] = 1
toothneighbor[12,11] = 1; toothneighbor[12,13] = 1
toothneighbor[13,12] = 1; toothneighbor[13,14] = 1
toothneighbor[14,13] = 1
toothneighbor[15,16] = 1
toothneighbor[16,15] = 1; toothneighbor[16,17] = 1
toothneighbor[17,16] = 1; toothneighbor[17,18] = 1
toothneighbor[18,17] = 1; toothneighbor[18,19] = 1
toothneighbor[19,18] = 1; toothneighbor[19,20] = 1
toothneighbor[20,19] = 1; toothneighbor[20,21] = 1
toothneighbor[21,20] = 1; toothneighbor[21,22] = 1
toothneighbor[22,21] = 1; toothneighbor[22,23] = 1
toothneighbor[23,22] = 1; toothneighbor[23,24] = 1
toothneighbor[24,23] = 1; toothneighbor[24,25] = 1
toothneighbor[25,24] = 1; toothneighbor[25,26] = 1
toothneighbor[26,25] = 1; toothneighbor[26,27] = 1
toothneighbor[27,26] = 1; toothneighbor[27,28] = 1
toothneighbor[28,27] = 1

Randomx1 = matrix(0,ncol=m,nrow=n*m)
Randomx2 = matrix(0,ncol=2,nrow=n*m)

for(i in 1:n)
{ ibegin = (i-1)*m+1; iend = i*m
Randomx1[ibegin:iend,] = diag(m)
for(j in 1:m){
  Randomx2[ibegin+j-1,1] = 1*is.element(j, incisor)+1*is.element(j, canine)
  Randomx2[ibegin+j-1,2] = 1*is.element(j, premolar)+1*is.element(j, molar)
}
}
#mcmc
library(BDgraph)
adj.g = toothneighbor
adj.g.wish = upper.tri(toothneighbor, diag = FALSE)
Dn = diag(apply(adj.g,1,sum))



###########################################

