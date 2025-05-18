# mvmr_local.R

# Usage: [VHat,betaHat,ci,betaSD] = mvmr_local(Y,D,Z,X,intercept,alpha,a0,Cbeta)
#input:
# Y: nx1 outcome vector
# D: nxq treatment matrix
# Z: nxpz candidate IVs
# X: nxpx corvariates
# intercept: whether or not introduce a intercept in linear regression
# alpha: confidence level
# a0: tuning parameter to adjust threshoulding level in first stage;
# Cbeta: bound for |beta_j|

#output:
# VHat: estimated exogenous IVs
# betaHat: estimated treatments effects
# ci: alpha-level confidence intervals for treatments effects
# betaSD: standard deviation of estimated treatments effects









Local.summary <- function(GamHat,gamHat,Sig_w,Omega_n,n,a0=1,Cbeta=1,alpha=0.05,sep=NULL) {
  # Include intercept
  pz=length(GamHat)
  q=length(gamHat)/pz
  p=nrow(SigHat)
  
  
  #=========== estimate V =================================================
  if(is.null(sep)){
    b=seq(-Cbeta,Cbeta,q*Cbeta/sqrt(n))# candidate values
  }
  else b=seq(-Cbeta,Cbeta,sep)
  
  candidates=rep(list(b), q)
  B=expand.grid(candidates)
  VFlag=matrix(0,nrow=pz,ncol=length(b)^q) # Voting matrix
  
  for(m in 1:length(b)^q){
    beta_m=as.numeric(B[m,])
    pi_m=GamHat-gamHat%*%beta_m
    mat_beta=kronecker(t(c(1,-beta_m)),diag(pz))
    thresh_m=sqrt(a0*log(n))*sqrt(diag(mat_beta%*%Omega_n%*%t(mat_beta))) # threshold of pi_bm
    VmHat=which(abs(pi_m)<=thresh_m) # valid set based on bm: V_bm
    VFlag[,m]=1:pz %in% VmHat
  }
  
  
  
  
  
  
  ## voting
  
  maxIV=max(apply(VFlag,2,sum))  #the maximum number of valid IVs
  BHatIndex=which(apply(VFlag,2,sum)==maxIV) #give all couples with the maximum number of valid IVs
  
  
  #select
  VHat=(1:pz)[VFlag[,which.max(apply(VFlag,2,sum))]>0]
  
  
  # Warning check: more than one BHatIndex
  if(length(BHatIndex)>1){
    max_columns <- VFlag[, BHatIndex]
    identical_columns <- apply(max_columns, 2, function(col) all(col == max_columns[, 1]))
    all_identical <- all(identical_columns)
    if(!all_identical){
      warning("VHat Warning: More than one cluster has the largest cluster size. This may be due to identification condition not being met. \n")
      #VHat=(1:pz)[apply(VFlag[,BHatIndex],1,sum)>0]
    }
    }
  
  # Error check: not enough IVs
  if(length(VHat) < q){
    warning("VHat Warning: No enough valid IVs estimated. This may be due to weak IVs or identification condition not being met. Use more robust methods.\n")
    warning("Defaulting to all relevant IVs being valid.\n")
    VHat = 1:pz
  }
  
  # ====== Obtain divw est, se, and ci==============================
  
  Sig.V=Sig_w[VHat,VHat]
  eigen_de<-eigen(Sig.V)
  Sig.V.half<-(eigen_de$vectors)%*%diag(sqrt(eigen_de$values))%*%t(eigen_de$vectors)
  
  #est of bias
  Iv=sort(as.vector(outer(VHat, 0:q, function(j,k) j+k*pz)))
  Theta_term<-lapply(VHat, function(j) {
    j0=which(VHat==j)   # 提取 w_i 行向量
    product=kronecker(diag(q+1),Sig.V.half[j0,])
  t(product) %*% Omega_n[Iv,Iv] %*% product
  })
  ThetaSum <- Reduce(`+`, Theta_term)
  
  #est of beta
  denom=t(gamHat[VHat,])%*%Sig.V%*%gamHat[VHat,]-ThetaSum[-1,-1]
  betaHat=solve(denom,t(gamHat[VHat,])%*%Sig.V%*%GamHat[VHat]-ThetaSum[-1,1])
  
  #est of sd
  #Iv=sort(as.vector(outer(VHat, 0:q, function(j,k) j+k*pz)))
  betabar=c(1,-betaHat)
  product2=kronecker(t(betabar),t(gamHat[VHat,])%*%Sig.V)
  num=product2%*%Omega_n[Iv,Iv]%*%t(product2)
  denomInv=solve(denom)
  H=denomInv%*%num%*%denomInv
  sd=sqrt(diag(as.matrix(H)))
  CI=cbind(betaHat-qnorm(1-alpha/2)*sd,betaHat+qnorm(1-alpha/2)*sd) # confidence interval
  
  return( return(list( VHat=VHat, betaHat=betaHat, CI=CI, betaSD=sd, BHatIndex=BHatIndex)))
}



