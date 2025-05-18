# mvmr_voting.R
# Two stage hard thresholding algorithm applied in MVMR individual data when IVs vote for each other
# Select the valid IVs and estimate multiple treatments effects. 
#
# Usage: [VHat,betaHat,ci,betaSD] = mvmr_voting(Y,D,Z,X,intercept,alpha,a0)
#input:
# Y: nx1 outcome vector
# D: nxq treatment matrix
# Z: nxpz candidate IVs
# X: nxpx corvariates
# intercept: whether or not introduce a intercept in linear regression
# alpha: confidence level
# a0: tuning parameter to adjust threshoulding level in first stage; -a0<= beta <= a0
# vote: whether use the voting method

#output:
# VHat: estimated valid IVs
# betaHat: estimated treatments effects
# ci: alpha-level confidence intervals for treatments effects
# betaSD: standard deviation of estimated treatments effects






Voting.summary <- function(GamHat,gamHat,Sig_w,Omega_n,n,a0=1,alpha=0.05) {
  # Include intercept
  pz=length(GamHat)
  q=length(gamHat)/pz
  p=nrow(SigHat)
  # 
  #========= estimate S =================================================
  # gamThresh=sqrt(diag(solve(SigHat))[1:pz] %*% t(diag(ThetaHat[-1,-1])))*sqrt(a0*log(n)/n)
  # #threshold of gam
  # SFlag=(abs(gamHat)>gamThresh)
  # SHat=which(apply(as.matrix(SFlag),1,sum)>0)

  
  SHat=1:pz
  
  #=========== estimate V =================================================
  # B=combn(1:pz,q) # B is a q*C(s,q) matrix
  B=combn(SHat,q)
  # SHat=1:pz
  
 
  SFlag=rep(NA,ncol(B))
  


  for(m in 1:ncol(B)){
    bm=B[,m]
    lambda=svd(gamHat[bm,])$d
    SinMIN=min(lambda) # the minmum singular value
    SinMAX=max(lambda) # the maximum singular value

    SinThresh=2*sqrt(a0*log(n)/n*log(pz))
    SFlag[m]=(SinMIN>SinThresh)

  }

  selected_cols <- which(SFlag == 1)
  
  SHat=sort(unique(as.vector(B[, selected_cols])))
  ##Error check
  if(length(SHat) < q){
    warning("VHat Warning: No enough relevant IVs estimated. This may be due to weak IVs or identification condition not being met. Use more robust methods.\n")
    warning("Defaulting to all IVs being relevant.\n")
    SHat = 1:pz
  }

  
  VFlag=matrix(0,nrow=length(SHat),ncol=ncol(B))
  for(m in selected_cols){
    bm=B[,m]
    beta_m=solve(gamHat[bm,],GamHat[bm])
    pi_m=GamHat[SHat]-gamHat[SHat,]%*%beta_m
    mat_beta=kronecker(t(c(1,-beta_m)),diag(pz)[SHat,]-gamHat[SHat,]%*%solve(gamHat[bm,],diag(pz)[bm,]))
    thresh_m=sqrt(a0*log(n))*sqrt(diag(mat_beta%*%Omega_n%*%t(mat_beta))) # threshold of pi_bm
    VmHat=which(abs(pi_m)<=thresh_m)# valid set based on bm: V_bm
    VFlag[,m]=(1:length(SHat)) %in% VmHat
  }

  
  
  maxIV=max(apply(VFlag,2,sum))  #the maximum number of valid IVs
  BHatIndex=which(apply(VFlag,2,sum)==maxIV) #give all couples with the maximum number of valid IVs
  
  
  ## voting
  # VHat=(1:pz)[VFlag[,which.max(apply(VFlag,2,sum))]>0]
  VHat=SHat[VFlag[,which.max(apply(VFlag,2,sum))]>0]
  
  # Warning check: more than one largest cluster
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
  
  #2SLS
  Sig.V=Sig_w[VHat,VHat]
  denom=t(gamHat[VHat,])%*%Sig.V%*%gamHat[VHat,]
  betaHat=solve(denom,t(gamHat[VHat,])%*%Sig.V%*%GamHat[VHat])
  
  

  #se
  Iv=sort(as.vector(outer(VHat, 0:q, function(j,k) j+k*pz)))
  product=kronecker(t(c(1,-betaHat)),t(gamHat[VHat,])%*%Sig.V)
  denomInv=solve(denom)
  add<-matrix(0,nrow = q,ncol = q)
  for (j in VHat) {
    Ij=outer(j, 0:q, function(j,k) j+k*pz)
    vec=Sig_w[j,j]*(Omega_n[Ij,Ij])[-1,-1]%*%betaHat
    add=add+vec%*%t(vec)
  }
  num=product%*%Omega_n[Iv,Iv]%*%t(product)+add
  H=denomInv%*%num%*%denomInv
  sd=sqrt(diag(as.matrix(H)))
  w=diag(Sig.V)
  if(setequal(VHat,1:pz)){
    res=summary(lm(GamHat ~ 0 + gamHat,  weights = w))
    sd=res$coefficients[,2]
  }
  CI=cbind(betaHat-qnorm(1-alpha/2)*sd,betaHat+qnorm(1-alpha/2)*sd) # confidence interval
  
  
  return(list( VHat=VHat, SHat=SHat, betaHat=betaHat, CI=CI, betaSD=sd, BHatIndex=BHatIndex))
  
}