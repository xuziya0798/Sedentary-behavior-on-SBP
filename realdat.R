library(mr.raps)
library(MendelianRandomization)
library(RobustIV)


set.seed(123)

phenos<-c('tv','wbmi')
phenos.Y<-'sbp'

#standard methods, mv_ivw, mv_egger

re.ivw<-NULL

for(pheno.Y in phenos.Y){
  for(pheno.D in phenos){
  dat<-readRDS(paste(pheno.D,'.',pheno.Y,'.rds',sep='')) 
  p<-length(dat$beta.out)
  
  res_ivw =TwoSampleMR::mr_ivw(dat$beta.exp, dat$beta.out, dat$se.exp, dat$se.out)
  
  re.ivw<-rbind(re.ivw, 
                  c(pheno.D,pheno.Y,p,res_ivw$b,res_ivw$se, 1))
  cat(pheno.D,pheno.Y, p,'\n')}
  
}
re.ivw

re.local<-NULL
tau0=1.6
for(pheno.Y in phenos.Y){
 for(pheno.D in phenos){
  dat<-readRDS(paste(pheno.D,'.',pheno.Y,'.rds',sep='')) 
  p<-length(dat$beta.out)
  
  prop.re =MR.local(beta.out=dat$beta.out, se.out=dat$se.out, 
                    beta.exp=dat$beta.exp, se.exp=dat$se.exp,
                    tau0=tau0,tau.s=tau0,skew.ind = F)

  re.local<-rbind(re.local, 
               c(pheno.D,pheno.Y,length(prop.re$clust.sel),prop.re$beta.divw,prop.re$sd.divw, prop.re$balanced.pleiotropy))
  cat(pheno.D,pheno.Y, p,'\n')
}
}
re.local

re.local0<-NULL
tau0=1.6

for(pheno.Y in phenos.Y){
 for(pheno.D in phenos){
  dat<-readRDS(paste(pheno.D,'.',pheno.Y,'.rds',sep='')) 
  p<-length(dat$beta.out)
  
  prop.re =MR.local(beta.out=dat$beta.out, se.out=dat$se.out, 
                    beta.exp=dat$beta.exp, se.exp=dat$se.exp,
                    tau0=tau0,tau.s=tau0,skew.ind = T)
  
  re.local0<-rbind(re.local0, 
                  c(pheno.D,pheno.Y,length(prop.re$clust.sel),prop.re$beta.divw,prop.re$sd.divw, prop.re$balanced.pleiotropy))
  cat(pheno.D,pheno.Y, p,'\n')
 }
}

re.local0


#mr.raps
re.raps<-NULL
for (pheno.Y in phenos.Y){
  for (pheno.D in phenos){
    dat<-readRDS(paste(pheno.D,'.',pheno.Y,'.rds',sep=''))  
    p<-length(dat$beta.out)
    raps.re=mr.raps(b_out=dat$beta.out, se_out=dat$se.out, 
                    b_exp=dat$beta.exp, se_exp=dat$se.exp,
                    over.dispersion = T)
    
    re.raps<-rbind(re.raps, 
                 c(pheno.D,pheno.Y,p,raps.re$beta.hat,raps.re$beta.se,1))
    cat(pheno.D,pheno.Y, p,'\n')
  }
}
re.raps

#mr.mbe
re.mbe<-NULL
for(pheno.Y in phenos.Y){
for(pheno.D in phenos){
  dat<-readRDS(paste(pheno.D,'.',pheno.Y,'.rds',sep=''))  
  p<-length(dat$beta.out)
  Sd.hat<-which(abs(dat$beta.exp/dat$se.exp)>=sqrt(2*log(p)))
  Mbe.re<-MendelianRandomization::mr_mbe(mr_input(bx = dat$beta.exp[Sd.hat], bxse = dat$se.exp[Sd.hat],
                                                  by = dat$beta.out[Sd.hat], byse = dat$se.out[Sd.hat]), 
                                         iterations=1000,weighting='unweighted')
  
  re.mbe<-rbind(re.mbe, 
               c(pheno.D,pheno.Y,Mbe.re@SNPs,Mbe.re@Estimate,Mbe.re@StdError,0))
  cat(pheno.D,pheno.Y, p,'\n')
}
}
re.mbe



re.2sample<-data.frame(rbind(re.ivw,re.raps,re.local,re.local0,re.mbe))
colnames(re.2sample)<-c('treatment','outcome', 'SNPs','est','sd','bp')



#add ci
z_score <- qnorm(0.975) 

re.2sample$CI_lower <- as.numeric(re.2sample$est) - z_score * as.numeric(re.2sample$sd)
re.2sample$CI_upper <- as.numeric(re.2sample$est) + z_score * as.numeric(re.2sample$sd)

#add p.val
z.value<-as.numeric(re.2sample$est)/as.numeric(re.2sample$sd)
re.2sample$p.value<-2 * (1 - pnorm(abs(z.value)))


methods<-c('ivw','raps','local','local0','mbe')
re.2sample$method <-rep(methods,each=4)
re.2sample

write.csv(re.2sample,"re.2sample.2out.csv",row.names = FALSE)





######## mvmr ################
mvdat<-readRDS('tv_wbmi.dbp.rds')


#mv-ivw
res.mvivw<-mv_multiple(dat, method = "IVW", se = "classic")
re.mvivw<-cbind(res.mvivw$result$b,
                res.mvivw$result$se,res.mvivw$result$nsnp)
#instrument_specific	using all instruments from all exposures (FALSE, default) 
re.mvivw


set.seed(123)
# mvmr-local


n=4e5

svd(mvdat$exposure_beta/mvdat$outcome_se)$d
vec_se<-c(mvdat$outcome_se,as.vector(mvdat$exposure_se))
Omega_n=diag(vec_se^2)
Sig_w=diag(1/mvdat$outcome_se^2)



res.mvlocal=Local.summary(GamHat=as.matrix(mvdat$outcome_beta),gamHat=mvdat$exposure_beta,
              Sig_w=Sig_w,Omega_n=Omega_n,n=n,sep=0.05) 
re.mvlocal<-cbind(res.mvlocal$betaHat,res.mvlocal$betaSD,rep(length(res.mvlocal$VHat),2))
re.mvlocal

res.voting=Voting.summary(GamHat=as.matrix(mvdat$outcome_beta),gamHat=mvdat$exposure_beta,
                          Sig_w=Sig_w,Omega_n=Omega_n,n=n) 
re.voting<-cbind(res.voting$betaHat,res.voting$betaSD,rep(length(res.voting$VHat),2))
re.voting

# mvmr-robust
res_robust = mvmr_robust(mvdat$exposure_beta, mvdat$outcome_beta,
                         mvdat$outcome_se, k.max = 1000, maxit.scale = 1000)
re.robust<-cbind(res_robust$coefficients,res_robust$se,rep(NA,2))
re.robust
# mvmr-lasso
res_lasso = mvmr_lasso(mvdat$exposure_beta, mvdat$outcome_beta,
                       mvdat$outcome_se)
re.lasso<-cbind(res_lasso$est_post,res_lasso$se_post,rep(length(res_lasso$v)))
re.lasso

# mr-median
res_median = mvmr_median(mvdat$exposure_beta,mvdat$exposure_se, mvdat$outcome_beta,
                         mvdat$outcome_se,boot = TRUE)
re.median<-cbind(res_median$coefficients,res_median$se,rep(NA,2))
re.median


##summary
re.3sample<-data.frame(rbind(re.mvlocal,re.voting,re.median))
colnames(re.3sample)<-c('est','sd','SNPs')


#add ci
z_score <- qnorm(0.975) 

re.3sample$CI_lower <- as.numeric(re.3sample$est) - z_score * as.numeric(re.3sample$sd)
re.3sample$CI_upper <- as.numeric(re.3sample$est) + z_score * as.numeric(re.3sample$sd)

#add p.val
z.value<-as.numeric(re.3sample$est)/as.numeric(re.3sample$sd)
re.3sample$p.value<-2 * (1 - pnorm(abs(z.value)))

methods<-c('mvlocal','voting','mvivw','robust',"median")
re.3sample$method <-rep(methods,each=2)

re.3sample

write.csv(re.3sample,"re.3sample.dbp.csv",row.names = TRUE)

