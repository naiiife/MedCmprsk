## Mediation analysis with competing events ##
## simulation

source("main.R")

gendata = function(n=500, a=NULL, m=NULL, h2=NULL, msample=FALSE){
  X1 = runif(n); X2 = runif(n); X3 = rbinom(n,1,.5); X4 = rbinom(n,1,.5)
  X = cbind(X1, X2, X3, X4)
  if (is.null(a)){
    A = rbinom(n, 1, plogis(-0.2+0.2*X1+0.2*X2+0.2*X3+0.2*X4))
    M = -0.5*X1 + X2 - 0.5*X3 + X4 + A + A*X2 + rnorm(n,1,1)
    h1 = 0.1*X1-0.1*X2-0.1*X3-0.6*M-0.2*A-0.1*M*A
    h2 = 0.1*X1+0.1*X2+0.1*X3-0.8*M-0.2*A-0.1*M*A
    h0 = 0.05*X1+0.05*X3-0.05*X4-0.05*M
    C = 5 + sqrt(-25*log(runif(n))/exp(h0))
    C[C>10] = 10
  } else{
    a0 = a[1]; a1 = a[2]; a2 = a[3]; A = rep(NA,n)
    if (!is.null(m)) {
      M = rep(m, n)
    } else {
      M = -0.5*X1 + X2 - 0.5*X3 + X4 + a0 + a0*X2 + rnorm(n,1,1)
    }
    h1 = 0.1*X1-0.1*X2-0.1*X3-0.6*M-0.2*a1-0.1*M*a1
    if (is.null(h2)) {
      h2 = 0.1*X1+0.1*X2+0.1*X3-0.8*M-0.2*a2-0.1*M*a2
    } else {
      h2 = log(0.25)*rep(1,n)
    }
    C = rep(10,n)
  }
  # cumhaz1: 0.2*t * exp(h1)
  # cumhaz2: 0.8*log(1+t) * exp(h2)
  # cumhaz0: 0.04*max(0,t-5)^2 * exp(h0)
  Time = status = rep(NA,n)
  for (i in 1:n){
    tol = 1; U = -log(runif(1))
    T_l = 0.001; T_u = 10.001; Tk = 4
    EH1 = exp(h1[i]); EH2 = exp(h2[i])
    while (tol>0.00001){
      Q = (0.2*Tk*EH1 + 0.8*log(1+Tk)*EH2)
      if (Q<U) {
        T_l = Tk; Tk = (Tk+T_u)/2; tol = T_u-T_l
      } else {
        T_u = Tk; Tk = (Tk+T_l)/2; tol = T_u-T_l
      }
    }
    prob1 = 0.2 * EH1
    prob2 = 0.8/(1+Tk) * EH2
    Time[i] = min(Tk, C[i])
    status[i] = (1+rbinom(1,1,prob2/(prob1+prob2))) * (C[i]>Tk)
  }
  if (msample){
    DM = rbinom(n, 1, 0.5+0.3*(status>0))
    M[DM==0] = NA
  }
  return(list(X=X,A=A,M=M,Time=Time,status=status))
}


n = 1000
tseq = 1:10; mseq = 0:6/2; h = 0.2

b = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(2026+b)
dat = gendata(n, msample=TRUE)
A = dat$A; X = dat$X; M = dat$M; Time = dat$Time; status = dat$status
fit = try(fit_eif_cf(A, Time, status, M, X, method='riesz'))
if ('try-error' %in% class(fit)) next
h2 = function(tt,x,m) {0.2*rep(1,nrow(x))%*%t(log(1+tt))}
natural = est_decompose(fit)
control = est_con(fit, m_seq=mseq, h=h)
control2 = est_con(fit, m_seq=mseq, h2=h2, h=h)

tt = fit$time
tindex = sapply(1:10, function(k) max(which(tt<=k)))
F_reg = fit$F_reg[,tindex]
F_eff = fit$F_eff[,tindex]
se = fit$se[,tindex]
Fnat_reg = Fnat_eff = senat = matrix(0, nrow=3, ncol=10)
Fnat_reg[1,] = natural$ME_reg[tindex]
Fnat_reg[2,] = natural$DE_reg[tindex]
Fnat_reg[3,] = natural$CE_reg[tindex]
Fnat_eff[1,] = natural$ME_eff[tindex]
Fnat_eff[2,] = natural$DE_eff[tindex]
Fnat_eff[3,] = natural$CE_eff[tindex]
senat[1,] = natural$se_ME[tindex]
senat[2,] = natural$se_DE[tindex]
senat[3,] = natural$se_CE[tindex]
Fcon_reg = control$Fm_reg[,tindex]
Fcon_eff = control$Fm_eff[,tindex]
secon = control$se[,tindex]
Fcon2_reg = control2$Fm_reg[,tindex]
Fcon2_eff = control2$Fm_eff[,tindex]
secon2 = control2$se[,tindex]
  
fit = list(F_reg=F_reg, F_eff=F_eff, se=se, 
           Fnat_reg=Fnat_reg, Fnat_eff=Fnat_eff, senat=senat,
           Fcon_reg=Fcon_reg, Fcon_eff=Fcon_eff, secon=secon, 
           Fcon2_reg=Fcon2_reg, Fcon2_eff=Fcon2_eff, secon2=secon2)
save(fit, file=paste0('fit',b,'.Rdata'))