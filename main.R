## Mediation analysis with competing events ##

library(statmod)
library(survival)

matchy = function(x, y, newx, exact=FALSE){
  if (exact) {
    newy = sapply(newx, function(l) y[max(which(x==l))])
  } else {
    newy = sapply(newx, function(l) y[max(which(x<=l))])
  }
  newy = as.numeric(newy)
  newy[is.na(newy)] = 0
  return(newy)
}

ngrid = 10
quad = gauss.quad(ngrid, "hermite")
b_list = quad$nodes
w_list = quad$weights/sqrt(pi)

riesz_basis <- function(M, X) {
  X <- as.matrix(X)
  cbind(`(Intercept)` = 1, M = M, M2 = M^2,
    X, MX = M * X, X2 = X^2)
}

fit_riesz_r10 <- function(M, X, A, D, wts, ps, train, lambda = 0.01) {
  X <- as.matrix(X)
  D <- as.integer(D)
  keep1 <- train & D == 1 & A == 1
  keep0 <- train & D == 1 & A == 0
  Phi1 <- riesz_basis(M[keep1], X[keep1, , drop = FALSE])
  Phi0 <- riesz_basis(M[keep0], X[keep0, , drop = FALSE])
  w1 <- wts[keep1] / ps[keep1]
  w0 <- wts[keep0] / ps[keep0]
  w1 <- w1 / sum(w1)
  w0 <- w0 / sum(w0)
  H <- crossprod(Phi0, Phi0 * w0)
  h <- colSums(Phi1 * w1)
  d <- ncol(Phi0)
  theta <- solve(H + lambda * diag(d), h)
  list(
    theta = theta,
    lambda = lambda,
    basis_names = colnames(Phi0),
    weight_summary = list(
      n_A1 = sum(keep1),
      n_A0 = sum(keep0),
      ess_A1 = sum(w1)^2 / sum(w1^2),
      ess_A0 = sum(w0)^2 / sum(w0^2)
    )
  )
}

predict_riesz_r10 <- function(fit, M_new, X_new, lower = 0.01, upper = 100) {
  rhat <- drop(riesz_basis(M_new, X_new) %*% fit$theta)
  rhat <- pmin(pmax(rhat, lower), upper)
  rhat
}

fit_eif_cf <- function(A, Time, status, M, X, method='riesz', event=1, monotone=TRUE){
  # fit proportional hazard models
  if (event!=1){
    status.pr = as.numeric(status==event)
    status.cp = as.numeric(status!=event & status!=0)
    status = status.pr + 2*status.cp
  }
  Time = round(Time, 4)
  dat = data.frame(M,A,X,Time,status)
  tt = sort(unique(c(0,Time)))
  mt = length(tt)
  D = !is.na(M)
  n = length(A)
  ss = (1:n) %in% sample(n, n/2)
  
  ## estimate sampling weight
  wts = rep(1,n)
  if (sum(D)<n){
    iwts = glm(D~., data=data.frame(A,X,Time,st=as.factor(status)), family='binomial', subset=ss)
    iwts = predict(iwts, newdata=data.frame(A,X,Time,st=as.factor(status))[!ss,], type='response')
    wts[!ss] = D[!ss]/iwts
    iwts = glm(D~., data=data.frame(A,X,Time,st=as.factor(status)), family='binomial', subset=!ss)
    iwts = predict(iwts, newdata=data.frame(A,X,Time,st=as.factor(status))[ss,], type='response')
    wts[ss] = D[ss]/iwts
  }
  wts = pmin(100, pmax(wts, 0.01))
  wts = wts * D
  
  ## fit nuisance models
  fit11_cf = fit01_cf = fit12_cf = fit02_cf = fit10_cf = fit00_cf = list(NA,NA)
  lam11_cf = lam01_cf = lam12_cf = lam02_cf = lam10_cf = lam00_cf = list(NA,NA)
  # fit model on ss, predict on 1ss
  {
    fit11 = coxph(Surv(Time, status==1) ~ ., data=data.frame(M,X), weights=wts, subset=(A==1&D==1&ss))
    fit01 = coxph(Surv(Time, status==1) ~ ., data=data.frame(M,X), weights=wts, subset=(A==0&D==1&ss))
    fit12 = coxph(Surv(Time, status==2) ~ ., data=data.frame(M,X), weights=wts, subset=(A==1&D==1&ss))
    fit02 = coxph(Surv(Time, status==2) ~ ., data=data.frame(M,X), weights=wts, subset=(A==0&D==1&ss))
    fit10 = coxph(Surv(Time, status==0) ~ ., data=data.frame(M,X), weights=wts, subset=(A==1&D==1&ss))
    fit00 = coxph(Surv(Time, status==0) ~ ., data=data.frame(M,X), weights=wts, subset=(A==0&D==1&ss))
    haz11 = rbind(c(0,0), basehaz(fit11, centered = FALSE))
    haz01 = rbind(c(0,0), basehaz(fit01, centered = FALSE))
    haz12 = rbind(c(0,0), basehaz(fit12, centered = FALSE))
    haz02 = rbind(c(0,0), basehaz(fit02, centered = FALSE))
    haz10 = rbind(c(0,0), basehaz(fit10, centered = FALSE))
    haz00 = rbind(c(0,0), basehaz(fit00, centered = FALSE))
    lam11 = diff(c(0, matchy(haz11$time, haz11$hazard, tt)))
    lam01 = diff(c(0, matchy(haz01$time, haz01$hazard, tt)))
    lam12 = diff(c(0, matchy(haz12$time, haz12$hazard, tt)))
    lam02 = diff(c(0, matchy(haz02$time, haz02$hazard, tt)))
    lam10 = diff(c(0, matchy(haz10$time, haz10$hazard, tt)))
    lam00 = diff(c(0, matchy(haz00$time, haz00$hazard, tt)))
    Xb11 = predict(fit11, data.frame(M,X)[!ss,], type='lp', reference='zero')
    Xb01 = predict(fit01, data.frame(M,X)[!ss,], type='lp', reference='zero')
    Xb12 = predict(fit12, data.frame(M,X)[!ss,], type='lp', reference='zero')
    Xb02 = predict(fit02, data.frame(M,X)[!ss,], type='lp', reference='zero')
    Xb10 = predict(fit10, data.frame(M,X)[!ss,], type='lp', reference='zero')
    Xb00 = predict(fit00, data.frame(M,X)[!ss,], type='lp', reference='zero')
    haz11_1ss = exp(Xb11) %*% t(lam11)
    haz01_1ss = exp(Xb01) %*% t(lam01)
    haz12_1ss = exp(Xb12) %*% t(lam12)
    haz02_1ss = exp(Xb02) %*% t(lam02)
    haz10_1ss = exp(Xb10) %*% t(lam10)
    haz00_1ss = exp(Xb00) %*% t(lam00)
    fit11_cf[[1]] = fit11; fit12_cf[[1]] = fit12; fit10_cf[[1]] = fit10
    fit01_cf[[1]] = fit01; fit02_cf[[1]] = fit02; fit00_cf[[1]] = fit00
    lam11_cf[[1]] = lam11; lam12_cf[[1]] = lam12; lam10_cf[[1]] = lam10
    lam01_cf[[1]] = lam01; lam02_cf[[1]] = lam02; lam00_cf[[1]] = lam00
  }
  # fit model on 1ss, predict on ss
  {
    fit11 = coxph(Surv(Time, status==1) ~ ., data=data.frame(M,X), weights=wts, subset=(A==1&D==1&!ss))
    fit01 = coxph(Surv(Time, status==1) ~ ., data=data.frame(M,X), weights=wts, subset=(A==0&D==1&!ss))
    fit12 = coxph(Surv(Time, status==2) ~ ., data=data.frame(M,X), weights=wts, subset=(A==1&D==1&!ss))
    fit02 = coxph(Surv(Time, status==2) ~ ., data=data.frame(M,X), weights=wts, subset=(A==0&D==1&!ss))
    fit10 = coxph(Surv(Time, status==0) ~ ., data=data.frame(M,X), weights=wts, subset=(A==1&D==1&!ss))
    fit00 = coxph(Surv(Time, status==0) ~ ., data=data.frame(M,X), weights=wts, subset=(A==0&D==1&!ss))
    haz11 = rbind(c(0,0), basehaz(fit11, centered = FALSE))
    haz01 = rbind(c(0,0), basehaz(fit01, centered = FALSE))
    haz12 = rbind(c(0,0), basehaz(fit12, centered = FALSE))
    haz02 = rbind(c(0,0), basehaz(fit02, centered = FALSE))
    haz10 = rbind(c(0,0), basehaz(fit10, centered = FALSE))
    haz00 = rbind(c(0,0), basehaz(fit00, centered = FALSE))
    lam11 = diff(c(0, matchy(haz11$time, haz11$hazard, tt)))
    lam01 = diff(c(0, matchy(haz01$time, haz01$hazard, tt)))
    lam12 = diff(c(0, matchy(haz12$time, haz12$hazard, tt)))
    lam02 = diff(c(0, matchy(haz02$time, haz02$hazard, tt)))
    lam10 = diff(c(0, matchy(haz10$time, haz10$hazard, tt)))
    lam00 = diff(c(0, matchy(haz00$time, haz00$hazard, tt)))
    Xb11 = predict(fit11, data.frame(M,X)[ss,], type='lp', reference='zero')
    Xb01 = predict(fit01, data.frame(M,X)[ss,], type='lp', reference='zero')
    Xb12 = predict(fit12, data.frame(M,X)[ss,], type='lp', reference='zero')
    Xb02 = predict(fit02, data.frame(M,X)[ss,], type='lp', reference='zero')
    Xb10 = predict(fit10, data.frame(M,X)[ss,], type='lp', reference='zero')
    Xb00 = predict(fit00, data.frame(M,X)[ss,], type='lp', reference='zero')
    haz11_ss = exp(Xb11) %*% t(lam11)
    haz01_ss = exp(Xb01) %*% t(lam01)
    haz12_ss = exp(Xb12) %*% t(lam12)
    haz02_ss = exp(Xb02) %*% t(lam02)
    haz10_ss = exp(Xb10) %*% t(lam10)
    haz00_ss = exp(Xb00) %*% t(lam00)
    fit11_cf[[2]] = fit11; fit12_cf[[2]] = fit12; fit10_cf[[2]] = fit10
    fit01_cf[[2]] = fit01; fit02_cf[[2]] = fit02; fit00_cf[[2]] = fit00
    lam11_cf[[2]] = lam11; lam12_cf[[2]] = lam12; lam10_cf[[2]] = lam10
    lam01_cf[[2]] = lam01; lam02_cf[[2]] = lam02; lam00_cf[[2]] = lam00
  }
  haz11 = haz01 = haz12 = haz02 = haz10 = haz00 = matrix(0,nrow=n,ncol=mt)
  haz11[ss,] = haz11_ss; haz11[!ss,] = haz11_1ss
  haz01[ss,] = haz01_ss; haz01[!ss,] = haz01_1ss
  haz12[ss,] = haz12_ss; haz12[!ss,] = haz12_1ss
  haz02[ss,] = haz02_ss; haz02[!ss,] = haz02_1ss
  haz10[ss,] = haz10_ss; haz10[!ss,] = haz10_1ss
  haz00[ss,] = haz00_ss; haz00[!ss,] = haz00_1ss
  
  # fit propensity score
  ps_a = rep(NA, n)
  fit_a = glm(A~., data=data.frame(A,X), family='binomial', subset=ss)
  ps_a[!ss] = predict(fit_a, newdata=data.frame(A,X)[!ss,], type='response')
  fit_a = glm(A~., data=data.frame(A,X), family='binomial', subset=!ss)
  ps_a[ss] = predict(fit_a, newdata=data.frame(A,X)[ss,], type='response')
  ps = ps_a*A + (1-ps_a)*(1-A)
  ps[A==1] = ps[A==1] * mean(A/ps)
  ps[A==0] = ps[A==0] * mean((1-A)/ps)
  
  # fit mediator density
  if (method=='ratio'){
    mx = lm(M ~ ., data=data.frame(A,X,AX=A*X), weights=wts, subset=(D==1&ss))
    res_m = rep(NA, n)
    res_m[D==1&ss] = residuals(mx)^2
    mx1 = predict(mx, newdata=data.frame(A=1,X,AX=X)[!ss,])
    mx0 = predict(mx, newdata=data.frame(A=0,X,AX=0*X)[!ss,])
    sx = lm(log(res_m+0.000001) ~ ., data=data.frame(A,X,AX=A*X), weights=wts, subset=(D==1&ss))
    s2 = mean(residuals(sx)^2)
    sx1 = exp(predict(sx, newdata=data.frame(A=1,X,AX=X)[!ss,]) + s2/2)
    sx0 = exp(predict(sx, newdata=data.frame(A=0,X,AX=0*X)[!ss,]) + s2/2)
    mx_1ss = cbind(mx0, mx1)
    sx_1ss = sqrt(cbind(sx0, sx1))
    mx = lm(M ~ ., data=data.frame(A,X,AX=A*X), weights=wts, subset=(D==1&!ss))
    res_m = rep(NA, n)
    res_m[D==1&!ss] = residuals(mx)^2
    mx1 = predict(mx, newdata=data.frame(A=1,X,AX=X)[ss,])
    mx0 = predict(mx, newdata=data.frame(A=0,X,AX=0*X)[ss,])
    sx = lm(log(res_m+0.000001) ~ ., data=data.frame(A,X,AX=A*X), weights=wts, subset=(D==1&!ss))
    s2 = mean(residuals(sx)^2)
    sx1 = exp(predict(sx, newdata=data.frame(A=1,X,AX=X)[ss,]) + s2/2)
    sx0 = exp(predict(sx, newdata=data.frame(A=0,X,AX=0*X)[ss,]) + s2/2)
    mx_ss = cbind(mx0, mx1)
    sx_ss = sqrt(cbind(sx0, sx1))
    mx = sx = matrix(0, nrow=n, ncol=2)
    mx[ss,] = mx_ss; mx[!ss,] = mx_1ss
    sx[ss,] = sx_ss; sx[!ss,] = sx_1ss
  } else if (method=='propensity') {
    ps_m = rep(NA, n)
    fit_m = glm(A~., data=data.frame(M,X,MX=M*X), family='binomial', weights=wts, subset=(D==1&ss))
    ps_m[D==1&!ss] = predict(fit_m, newdata=data.frame(M,X,MX=M*X)[D==1&!ss,], type='response')
    fit_m = glm(A~., data=data.frame(M,X,MX=M*X), family='binomial', weights=wts, subset=(D==1&!ss))
    ps_m[D==1&ss] = predict(fit_m, newdata=data.frame(M,X,MX=M*X)[D==1&ss,], type='response')
    ps_a = cbind(1-ps_a, ps_a)
    ps_m = cbind(1-ps_m, ps_m)/ps_a
  } else if (method=='riesz'){
    r10 = rep(NA, n)
    fit_r10_ss = fit_riesz_r10(M, X, A, D, wts, ps, ss)
    r10[!ss] = predict_riesz_r10(fit_r10_ss, M_new=M[!ss], X_new=X[!ss,,drop = FALSE])
    fit_r10_1ss = fit_riesz_r10(M, X, A, D, wts, ps, !ss)
    r10[ss] = predict_riesz_r10(fit_r10_1ss, M_new=M[ss], X_new=X[ss,,drop = FALSE])
  }
  
  # observed survival functions
  haz1 = haz11*A + haz01*(1-A)
  haz2 = haz12*A + haz02*(1-A)
  haz0 = cbind(0, haz10*A + haz00*(1-A))[,1:ncol(haz00)]
  S = t(exp(apply(-haz1-haz2-haz0, 1, cumsum)))
  S[S<0.01] = 0.01
  Yt = sapply(tt, function(l) Time>=l)
  dNt = sapply(tt, function(l) Time==l)
  dG1t = (dNt * (status==1) - Yt * haz1) / S
  dG2t = (dNt * (status==2) - Yt * haz2) / S
  G1t = t(apply(dG1t, 1, cumsum))
  G2t = t(apply(dG2t, 1, cumsum))
  
  # EIF
  # a = (a0, a1, a2)
  n = length(A)
  a_seq = rbind(c(0,0,0), c(1,0,0), c(1,1,0), c(1,1,1))
  F_reg = F_eff = se_eff = matrix(NA, nrow=4, ncol=length(tt))
  a_names = apply(a_seq, 1, paste0, collapse="")
  rownames(F_reg) = rownames(F_eff) = rownames(se_eff) = a_names
  eif = array(NA, dim=c(n,length(tt),4), dimnames=list(NULL,round(tt,4),a_names))
  for (k in 1:4){
    a = a_seq[k,]
    a0 = a[1]; a1 = a[2]; a2 = a[3]
    if (a1==1){
      fit1 = fit11_cf; haz1 = haz11; lam1 = lam11_cf
    } else {
      fit1 = fit01_cf; haz1 = haz01; lam1 = lam01_cf
    }
    if (a2==1){
      fit2 = fit12_cf; haz2 = haz12; lam2 = lam12_cf
    } else {
      fit2 = fit02_cf; haz2 = haz02; lam2 = lam02_cf
    }
    S = exp(t(apply(-haz1-haz2, 1, cumsum)))
    rm01 = rm02 = rep(1, n)
    if (method=='ratio'){
      if (a0!=a1) rm01 = dnorm(M,mx[,a0+1],sx[,a0+1])/dnorm(M,mx[,a1+1],sx[,a1+1])
      if (a0!=a2) rm02 = dnorm(M,mx[,a0+1],sx[,a0+1])/dnorm(M,mx[,a2+1],sx[,a2+1])
    } else if (method=='propensity') {
      if (a0!=a1) rm01 = ps_m[,a0+1]/ps_m[,a1+1]
      if (a0!=a2) rm02 = ps_m[,a0+1]/ps_m[,a2+1]
    } else if (method=='riesz'){
      if (a0!=a1) rm01 = (a1==1)/r10 + (a1==0)*r10
      if (a0!=a2) rm02 = (a2==1)/r10 + (a2==0)*r10
    }
    rm01[rm01<0.01] = 0.01; rm02[rm02<0.01] = 0.01
    rm01[rm01>99] = 99; rm02[rm02>99] = 99
    term1 = (A==a1)/ps * rm01 * (dG1t - G1t * haz1)
    term2 = -(A==a2)/ps * rm02 * G2t * haz1
    term3 = (A==a0)/ps * haz1
    
    idx0_ss  = (A==a0 & D==1 & ss)
    idx0_1ss = (A==a0 & D==1 & !ss)
    gamma_obs = matrix(NA, nrow=n, ncol=mt)
    h_hat = matrix(0, nrow=n, ncol=mt)
    Xdf = data.frame(X)
    Xb1o = predict(fit1[[1]], data.frame(M,X)[idx0_ss,], type='lp', reference='zero')
    Xb2o = predict(fit2[[1]], data.frame(M,X)[idx0_ss,], type='lp', reference='zero')
    h1o = exp(Xb1o) %*% t(lam1[[1]])
    h2o = exp(Xb2o) %*% t(lam2[[1]])
    gamma_obs[idx0_ss,] = exp(t(apply(-h1o-h2o, 1, cumsum))) * h1o
    fit_h = lm(gamma_obs ~ ., data=Xdf, weights=wts, subset=idx0_ss)
    h_hat[!ss,] = predict(fit_h, newdata=Xdf[!ss,,drop=FALSE])
    Xb1o = predict(fit1[[2]], data.frame(M,X)[idx0_1ss,], type='lp', reference='zero')
    Xb2o = predict(fit2[[2]], data.frame(M,X)[idx0_1ss,], type='lp', reference='zero')
    h1o = exp(Xb1o) %*% t(lam1[[2]])
    h2o = exp(Xb2o) %*% t(lam2[[2]])
    gamma_obs[idx0_1ss,] = exp(t(apply(-h1o-h2o, 1, cumsum))) * h1o
    fit_h = lm(gamma_obs ~ ., data=Xdf, weights=wts, subset=idx0_1ss)
    h_hat[ss,] = predict(fit_h, newdata=Xdf[ss,,drop=FALSE])
    h_hat[h_hat<0] = 0
    dF = h_hat
    F_reg[k,] = rowMeans(apply(dF, 1, cumsum), na.rm=TRUE)
    term4 = (1-(A==a0)/ps) * dF
    if (a0==0&a1==0&a2==0){
      phi0 = term1 + term2; phi00 = term1
    }
    if (a0==1&a1==1&a2==1){
      phi1 = term1 + term2; phi10 = term1
    }
    dF = S * (term1 + term2 + term3) + term4
    Fx = t(apply(dF, 1, cumsum))
    Fx[D==0,] = 0
    EFx = matrix(0,nrow=n,ncol=mt)
    fit_phi = lm(Fx ~ ., data=data.frame(A,X,Time,st=as.factor(status)), weights=wts, subset=(D==1&ss))
    EFx[!ss,] = predict(fit_phi, newdata=data.frame(A,X,Time,st=as.factor(status))[!ss,])
    fit_phi = lm(Fx ~ ., data=data.frame(A,X,Time,st=as.factor(status)), weights=wts, subset=(D==1&!ss))
    EFx[ss,] = predict(fit_phi, newdata=data.frame(A,X,Time,st=as.factor(status))[ss,])
    Fx = wts*Fx + (1-wts)*EFx
    F_effk = pmin(1,pmax(0,colMeans(Fx, na.rm=TRUE)))
    if (monotone){
      F_effk = isoreg(F_effk)$yf
    } 
    eifk = t(t(Fx) - F_effk)
    se_effk = sqrt(apply(eifk^2, 2, mean, na.rm=TRUE)) / sqrt(n-1)
    F_eff[k,] = F_effk
    se_eff[k,] = se_effk
    eif[,,k] = eifk
  }
  colnames(F_reg) = colnames(F_eff) = round(tt,4)
  
  return(list(data=dat, time=tt, phi0=phi0, phi1=phi1, phi00=phi00, phi10=phi10, 
              eif=eif, wts=wts, ps=ps,
              lam01=lam01_cf, lam02=lam02_cf, fit01=fit01_cf, fit02=fit02_cf,
              lam11=lam11_cf, lam12=lam12_cf, fit11=fit11_cf, fit12=fit12_cf,
              F_reg=F_reg, F_eff=F_eff, se=se_eff, cf=TRUE, ss=ss))
}


## Estimate decomposed effects
est_decompose <- function(fit){
  n = nrow(fit$data)
  tt = fit$time
  ctt = diff(c(0,tt))
  eif = fit$eif
  F_reg = fit$F_reg
  F_eff = fit$F_eff
  TE_reg = F_reg["111",] - F_reg["000",]
  TE_eff = F_eff["111",] - F_eff["000",]
  eif_TE = eif[,,"111"]-eif[,,"000"]
  se_TE = apply(eif_TE, 2, sd, na.rm=TRUE) / sqrt(n)
  RMSTL_TE = sum(TE_eff*ctt, na.rm=TRUE)
  RMSTL_TE_se = sd(colSums(t(eif_TE)*ctt, na.rm=TRUE))/sqrt(n)
  p_TE = 2*pnorm(-abs(RMSTL_TE/RMSTL_TE_se))
  ME_reg = F_reg["100",] - F_reg["000",]
  ME_eff = F_eff["100",] - F_eff["000",]
  eif_ME = eif[,,"100"]-eif[,,"000"]
  se_ME = apply(eif_ME, 2, sd, na.rm=TRUE) / sqrt(n)
  RMSTL_ME = sum(ME_eff*ctt, na.rm=TRUE)
  RMSTL_ME_se = sd(colSums(t(eif_ME)*ctt, na.rm=TRUE))/sqrt(n)
  p_ME = 2*pnorm(-abs(RMSTL_ME/RMSTL_ME_se))
  DE_reg = F_reg["110",] - F_reg["100",]
  DE_eff = F_eff["110",] - F_eff["100",]
  eif_DE = eif[,,"110"]-eif[,,"100"]
  se_DE = apply(eif_DE, 2, sd, na.rm=TRUE) / sqrt(n)
  RMSTL_DE = sum(DE_eff*ctt, na.rm=TRUE)
  RMSTL_DE_se = sd(colSums(t(eif_DE)*ctt, na.rm=TRUE))/sqrt(n)
  p_DE = 2*pnorm(-abs(RMSTL_DE/RMSTL_DE_se))
  CE_reg = F_reg["111",] - F_reg["110",]
  CE_eff = F_eff["111",] - F_eff["110",]
  eif_CE = eif[,,"111"]-eif[,,"110"]
  se_CE = apply(eif_CE, 2, sd, na.rm=TRUE) / sqrt(n)
  RMSTL_CE = sum(CE_eff*ctt, na.rm=TRUE)
  RMSTL_CE_se = sd(colSums(t(eif_CE)*ctt, na.rm=TRUE))/sqrt(n)
  p_CE = 2*pnorm(-abs(RMSTL_CE/RMSTL_CE_se))
  
  return(list(time=tt, F_reg=F_reg, F_eff=F_eff,
              TE_reg=TE_reg, TE_eff=TE_eff, se_TE=se_TE, RMSTL_TE=RMSTL_TE, p_TE=p_TE,
              ME_reg=ME_reg, ME_eff=ME_eff, se_ME=se_ME, RMSTL_ME=RMSTL_ME, p_ME=p_ME,
              DE_reg=DE_reg, DE_eff=DE_eff, se_DE=se_DE, RMSTL_DE=RMSTL_DE, p_DE=p_DE,
              CE_reg=CE_reg, CE_eff=CE_eff, se_CE=se_CE, RMSTL_CE=RMSTL_CE, p_CE=p_CE))
}

GK <- function(u,h=1) {
  (1/sqrt(2*pi)) * exp(-((u/h)^2)/2) / h
}

fit_riesz_controlled <- function(m0, h, M, X, A, D, wts, ps, a, train, lambda = 0.01) {
  X <- as.matrix(X)
  D <- as.integer(D)
  keep_a <- train & D == 1 & A == a
  Phi_a <- riesz_basis(M[keep_a], X[keep_a,,drop = FALSE])
  w_a <- wts[keep_a] / ps[keep_a]
  w_a <- w_a / sum(w_a)
  H <- crossprod(Phi_a, Phi_a * w_a)
  g <- rep(0, ncol(Phi_a))
  X_a <- X[keep_a, , drop = FALSE]
  for (b in 1:ngrid) {
    m_b <- m0 + sqrt(2) * h * b_list[b]
    Phi_b <- riesz_basis(rep(m_b, nrow(X_a)), X_a)
    g <- g + w_list[b] * colSums(Phi_b * w_a)
  }
  d <- ncol(Phi_a)
  theta <- solve(H + lambda * diag(d), g)
  list(
    theta = theta,
    lambda = lambda,
    a = a,
    m0 = m0,
    h = h,
    n_train = sum(keep_a),
    ess_train = sum(w_a)^2 / sum(w_a^2)
  )
}

predict_riesz_controlled <- function(fit, M_new, X_new, lower = 1e-5, upper = 99) {
  alpha_hat <- drop(riesz_basis(M_new, X_new) %*% fit$theta)
  alpha_hat <- pmin(pmax(alpha_hat, lower), upper)
  alpha_hat
}

est_con <- function(fit, a=0, m_seq=NULL, h2=NULL, time=NULL, h=0.5, monotone=TRUE){
  n = nrow(fit$data)
  eif = fit$eif; wts = fit$wts; ps = fit$ps; ss = fit$ss
  tt = fit$time; mt = length(tt)
  if (a==0){
    fit1 = fit$fit01; fit2 = fit$fit02
    lam1 = fit$lam01; lam2 = fit$lam02
    phi0 = fit$phi0; phi00 = fit$phi00
  } else if (a==1) {
    fit1 = fit$fit11; fit2 = fit$fit12
    lam1 = fit$lam11; lam2 = fit$lam12
    phi0 = fit$phi1; phi00 = fit$phi10
  }
  A = fit$data$A
  X = as.matrix(fit$data[, grepl("^X", names(fit$data)) | names(fit$data) == "X"])
  M = fit$data$M; D = !is.na(M)
  if (is.null(m_seq)){
    M_quantile = quantile(na.omit(fit$data$M), c(0.1, 0.9), na.rm=TRUE)
    m_seq = seq(M_quantile[1], M_quantile[2], length = 11)
  }
  m_length = length(m_seq)
  Fm_reg = Fm_eff = se = matrix(NA, ncol=length(tt), nrow=m_length)
  idx0_ss  = which(A==a & D==1 & ss)
  idx0_1ss = which(A==a & D==1 & !ss)
  for (k in 1:m_length){
    m = m_seq[k]
    alpha0 = rep(0, n)
    fit_alpha_ss = fit_riesz_controlled(m, h, M, X, A, D, wts, ps, a, ss)
    alpha0[!ss] = predict_riesz_controlled(fit_alpha_ss, M_new=M[!ss], X_new=X[!ss,,drop=FALSE])
    fit_alpha_1ss = fit_riesz_controlled(m, h, M, X, A, D, wts, ps, a, !ss)
    alpha0[ss] = predict_riesz_controlled(fit_alpha_1ss, M_new=M[ss], X_new=X[ss,,drop=FALSE])
    haz1 = haz2 = matrix(0, nrow=n, ncol=mt)
    Xb1 = predict(fit1[[1]], data.frame(M,X)[!ss,], type='lp', reference='zero')
    haz1[!ss,] = exp(Xb1) %*% t(lam1[[1]])
    Xb1 = predict(fit1[[2]], data.frame(M,X)[ss,], type='lp', reference='zero')
    haz1[ss,] = exp(Xb1) %*% t(lam1[[2]])
    if (is.null(h2)){
      Xb2 = predict(fit2[[1]], data.frame(M,X)[!ss,], type='lp', reference='zero')
      haz2[!ss,] = exp(Xb2) %*% t(lam2[[1]])
      Xb2 = predict(fit2[[2]], data.frame(M,X)[ss,], type='lp', reference='zero')
      haz2[ss,] = exp(Xb2) %*% t(lam2[[2]])
      term1 = alpha0 * phi0
    } else {
      haz2 = t(apply(cbind(0,h2(tt,X,M)),1,diff))
      term1 = alpha0 * phi00
    }
    term1 = term1 * exp(t(apply(-haz1-haz2, 1, cumsum)))
    term2 = 0
    for (b in 1:ngrid){
      m_b = m + sqrt(2)*h*b_list[b]
      haz1 = haz2 = matrix(0, nrow=n, ncol=mt)
      Xb1 = predict(fit1[[1]], data.frame(M=m_b, X)[!ss,], type='lp', reference='zero')
      haz1[!ss,] = exp(Xb1) %*% t(lam1[[1]])
      Xb1 = predict(fit1[[2]], data.frame(M=m_b, X)[ss,], type='lp', reference='zero')
      haz1[ss,] = exp(Xb1) %*% t(lam1[[2]])
      if (is.null(h2)) {
        Xb2 = predict(fit2[[1]], data.frame(M=m_b, X)[!ss,], type='lp', reference='zero')
        haz2[!ss,] = exp(Xb2) %*% t(lam2[[1]])
        Xb2 = predict(fit2[[2]], data.frame(M=m_b, X)[ss,], type='lp', reference='zero')
        haz2[ss,] = exp(Xb2) %*% t(lam2[[2]])
      } else {
        haz2 = t(apply(cbind(0,h2(tt,X,m_b)),1,diff))
      }
      term2 = term2 + exp(t(apply(-haz1 - haz2, 1, cumsum))) * haz1 * w_list[b]
    }
    Fm_reg[k,] = rowMeans(apply(term2, 1, cumsum))
    term = t(apply(term1 + term2, 1, cumsum))
    term[D==0,] = 0
    Eterm = matrix(0, nrow=n, ncol=mt)
    fit_phi = lm(term ~ ., data=fit$data[,-1], weights=wts, subset=(D==1 & ss))
    Eterm[!ss,] = predict(fit_phi, newdata=fit$data[!ss,-1])
    fit_phi = lm(term ~ ., data=fit$data[,-1], weights=wts, subset=(D==1 & !ss))
    Eterm[ss,] = predict(fit_phi, newdata=fit$data[ss,-1])
    term = wts*term + (1-wts)*Eterm
    Fm_eff[k,] = pmin(1,pmax(0,colMeans(term, na.rm=TRUE)))
    if (monotone){
      Fm_eff[k,] = isoreg(Fm_eff[k,])$yf
    } 
    se[k,] = sqrt(apply(term^2, 2, mean, na.rm=TRUE)) / sqrt(n-1)
  }
  colnames(Fm_reg) = colnames(Fm_eff) = colnames(se) = round(tt, 4)
  rownames(Fm_reg) = rownames(Fm_eff) = rownames(se) = round(m_seq, 4)
  if (!is.null(time)){
    idx = max(which(tt <= time))
    Fm_reg = Fm_reg[, idx]; Fm_eff = Fm_eff[, idx]; se = se[, idx]
  }
  return(list(time=tt, a=a, m_seq=m_seq, Fm_reg=Fm_reg, Fm_eff=Fm_eff, se=se))
}
