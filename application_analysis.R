## Application
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
w_list = quad$weights/sum(quad$weights)

riesz_basis <- function(M, X, keep = NULL) {
  X <- as.matrix(X)
  # An unadjusted analysis supplies a zero-column X; the basis is then (1, M, M^2).
  # Handle it separately because paste0(NULL, '^2') would yield a spurious name.
  if (ncol(X) == 0L) {
    nm <- c('(Intercept)', 'M', 'M2')
    B <- if (length(M) == 0L) matrix(numeric(0), 0L, 3L, dimnames = list(NULL, nm))
         else { b <- cbind(1, M, M^2); colnames(b) <- nm; b }
    if (!is.null(keep)) B <- B[, keep, drop = FALSE]
    return(B)
  }
  if (is.null(colnames(X))) colnames(X) <- paste0('X', seq_len(ncol(X)))
  xn <- colnames(X)
  nm <- c('(Intercept)', 'M', 'M2', xn, paste0('M:', xn), paste0(xn, '^2'))
  if (length(M) == 0L) {
    # cbind() would recycle the scalar intercept against zero-length columns
    B <- matrix(numeric(0), nrow = 0L, ncol = length(nm), dimnames = list(NULL, nm))
  } else {
    B <- cbind(1, M, M^2, X, M * X, X^2)
    colnames(B) <- nm
  }
  if (!is.null(keep)) B <- B[, keep, drop = FALSE]
  B
}

# Maximal set of linearly independent basis columns via QR with pivoting.
# A binary covariate makes X^2 identical to X, so the raw basis is rank deficient;
# dqrdc2 keeps the earlier columns and moves the dependent ones to the end.
riesz_keep <- function(Phi, tol = 1e-7) {
  q <- qr(Phi, tol = tol)
  sort(q$pivot[seq_len(q$rank)])
}

fit_riesz_r10 <- function(M, X, A, D, wts, ps, train, lambda = 0.01) {
  X <- as.matrix(X)
  D <- as.integer(D)
  keep1 <- train & D == 1 & A == 1
  keep0 <- train & D == 1 & A == 0
  Phi1 <- riesz_basis(M[keep1], X[keep1, , drop = FALSE])
  Phi0 <- riesz_basis(M[keep0], X[keep0, , drop = FALSE])
  # Fit and predict must share one basis subset, so pick it on the pooled training design.
  cols <- riesz_keep(rbind(Phi1, Phi0))
  Phi1 <- Phi1[, cols, drop = FALSE]
  Phi0 <- Phi0[, cols, drop = FALSE]
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
    keep = cols,
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
  rhat <- drop(riesz_basis(M_new, X_new, keep = fit$keep) %*% fit$theta)
  rhat <- pmin(pmax(rhat, lower), upper)
  rhat
}

validate_cf_split <- function(ss, n) {
  if (is.null(ss)) {
    ss <- seq_len(n) %in% sample.int(n, floor(n / 2))
  }
  if (!is.logical(ss) || length(ss) != n || anyNA(ss) ||
      !any(ss) || all(ss)) {
    stop("ss must be a logical vector with both folds present and no missing values.")
  }
  ss
}

## The sampling weight plays two separate roles: it restores the population in the
## nuisance regressions, and it scales the augmentation psi = Ehat + wts*(phi - Ehat).
## wts_cap bounds both at once and therefore corrupts the nuisance fits; aug_cap
## bounds only the multiplier in the augmentation step, leaving the regressions with
## the full weights. Shrinking aug_cap pulls the estimate toward the plain regression,
## which is a bias-variance trade-off rather than a correction.
## sample_df supplies the covariates of the two-phase sampling model, defaulting to the
## analysis (Time, status). When the outcome has been administratively truncated, pass
## the untruncated follow-up and case indicator: sampling was driven by case status over
## the full follow-up, and once that variable is recoded the weights stop restoring the
## cohort (here truncation at 365 days inflated the CIF by a factor of ~4.7).
## wts overrides the estimated two-phase weight. Supplying a constant vector drops the
## inverse-sampling correction and treats the measured subcohort as if it were a simple
## random sample; it is still multiplied by D, since a subject without a mediator cannot
## enter any model that conditions on M.
## target selects the population the counterfactual risk is averaged over; NULL uses
## everyone. Restricting it to the treated arm avoids asking the augmentation regression,
## which is fit only on mediator-measured (hence treated) subjects, to extrapolate to an
## arm where A has no support. Randomisation makes the two targets agree in expectation.
## aug_time = FALSE drops Time from the augmentation regression E[phi | A, X, Time, status];
## the choice is stored on the fit and inherited by est_con. The sampling-weight model
## is unaffected.
fit_eif_cf <- function(A, Time, status, M, X, method='riesz', event=1,
                       monotone=TRUE, ss=NULL, wts_cap=100, aug_cap=Inf,
                       sample_df=NULL, wts=NULL, target=NULL, aug_time=TRUE,
                       a_seq=rbind(c(0,0,0), c(1,0,0), c(1,1,0), c(1,1,1))){
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
  ss = validate_cf_split(ss, n)
  
  ## estimate sampling weight
  if (is.null(wts)) {
    wts = rep(1,n)
    if (sum(D)<n){
      sdf = if (is.null(sample_df)) data.frame(A, X, Time, st=as.factor(status)) else
            data.frame(A, X, sample_df)
      iwts = glm(D~., data=sdf, family='binomial', subset=ss)
      wts[!ss] = D[!ss]/predict(iwts, newdata=sdf[!ss,,drop=FALSE], type='response')
      iwts = glm(D~., data=sdf, family='binomial', subset=!ss)
      wts[ss] = D[ss]/predict(iwts, newdata=sdf[ss,,drop=FALSE], type='response')
    }
    wts = pmin(wts_cap, pmax(wts, 0.01))
  } else {
    if (!is.numeric(wts) || length(wts) != n || any(!is.finite(wts)) || any(wts < 0))
      stop('wts must be a finite non-negative vector of length n.')
  }
  wts = wts * D
  
  ## A cause-specific model is unavailable when its training subset is empty, e.g. an
  ## arm in which the mediator was never measured. Its hazard is then taken as zero,
  ## which is harmless only for targets that give that arm zero weight; arm_ok below
  ## records which counterfactual arms remain estimable.
  safe_cox <- function(sub, ev) {
    if (!any(sub)) return(NULL)
    coxph(Surv(Time, ev) ~ ., data = data.frame(M, X), weights = wts, subset = sub)
  }
  safe_lam <- function(fo) {
    if (is.null(fo)) return(rep(0, mt))
    bh <- rbind(c(0, 0), basehaz(fo, centered = FALSE))
    diff(c(0, matchy(bh$time, bh$hazard, tt)))
  }
  safe_haz <- function(fo, lam, idx) {
    if (is.null(fo)) return(matrix(0, sum(idx), mt))
    exp(predict(fo, data.frame(M, X)[idx, , drop = FALSE],
                type = 'lp', reference = 'zero')) %*% t(lam)
  }


  ## fit nuisance models
  fit11_cf = fit01_cf = fit12_cf = fit02_cf = fit10_cf = fit00_cf = list(NA,NA)
  lam11_cf = lam01_cf = lam12_cf = lam02_cf = lam10_cf = lam00_cf = list(NA,NA)
  # fit model on ss, predict on 1ss
  {
    fit11 = safe_cox(A==1&D==1&ss, status==1)
    fit01 = safe_cox(A==0&D==1&ss, status==1)
    fit12 = safe_cox(A==1&D==1&ss, status==2)
    fit02 = safe_cox(A==0&D==1&ss, status==2)
    fit10 = safe_cox(A==1&D==1&ss, status==0)
    fit00 = safe_cox(A==0&D==1&ss, status==0)
    lam11 = safe_lam(fit11); lam01 = safe_lam(fit01); lam12 = safe_lam(fit12)
    lam02 = safe_lam(fit02); lam10 = safe_lam(fit10); lam00 = safe_lam(fit00)
    haz11_1ss = safe_haz(fit11, lam11, !ss)
    haz01_1ss = safe_haz(fit01, lam01, !ss)
    haz12_1ss = safe_haz(fit12, lam12, !ss)
    haz02_1ss = safe_haz(fit02, lam02, !ss)
    haz10_1ss = safe_haz(fit10, lam10, !ss)
    haz00_1ss = safe_haz(fit00, lam00, !ss)
    # single-bracket assignment stores a NULL model instead of deleting the slot
    fit11_cf[1] = list(fit11); fit12_cf[1] = list(fit12); fit10_cf[1] = list(fit10)
    fit01_cf[1] = list(fit01); fit02_cf[1] = list(fit02); fit00_cf[1] = list(fit00)
    lam11_cf[[1]] = lam11; lam12_cf[[1]] = lam12; lam10_cf[[1]] = lam10
    lam01_cf[[1]] = lam01; lam02_cf[[1]] = lam02; lam00_cf[[1]] = lam00
  }
  # fit model on 1ss, predict on ss
  {
    fit11 = safe_cox(A==1&D==1&!ss, status==1)
    fit01 = safe_cox(A==0&D==1&!ss, status==1)
    fit12 = safe_cox(A==1&D==1&!ss, status==2)
    fit02 = safe_cox(A==0&D==1&!ss, status==2)
    fit10 = safe_cox(A==1&D==1&!ss, status==0)
    fit00 = safe_cox(A==0&D==1&!ss, status==0)
    lam11 = safe_lam(fit11); lam01 = safe_lam(fit01); lam12 = safe_lam(fit12)
    lam02 = safe_lam(fit02); lam10 = safe_lam(fit10); lam00 = safe_lam(fit00)
    haz11_ss = safe_haz(fit11, lam11, ss)
    haz01_ss = safe_haz(fit01, lam01, ss)
    haz12_ss = safe_haz(fit12, lam12, ss)
    haz02_ss = safe_haz(fit02, lam02, ss)
    haz10_ss = safe_haz(fit10, lam10, ss)
    haz00_ss = safe_haz(fit00, lam00, ss)
    fit11_cf[2] = list(fit11); fit12_cf[2] = list(fit12); fit10_cf[2] = list(fit10)
    fit01_cf[2] = list(fit01); fit02_cf[2] = list(fit02); fit00_cf[2] = list(fit00)
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
  ps_a[!ss] = predict(fit_a, newdata=data.frame(A,X)[!ss,,drop=FALSE], type='response')
  fit_a = glm(A~., data=data.frame(A,X), family='binomial', subset=!ss)
  ps_a[ss] = predict(fit_a, newdata=data.frame(A,X)[ss,,drop=FALSE], type='response')
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
    mediator_logvar = mean(residuals(sx)^2)
    sx1 = exp(predict(sx, newdata=data.frame(A=1,X,AX=X)[!ss,]) + mediator_logvar/2)
    sx0 = exp(predict(sx, newdata=data.frame(A=0,X,AX=0*X)[!ss,]) + mediator_logvar/2)
    mx_1ss = cbind(mx0, mx1)
    sx_1ss = sqrt(cbind(sx0, sx1))
    mx = lm(M ~ ., data=data.frame(A,X,AX=A*X), weights=wts, subset=(D==1&!ss))
    res_m = rep(NA, n)
    res_m[D==1&!ss] = residuals(mx)^2
    mx1 = predict(mx, newdata=data.frame(A=1,X,AX=X)[ss,])
    mx0 = predict(mx, newdata=data.frame(A=0,X,AX=0*X)[ss,])
    sx = lm(log(res_m+0.000001) ~ ., data=data.frame(A,X,AX=A*X), weights=wts, subset=(D==1&!ss))
    mediator_logvar = mean(residuals(sx)^2)
    sx1 = exp(predict(sx, newdata=data.frame(A=1,X,AX=X)[ss,]) + mediator_logvar/2)
    sx0 = exp(predict(sx, newdata=data.frame(A=0,X,AX=0*X)[ss,]) + mediator_logvar/2)
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
  # a = (a0, a1, a2); a_seq lists the counterfactual arms to estimate, one per row.
  # (0,0,0) and (1,1,1) must be included: est_con uses their phi terms.
  n = length(A)
  a_seq = as.matrix(a_seq)
  n_arm = nrow(a_seq)
  F_reg = F_eff = se_eff = matrix(0, nrow=n_arm, ncol=length(tt))
  a_names = apply(a_seq, 1, paste0, collapse="")
  rownames(F_reg) = rownames(F_eff) = rownames(se_eff) = a_names
  eif = array(0, dim=c(n,length(tt),n_arm), dimnames=list(NULL,round(tt,4),a_names))
  tg = if (is.null(target)) rep(TRUE, n) else target
  if (!is.logical(tg) || length(tg) != n || anyNA(tg) || !any(tg))
    stop('target must be a logical vector of length n with at least one TRUE.')
  ntg = sum(tg)
  arm_ok = rep(TRUE, n_arm)
  for (k in 1:n_arm){
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
    arm_ok[k] = !is.null(fit1[[1]]) && !is.null(fit1[[2]]) &&
                !is.null(fit2[[1]]) && !is.null(fit2[[2]]) &&
                any(idx0_ss) && any(idx0_1ss)
    gamma_obs = matrix(NA, nrow=n, ncol=mt)
    h_hat = matrix(0, nrow=n, ncol=mt)
    Xdf = data.frame(X)
    # With no covariates the outcome regression reduces to an intercept
    lm_gamma = function(sub) if (ncol(Xdf) == 0L)
      lm(gamma_obs ~ 1, weights=wts, subset=sub) else
      lm(gamma_obs ~ ., data=Xdf, weights=wts, subset=sub)
    if (arm_ok[k]) {
      h1o = safe_haz(fit1[[1]], lam1[[1]], idx0_ss)
      h2o = safe_haz(fit2[[1]], lam2[[1]], idx0_ss)
      gamma_obs[idx0_ss,] = exp(t(apply(-h1o-h2o, 1, cumsum))) * h1o
      fit_h = lm_gamma(idx0_ss)
      h_hat[!ss,] = predict(fit_h, newdata=Xdf[!ss,,drop=FALSE])
      h1o = safe_haz(fit1[[2]], lam1[[2]], idx0_1ss)
      h2o = safe_haz(fit2[[2]], lam2[[2]], idx0_1ss)
      gamma_obs[idx0_1ss,] = exp(t(apply(-h1o-h2o, 1, cumsum))) * h1o
      fit_h = lm_gamma(idx0_1ss)
      h_hat[ss,] = predict(fit_h, newdata=Xdf[ss,,drop=FALSE])
      h_hat[h_hat<0] = 0
    }
    dF = h_hat
    F_reg[k,] = rowMeans(apply(dF[tg,,drop=FALSE], 1, cumsum), na.rm=TRUE)
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
    aug_df = if (aug_time) data.frame(A,X,Time,st=as.factor(status)) else
             data.frame(A,X,st=as.factor(status))
    fit_phi = lm(Fx ~ ., data=aug_df, weights=wts, subset=(D==1&ss))
    EFx[!ss,] = predict(fit_phi, newdata=aug_df[!ss,,drop=FALSE])
    fit_phi = lm(Fx ~ ., data=aug_df, weights=wts, subset=(D==1&!ss))
    EFx[ss,] = predict(fit_phi, newdata=aug_df[ss,,drop=FALSE])
    Fx = EFx + pmin(wts, aug_cap) * (Fx - EFx)
    F_effk = pmin(1,pmax(0,colMeans(Fx[tg,,drop=FALSE], na.rm=TRUE)))
    if (monotone){
      F_effk = isoreg(F_effk)$yf
    } 
    eifk = t(t(Fx) - F_effk)
    se_effk = sqrt(apply(eifk[tg,,drop=FALSE]^2, 2, mean, na.rm=TRUE)) / sqrt(ntg-1)
    F_eff[k,] = F_effk
    se_eff[k,] = se_effk
    eif[,,k] = eifk
  }
  names(arm_ok) = a_names
  if (!all(arm_ok)) {
    warning('Arms ', paste(a_names[!arm_ok], collapse=', '),
            ' are not estimable with the available mediator data; set to NA.')
    F_reg[!arm_ok,] = NA; F_eff[!arm_ok,] = NA; se_eff[!arm_ok,] = NA
    eif[,,!arm_ok] = NA
  }
  colnames(F_reg) = colnames(F_eff) = round(tt,4)
  
  return(list(data=dat, time=tt, phi0=phi0, phi1=phi1, phi00=phi00, phi10=phi10, 
              eif=eif, wts=wts, ps=ps,
              lam01=lam01_cf, lam02=lam02_cf, fit01=fit01_cf, fit02=fit02_cf,
              lam11=lam11_cf, lam12=lam12_cf, fit11=fit11_cf, fit12=fit12_cf,
              F_reg=F_reg, F_eff=F_eff, se=se_eff, cf=TRUE, ss=ss, arm_ok=arm_ok,
              aug_time=aug_time))
}

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
              TE_reg=TE_reg, TE_eff=TE_eff, se_TE=se_TE, RMSTL_TE=RMSTL_TE, RMSTL_TE_se=RMSTL_TE_se, p_TE=p_TE,
              ME_reg=ME_reg, ME_eff=ME_eff, se_ME=se_ME, RMSTL_ME=RMSTL_ME, RMSTL_ME_se=RMSTL_ME_se, p_ME=p_ME,
              DE_reg=DE_reg, DE_eff=DE_eff, se_DE=se_DE, RMSTL_DE=RMSTL_DE, RMSTL_DE_se=RMSTL_DE_se, p_DE=p_DE,
              CE_reg=CE_reg, CE_eff=CE_eff, se_CE=se_CE, RMSTL_CE=RMSTL_CE, RMSTL_CE_se=RMSTL_CE_se, p_CE=p_CE))
}

## Effects along an arbitrary sequence of arms, path[1] -> path[2] -> ... -> path[K]:
## effect j is F[path[j+1]] - F[path[j]], so the effects sum to the total effect
## F[path[K]] - F[path[1]]. names(path)[-1] name the effects; the result uses the field
## names of est_decompose (TE_eff, se_TE, RMSTL_TE, p_TE, ...) so plot_effects can use it.
est_path_effects <- function(fit, path) {
  n = nrow(fit$data)
  ctt = diff(c(0, fit$time))
  out = list(time = fit$time, F_eff = fit$F_eff)
  one = function(to, from, nm) {
    est = fit$F_eff[to,] - fit$F_eff[from,]
    eif_d = fit$eif[,,to] - fit$eif[,,from]
    rm_se = sd(colSums(t(eif_d)*ctt, na.rm=TRUE)) / sqrt(n)
    rm = sum(est*ctt, na.rm=TRUE)
    setNames(list(est, apply(eif_d, 2, sd, na.rm=TRUE) / sqrt(n), rm, rm_se,
                  2*pnorm(-abs(rm/rm_se))),
             paste0(c('', 'se_', 'RMSTL_', 'RMSTL_', 'p_'), nm,
                    c('_eff', '', '', '_se', '')))
  }
  K = length(path)
  out = c(out, one(path[K], path[1], 'TE'))
  for (j in 2:K) out = c(out, one(path[j], path[j-1], names(path)[j]))
  out
}

GK <- function(u,h=1) {
  (1/sqrt(2*pi)) * exp(-((u/h)^2)/2) / h
}

fit_riesz_controlled <- function(m0, h, M, X, A, D, wts, ps, a, train, lambda = 0.01) {
  X <- as.matrix(X)
  D <- as.integer(D)
  keep_a <- train & D == 1 & A == a
  Phi_a <- riesz_basis(M[keep_a], X[keep_a,,drop = FALSE])
  cols <- riesz_keep(Phi_a)
  Phi_a <- Phi_a[, cols, drop = FALSE]
  w_a <- wts[keep_a] / ps[keep_a]
  w_a <- w_a / sum(w_a)
  H <- crossprod(Phi_a, Phi_a * w_a)
  g <- rep(0, ncol(Phi_a))
  X_a <- X[keep_a, , drop = FALSE]
  for (b in 1:ngrid) {
    m_b <- m0 + sqrt(2) * h * b_list[b]
    Phi_b <- riesz_basis(rep(m_b, nrow(X_a)), X_a, keep = cols)
    g <- g + w_list[b] * colSums(Phi_b * w_a)
  }
  d <- ncol(Phi_a)
  theta <- solve(H + lambda * diag(d), g)
  list(
    theta = theta,
    keep = cols,
    lambda = lambda,
    a = a,
    m0 = m0,
    h = h,
    n_train = sum(keep_a),
    ess_train = sum(w_a)^2 / sum(w_a^2)
  )
}

predict_riesz_controlled <- function(fit, M_new, X_new, lower = 1e-5, upper = 99) {
  alpha_hat <- drop(riesz_basis(M_new, X_new, keep = fit$keep) %*% fit$theta)
  alpha_hat <- pmin(pmax(alpha_hat, lower), upper)
  alpha_hat
}

## aug_strat = TRUE fits the augmentation regression separately for cases and
## non-cases, on X only. With the default additive model, status can only shift the
## intercept, so the Time slope is learned from cases alone (sampled non-cases barely
## vary in Time); unsampled non-cases who were censored early then inherit the case
## slope, and the estimated CIF turns into an inverted U.
est_con <- function(fit, a=0, m_seq=NULL, h2=NULL, time=NULL, h=0.5, monotone=TRUE,
                    aug_cap=Inf, target=NULL, aug_strat=FALSE,
                    aug_time=if (is.null(fit$aug_time)) TRUE else fit$aug_time){
  n = nrow(fit$data)
  tg = if (is.null(target)) rep(TRUE, n) else target
  if (!is.logical(tg) || length(tg) != n || anyNA(tg) || !any(tg))
    stop('target must be a logical vector of length n with at least one TRUE.')
  ntg = sum(tg)
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
  X = as.matrix(fit$data[, grepl("^X", names(fit$data)) | 
                           names(fit$data) == "X" |
                           colnames(fit$data)%in%Xname, drop=FALSE])
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
    Xb1 = predict(fit1[[1]], data.frame(M,X)[!ss,,drop=FALSE], type='lp', reference='zero')
    haz1[!ss,] = exp(Xb1) %*% t(lam1[[1]])
    Xb1 = predict(fit1[[2]], data.frame(M,X)[ss,,drop=FALSE], type='lp', reference='zero')
    haz1[ss,] = exp(Xb1) %*% t(lam1[[2]])
    if (is.null(h2)){
      Xb2 = predict(fit2[[1]], data.frame(M,X)[!ss,,drop=FALSE], type='lp', reference='zero')
      haz2[!ss,] = exp(Xb2) %*% t(lam2[[1]])
      Xb2 = predict(fit2[[2]], data.frame(M,X)[ss,,drop=FALSE], type='lp', reference='zero')
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
      Xb1 = predict(fit1[[1]], data.frame(M=m_b, X)[!ss,,drop=FALSE], type='lp', reference='zero')
      haz1[!ss,] = exp(Xb1) %*% t(lam1[[1]])
      Xb1 = predict(fit1[[2]], data.frame(M=m_b, X)[ss,,drop=FALSE], type='lp', reference='zero')
      haz1[ss,] = exp(Xb1) %*% t(lam1[[2]])
      if (is.null(h2)) {
        Xb2 = predict(fit2[[1]], data.frame(M=m_b, X)[!ss,,drop=FALSE], type='lp', reference='zero')
        haz2[!ss,] = exp(Xb2) %*% t(lam2[[1]])
        Xb2 = predict(fit2[[2]], data.frame(M=m_b, X)[ss,,drop=FALSE], type='lp', reference='zero')
        haz2[ss,] = exp(Xb2) %*% t(lam2[[2]])
      } else {
        haz2 = t(apply(cbind(0,h2(tt,X,m_b)),1,diff))
      }
      term2 = term2 + exp(t(apply(-haz1 - haz2, 1, cumsum))) * haz1 * w_list[b]
    }
    Fm_reg[k,] = rowMeans(apply(term2[tg,,drop=FALSE], 1, cumsum))
    term = t(apply(term1 + term2, 1, cumsum))
    term[D==0,] = 0
    Eterm = matrix(0, nrow=n, ncol=mt)
    if (aug_strat) {
      Xa = data.frame(X)
      is_case = fit$data$status > 0
      for (fold in list(ss, !ss)) for (cs in c(FALSE, TRUE)) {
        tr_i = D & fold & (is_case == cs)
        te_i = !fold & (is_case == cs)
        if (!any(te_i) || !any(tr_i)) next
        fit_phi = if (ncol(Xa) == 0L) lm(term ~ 1, weights=wts, subset=tr_i) else
          lm(term ~ ., data=Xa, weights=wts, subset=tr_i)
        Eterm[te_i,] = predict(fit_phi, newdata=Xa[te_i,,drop=FALSE])
      }
    } else {
      aug_df = fit$data[, setdiff(names(fit$data), c('M', if (!aug_time) 'Time')), drop=FALSE]
      # status as a factor, matching the augmentation regression in fit_eif_cf
      aug_df$status = as.factor(aug_df$status)
      fit_phi = lm(term ~ ., data=aug_df, weights=wts, subset=(D==1 & ss))
      Eterm[!ss,] = predict(fit_phi, newdata=aug_df[!ss,,drop=FALSE])
      fit_phi = lm(term ~ ., data=aug_df, weights=wts, subset=(D==1 & !ss))
      Eterm[ss,] = predict(fit_phi, newdata=aug_df[ss,,drop=FALSE])
    }
    term = Eterm + pmin(wts, aug_cap) * (term - Eterm)
    Fm_eff[k,] = pmin(1,pmax(0,colMeans(term[tg,,drop=FALSE], na.rm=TRUE)))
    if (monotone){
      Fm_eff[k,] = isoreg(Fm_eff[k,])$yf
    } 
    se[k,] = sqrt(apply(term[tg,,drop=FALSE]^2, 2, mean, na.rm=TRUE)) / sqrt(ntg-1)
  }
  colnames(Fm_reg) = colnames(Fm_eff) = colnames(se) = round(tt, 4)
  rownames(Fm_reg) = rownames(Fm_eff) = rownames(se) = round(m_seq, 4)
  if (!is.null(time)){
    idx = max(which(tt <= time))
    Fm_reg = Fm_reg[, idx]; Fm_eff = Fm_eff[, idx]; se = se[, idx]
  }
  return(list(time=tt, a=a, m_seq=m_seq, Fm_reg=Fm_reg, Fm_eff=Fm_eff, se=se))
}


## Synthetic data
## Data restrictions rule out synthesizing COV2008 records directly. synth_cov2008.csv is
## generated by simu_cov2008.R from a simu0-style model calibrated to the 74 trial
## participants with an observed mediator: subsample Cox coefficients (so the trt/ctr
## hazard ratio is the observed one), constant baselines matched to the pooled day-315
## CIF, and the subsample covariate and mediator distributions.
## M is NA outside the two-phase sample, measured with probability 0.3 for cases and 0.2
## for non-cases.
library(dplyr)
library(ggplot2)

setwd('D:/Data/COVID_Vaccine')
synth = read.csv('synth_cov2008.csv')
Xname = c('Vacc','Sex','Age')
# Remark: the effect of trt on M is very weak. Vacc is highly correlated with M!

df <- rbind(
  data.frame(M = synth$M[synth$trt == 1], group = "Treated"),
  data.frame(M = synth$M[synth$ctr == 1], group = "Control")
)
p_box <- ggplot(df[!is.na(df$M), ], aes(x = group, y = M)) +
  labs(x = 'Group', y = expression('log'[10] * ' level'), 
       title = expression('Boxplot of log'[10] *
                            ' CD4+ T cells expressing IFN-'*gamma*' and/or IL-2')) +
  theme_bw() +
  geom_boxplot()
ggsave(file.path('app1', 'app1_boxplot_m.pdf'), p_box, width = 400 / 72, height = 247 / 72)

## Analysis
## Figures follow the paper and its supplement. For event k (1 = asymptomatic, 2 =
## symptomatic first infection; the other type is the competing event):
##   app1_cif{k}, app1_decompose{k}, app1_cif{k}_control, app1_cif{k}_control2
##     -> Fig. 3 (k = 1) / Fig. S10 (k = 2), panels (a)-(d)
##   app1_curve{k}_control, app1_curve{k}_control2 -> Fig. S9 / S12, panels (a)-(b)
##   app1_sensitivity{k} -> Fig. S8 / S11 (saved after the sensitivity analysis below)
## y-limits are taken from the data so that no confidence band is cut off.
event_names <- c('asymptomatic', 'symptomatic')
marker_lab <- expression('log'[10] * ' CD4+ T cells expressing IFN-'*gamma*' and/or IL-2')
cif_cols <- c('(0,0,0)' = 'black', '(1,0,0)' = 'brown',
              '(1,1,0)' = 'darkorange', '(1,1,1)' = 'lightblue')
# Marker grids stay within the 10%-90% quantiles of observed M in the treated arm
m_q <- quantile(synth$M[synth$trt == 1], c(0.1, 0.9), na.rm = TRUE)

# Smallest multiple of step that covers every value of x
ceil_to <- function(x, step) step * ceiling(max(x, na.rm = TRUE) / step - 1e-8)

save_app1 <- function(p, name, width = 5, height = 10 / 3) {
  ggsave(file.path('app1', paste0(name, '.pdf')), p, width = width, height = height)
}

# Monotone CIF curves with pointwise 95% bands clipped to [0, 1]; se at t = 0 is 0
cif_band <- function(Fmat, semat, labels, levels = labels) {
  time_vec <- c(0, as.numeric(colnames(Fmat)))
  Fm <- t(apply(cbind(0, Fmat), 1, function(f) isoreg(f)$yf))
  SE <- cbind(0, semat)
  data.frame(time = rep(time_vec, times = nrow(Fm)),
             CIF = as.vector(t(Fm)),
             lower = pmax(0, as.vector(t(Fm - 1.96 * SE))),
             upper = pmin(1, as.vector(t(Fm + 1.96 * SE))),
             group = factor(rep(labels, each = length(time_vec)), levels = levels))
}

plot_cif <- function(df, cols, legend_name, panel, legend_pos) {
  ggplot(df, aes(x = time, y = CIF, colour = group)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 0.7) +
    scale_colour_manual(values = cols, name = legend_name,
                        guide = guide_legend(position = 'inside')) +
    scale_fill_manual(values = cols, guide = 'none') +
    coord_cartesian(ylim = c(0, ceil_to(df$upper, 0.1))) +
    labs(x = 'Study time (days)', y = 'Cumulative incidence', title = panel) +
    theme_bw() +
    theme(legend.position.inside = legend_pos,
          legend.background = element_rect(fill = alpha('white', 0.6)))
}

plot_effects <- function(natural, panel) {
  labs3 <- c('Mediation effect', 'Direct effect', 'Competing effect')
  est <- rbind(natural$ME_eff, natural$DE_eff, natural$CE_eff)
  se <- rbind(natural$se_ME, natural$se_DE, natural$se_CE)
  time_vec <- c(0, natural$time)
  df <- do.call(rbind, lapply(1:3, function(i) data.frame(
    time = time_vec, value = c(0, est[i, ]),
    lower = c(0, est[i, ] - 1.96 * se[i, ]), upper = c(0, est[i, ] + 1.96 * se[i, ]),
    effect = factor(labs3[i], levels = labs3))))
  # common (asymmetric) range across the three panels, rounded out to 0.02
  ylim <- c(-ceil_to(-c(df$lower, 0), 0.02), ceil_to(c(df$upper, 0), 0.02))
  ggplot(df, aes(x = time, y = value)) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line(linewidth = 0.6) +
    scale_y_continuous(breaks = seq(ylim[1], ylim[2], by = 0.02)) +
    coord_cartesian(ylim = ylim) +
    facet_wrap(~ effect, nrow = 1) +
    labs(x = 'Study time (days)', y = 'Treatment effect', title = panel) +
    theme_bw() +
    theme(strip.text = element_text(face = 'bold'))
}

# Day-tpoint controlled CIF against the marker; log10 axis labelled on the 10^m scale
plot_marker <- function(control, panel, tpoint = 315) {
  tp <- as.character(max(control$time[control$time <= tpoint]))
  df <- data.frame(marker = 10^control$m_seq, value = control$Fm_eff[, tp],
                   lower = pmax(0, control$Fm_eff[, tp] - 1.96 * control$se[, tp]),
                   upper = pmin(1, control$Fm_eff[, tp] + 1.96 * control$se[, tp]))
  ggplot(df, aes(x = marker, y = value)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line(linewidth = 0.6) +
    scale_x_log10(breaks = df$marker, labels = formatC(df$marker, format = 'e', digits = 1)) +
    coord_cartesian(ylim = c(0, ceil_to(df$upper, 0.1))) +
    labs(x = marker_lab, y = 'Cumulative incidence', title = panel) +
    theme_bw()
}

# Rows of an est_con result at the requested grid points
subset_con <- function(con, m) {
  i <- match(round(m, 8), round(con$m_seq, 8))
  list(time = con$time, m_seq = con$m_seq[i], Fm_reg = con$Fm_reg[i, , drop = FALSE],
       Fm_eff = con$Fm_eff[i, , drop = FALSE], se = con$se[i, , drop = FALSE])
}

run_app1 <- function(event, seed = 2025) {
  set.seed(seed)
  fit <- fit_eif_cf(synth$trt, synth$Time, synth$status, synth$M, synth[, Xname],
                    method = 'riesz', event = event, aug_time = FALSE)
  natural <- est_decompose(fit)

  # Cumulative hazard of the competing event under active treatment, halved; Cox
  # coefficients and baseline averaged over the two folds (n x length(fit$time)).
  h2cox <- (coef(fit$fit12[[1]]) + coef(fit$fit12[[2]])) / 2
  Lam2 <- cumsum((fit$lam12[[1]] + fit$lam12[[2]]) / 2)
  h2 <- function(tt, x, m) {
    stopifnot(length(tt) == length(Lam2))
    lp <- drop(h2cox['M'] * m + as.matrix(x[, names(h2cox)[-1], drop = FALSE]) %*% h2cox[-1])
    exp(lp) %o% Lam2 / 2
  }
  m3 <- seq(m_q[1], m_q[2], length = 3)
  m4 <- seq(m_q[1], m_q[2], length = 4)
  m_all <- sort(unique(c(m3, m4)))
  con_nat <- est_con(fit, a = 1, m_seq = m_all, h = 0.2)
  con_half <- est_con(fit, a = 1, m_seq = m_all, h2 = h2, h = 0.2)

  m_groups <- paste('m =', formatC(10^m3, format = 'e', digits = 1))
  m_cols <- setNames(c('black', 'brown', 'darkorange'), m_groups)
  con_band <- function(con) {
    s <- subset_con(con, m3)
    cif_band(s$Fm_eff, s$se, m_groups, levels = rev(m_groups))
  }
  plots <- list(
    cif = plot_cif(cif_band(fit$F_eff, fit$se, names(cif_cols)), cif_cols,
                   expression('(' * a[0] * ',' * a[1] * ',' * a[2] * ')'), '(a)', c(0.15, 0.7)),
    decompose = plot_effects(natural, '(b)'),
    cif_control = plot_cif(con_band(con_nat), m_cols, 'Marker', '(c)', c(0.15, 0.75)),
    cif_control2 = plot_cif(con_band(con_half), m_cols, 'Marker', '(d)', c(0.15, 0.75)),
    curve_control = plot_marker(subset_con(con_nat, m4), '(a)'),
    curve_control2 = plot_marker(subset_con(con_half, m4), '(b)')
  )
  list(event = event, fit = fit, natural = natural, con_nat = con_nat, con_half = con_half,
       p.vals = c(TE = natural$p_TE, ME = natural$p_ME, DE = natural$p_DE, CE = natural$p_CE),
       RMSTL = c(TE = natural$RMSTL_TE, ME = natural$RMSTL_ME,
                 DE = natural$RMSTL_DE, CE = natural$RMSTL_CE),
       plots = plots)
}

app1_res <- lapply(1:2, run_app1)
for (k in 1:2) {
  pl <- app1_res[[k]]$plots
  save_app1(pl$cif, paste0('app1_cif', k))
  save_app1(pl$decompose, paste0('app1_decompose', k))
  save_app1(pl$cif_control, paste0('app1_cif', k, '_control'))
  save_app1(pl$cif_control2, paste0('app1_cif', k, '_control2'))
  save_app1(pl$curve_control, paste0('app1_curve', k, '_control'))
  save_app1(pl$curve_control2, paste0('app1_curve', k, '_control2'))
}
sapply(app1_res, function(r) round(rbind(RMSTL = r$RMSTL, p = r$p.vals), 3), simplify = FALSE)

## Sensitivity analysis: alternative ordering of the decomposition. The total effect
## F^(111) - F^(000) is split along (000) -> (010) -> (011) -> (111):
##   direct effect     F^(010) - F^(000)  (primary-event hazard switched to treatment)
##   competing effect  F^(011) - F^(010)  (then the competing-event hazard)
##   mediation effect  F^(111) - F^(011)  (then the mediator distribution)
## The main analysis switches the mediator first, (000) -> (100) -> (110) -> (111). The
## same seed gives the same cross-fitting split, so (000), (111) and the total effect
## equal the main analysis. Saved as app1_cif{k}_alt (panel a) and app1_decompose{k}_alt
## (panel b).
## A second ordering, (000) -> (100) -> (101) -> (111), switches the mediator, then the
## competing-event hazard, then the primary-event hazard:
##   mediation effect  F^(100) - F^(000)  (identical to the main analysis)
##   competing effect  F^(101) - F^(100)
##   direct effect     F^(111) - F^(101)
## Saved as app1_cif{k}_alt2 and app1_decompose{k}_alt2.
## path: arm labels in order; names(path)[-1] name the effect of each step (DE, CE, ME).
run_app1_alt <- function(event, path, seed = 2025) {
  a_seq <- do.call(rbind, lapply(strsplit(path, ''), as.numeric))
  labs <- paste0('(', sapply(strsplit(path, ''), paste, collapse = ','), ')')
  cols <- setNames(c('black', 'brown', 'darkorange', 'lightblue'), labs)
  set.seed(seed)
  fit <- fit_eif_cf(synth$trt, synth$Time, synth$status, synth$M, synth[, Xname],
                    method = 'riesz', event = event, aug_time = FALSE, a_seq = a_seq)
  eff <- est_path_effects(fit, path)
  ord <- c('TE', names(path)[-1])
  list(event = event, path = path, fit = fit, effects = eff,
       RMSTL = sapply(ord, function(e) eff[[paste0('RMSTL_', e)]]),
       p.vals = sapply(ord, function(e) eff[[paste0('p_', e)]]),
       plots = list(
         cif = plot_cif(cif_band(fit$F_eff, fit$se, labs), cols,
                        expression('(' * a[0] * ',' * a[1] * ',' * a[2] * ')'), '(a)',
                        c(0.15, 0.7)),
         decompose = plot_effects(eff, '(b)')))
}
alt_paths <- list(alt  = c('000', DE = '010', CE = '011', ME = '111'),
                  alt2 = c('000', ME = '100', CE = '101', DE = '111'))
app1_alt <- lapply(alt_paths, function(path) lapply(1:2, run_app1_alt, path = path))
for (nm in names(alt_paths)) for (k in 1:2) {
  save_app1(app1_alt[[nm]][[k]]$plots$cif, paste0('app1_cif', k, '_', nm))
  save_app1(app1_alt[[nm]][[k]]$plots$decompose, paste0('app1_decompose', k, '_', nm))
}
lapply(app1_alt, function(a) sapply(a, function(r) round(rbind(delay = -r$RMSTL, p = r$p.vals), 3),
                                    simplify = FALSE))

## Interaction between the primary- and competing-event pathways under treated mediator:
##   I(t) = F^(111) - F^(110) - F^(101) + F^(100)
## It is the difference between the competing effects of the main order, F^(111) - F^(110),
## and of the second alternative order, F^(101) - F^(100) (equivalently between their
## direct effects), so I = 0 means the main/alt2 split does not depend on the order.
## SE from the influence functions of the four arms, all from one fit (same split).
## Saved as app1_interaction{k}: (a) the four CIFs, (b) I(t) with 95% band.
est_contrast <- function(fit, w) {
  n <- nrow(fit$data); ctt <- diff(c(0, fit$time))
  est <- drop(crossprod(w, fit$F_eff[names(w), , drop = FALSE]))
  eif_c <- Reduce(`+`, lapply(names(w), function(a) w[[a]] * fit$eif[, , a]))
  rm <- sum(est * ctt)
  rm_se <- sd(colSums(t(eif_c) * ctt)) / sqrt(n)
  list(time = fit$time, est = est, se = apply(eif_c, 2, sd) / sqrt(n),
       RMSTL = rm, RMSTL_se = rm_se, p = 2 * pnorm(-abs(rm / rm_se)))
}
int_w <- c('111' = 1, '110' = -1, '101' = -1, '100' = 1)
int_cols <- c('(1,0,0)' = 'black', '(1,1,0)' = 'brown', '(1,0,1)' = 'darkorange',
              '(1,1,1)' = 'lightblue')
run_app1_int <- function(event, seed = 2025) {
  set.seed(seed)
  # (000) is needed internally by fit_eif_cf; it is not plotted
  fit <- fit_eif_cf(synth$trt, synth$Time, synth$status, synth$M, synth[, Xname],
                    method = 'riesz', event = event, aug_time = FALSE,
                    a_seq = rbind(c(0,0,0), c(1,0,0), c(1,1,0), c(1,0,1), c(1,1,1)))
  ic <- est_contrast(fit, int_w)
  arms <- c('100', '110', '101', '111')
  p_a <- plot_cif(cif_band(fit$F_eff[arms, ], fit$se[arms, ], names(int_cols)), int_cols,
                  expression('(' * a[0] * ',' * a[1] * ',' * a[2] * ')'), '(a)', c(0.15, 0.7))
  df <- data.frame(time = c(0, ic$time), value = c(0, ic$est),
                   lower = c(0, ic$est - 1.96 * ic$se), upper = c(0, ic$est + 1.96 * ic$se))
  ylim <- c(-ceil_to(-c(df$lower, 0), 0.01), ceil_to(c(df$upper, 0), 0.01))
  p_b <- ggplot(df, aes(x = time, y = value)) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line(linewidth = 0.6) +
    coord_cartesian(ylim = ylim) +
    labs(x = 'Study time (days)', y = 'Interaction', title = '(b)') +
    theme_bw()
  list(event = event, fit = fit, interaction = ic, plots = list(cif = p_a, interaction = p_b))
}
app1_int <- lapply(1:2, run_app1_int)
for (k in 1:2) {
  save_app1(app1_int[[k]]$plots$cif, paste0('app1_interaction', k, '_cif'))
  save_app1(app1_int[[k]]$plots$interaction, paste0('app1_interaction', k))
}
sapply(app1_int, function(r) with(r$interaction, round(c(
  I_end = tail(est, 1), se_end = tail(se, 1), delay = -RMSTL, delay_se = RMSTL_se, p = p), 4)))

## Sensitivity analysis with a log-normal frailty
# |s| is the frailty SD. s > 0 shares the frailty across both causes; s < 0 lets the
# same latent draw raise the cause-1 hazard and lower the cause-2 hazard.
# Retains the existing prior-quadrature correction and baseline-hazard approximation.
# These numerical checks do not establish an observed-data EIF under latent frailty.
fit_eif_cf_f <- function(A, Time, status, M, X, method='riesz', event=1,
                        monotone=TRUE, s=0.1, ss=NULL, n_quad=10L, aug_time=TRUE){
  method <- match.arg(method, c('riesz', 'ratio', 'propensity'))
  if (length(s) != 1L || !is.numeric(s) || !is.finite(s)) {
    stop("s must be a finite scalar; |s| is the frailty SD.")
  }
  if (length(n_quad) != 1L || !is.numeric(n_quad) || !is.finite(n_quad) ||
      n_quad < 1L || n_quad != floor(n_quad)) {
    stop("n_quad must be a positive integer.")
  }
  X <- as.matrix(X)
  n <- length(A)
  if (n < 4L || length(Time) != n || length(status) != n ||
      length(M) != n || nrow(X) != n || !is.numeric(X) ||
      anyNA(A) || !all(A %in% c(0, 1)) || any(!is.finite(Time)) ||
      any(Time < 0) || anyNA(status) || any(!is.finite(X)) ||
      any(!is.finite(M[!is.na(M)]))) {
    stop("Check input lengths, binary A, numeric complete X, and finite times/observed M.")
  }
  frailty_var <- s^2
  # The quadrature is over a frailty with SD |s|; the sign only sets its cause-2 loading.
  cause2_sign <- if (s < 0) -1 else 1
  if (s == 0) {
    frailty_nodes <- 0
    quad_weights <- 1
  } else {
    if (!requireNamespace('coxme', quietly = TRUE)) stop("Package 'coxme' is required for s != 0.")
    quadrature <- statmod::gauss.quad(n_quad, 'hermite')
    frailty_nodes <- sqrt(2) * abs(s) * quadrature$nodes
    quad_weights <- quadrature$weights / sum(quadrature$weights)
  }
  if (!is.finite(frailty_var) || any(!is.finite(exp(frailty_nodes)))) {
    stop("|s| is too large for stable frailty quadrature.")
  }
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
  ss = validate_cf_split(ss, n)
  
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
  
  ## Fit cause-specific coefficients using the same design in coxme and coxph.
  MXdf <- data.frame(M, X)
  MXn <- names(MXdf)
  p <- ncol(MXdf)
  df_c <- cbind(rbind(MXdf, MXdf), rbind(MXdf, MXdf))
  df_c[n + seq_len(n), seq_len(p)] <- 0
  df_c[seq_len(n), p + seq_len(p)] <- 0
  names(df_c) <- c(paste0(MXn, '1'), MXn)
  if (anyDuplicated(names(df_c))) stop("Covariate names collide in the cause-specific design.")
  predictor_names <- names(df_c)
  reserved_names <- c('Time2', 'status2', 'cause', 'id', 'case_weight', 'frailty_load')
  if (any(predictor_names %in% reserved_names)) {
    stop('Covariate names collide with reserved survival-model fields.')
  }
  df_c$Time2 <- rep(Time, 2)
  df_c$status2 <- c(as.integer(status == 1), as.integer(status == 2))
  df_c$cause <- factor(rep(c(1, 2), each = n))
  df_c$id <- rep(seq_len(n), 2)
  df_c$case_weight <- rep(wts, 2)
  # Loading +1 on cause 1 and cause2_sign on cause 2; (load | id) is a random slope,
  # so it reduces to the original (1 | id) shared frailty whenever cause2_sign is +1.
  df_c$frailty_load <- rep(c(1, cause2_sign), each = n)
  rhs <- paste(sprintf('`%s`', predictor_names), collapse = ' + ')
  random_term <- if (cause2_sign == 1) '(1 | id)' else '(frailty_load | id)'
  form <- as.formula(paste('Surv(Time2, status2) ~', rhs, '+ strata(cause) +', random_term))
  form2 <- as.formula(paste('Surv(Time2, status2) ~', rhs, '+ strata(cause)'))
  folds <- list(ss, !ss)
  event_fits <- list(vector('list', 2), vector('list', 2))
  event_lam <- list(list(vector('list', 2), vector('list', 2)),
                    list(vector('list', 2), vector('list', 2)))
  event_haz <- list(list(matrix(0, n, mt), matrix(0, n, mt)),
                    list(matrix(0, n, mt), matrix(0, n, mt)))
  censor_haz <- list(matrix(0, n, mt), matrix(0, n, mt))
  hazard_increments <- function(bh) {
    cumulative <- c(0, bh$hazard)[findInterval(tt, bh$time) + 1L]
    diff(c(0, cumulative))
  }
  for (fold in seq_along(folds)) {
    train <- folds[[fold]]
    test <- which(!train)
    for (a in 0:1) {
      keep <- train & A == a & D
      if (!any(keep) || any(tabulate(status[keep] + 1L, nbins = 3L) == 0L)) {
        stop(sprintf('Fold %d, A=%d needs observed M and each event/censoring type.', fold, a))
      }
      train_data <- df_c[rep(keep, 2), , drop = FALSE]
      if (s == 0) {
        event_fit <- coxph(form2, data=train_data, weights=case_weight, model=TRUE)
      } else {
        frailty_fit <- coxme::coxme(form, data=train_data, weights=case_weight,
                                    vfixed=frailty_var, ties='efron')
        beta_all <- coxme::fixef(frailty_fit)[predictor_names]
        if (any(!is.finite(beta_all))) stop('Non-finite frailty coefficients.')
        # Preserve the original fixed-coefficient baseline approximation (no BLUP offset).
        event_fit <- coxph(form2, data=train_data, weights=case_weight,
                           init=unname(beta_all), control=coxph.control(iter.max=0), model=TRUE)
      }
      event_fits[[a + 1L]][[fold]] <- event_fit
      bh <- basehaz(event_fit, centered=FALSE)
      for (cause in 1:2) {
        lam <- hazard_increments(bh[as.character(bh$strata) == as.character(cause), , drop=FALSE])
        lp <- predict(event_fit, df_c[(cause - 1L) * n + test, , drop=FALSE],
                      type='lp', reference='zero')
        event_lam[[a + 1L]][[cause]][[fold]] <- lam
        event_haz[[a + 1L]][[cause]][test, ] <- exp(lp) %*% t(lam)
      }
      censor_fit <- coxph(Surv(Time, status == 0) ~ ., data=MXdf, weights=wts, subset=keep)
      lam <- hazard_increments(basehaz(censor_fit, centered=FALSE))
      lp <- predict(censor_fit, MXdf[test, , drop=FALSE], type='lp', reference='zero')
      censor_haz[[a + 1L]][test, ] <- exp(lp) %*% t(lam)
    }
  }
  fit0_cf <- event_fits[[1]]; fit1_cf <- event_fits[[2]]
  lam01_cf <- event_lam[[1]][[1]]; lam02_cf <- event_lam[[1]][[2]]
  lam11_cf <- event_lam[[2]][[1]]; lam12_cf <- event_lam[[2]][[2]]
  haz01 <- event_haz[[1]][[1]]; haz02 <- event_haz[[1]][[2]]
  haz11 <- event_haz[[2]][[1]]; haz12 <- event_haz[[2]][[2]]
  haz00 <- censor_haz[[1]]; haz10 <- censor_haz[[2]]
  
  # fit propensity score
  ps_a = rep(NA, n)
  fit_a = glm(A~., data=data.frame(A,X), family='binomial', subset=ss)
  ps_a[!ss] = predict(fit_a, newdata=data.frame(A,X)[!ss,,drop=FALSE], type='response')
  fit_a = glm(A~., data=data.frame(A,X), family='binomial', subset=!ss)
  ps_a[ss] = predict(fit_a, newdata=data.frame(A,X)[ss,,drop=FALSE], type='response')
  ps = ps_a*A + (1-ps_a)*(1-A)
  ps[A==1] = ps[A==1] * mean(A/ps)
  ps[A==0] = ps[A==0] * mean((1-A)/ps)
  Xdf = data.frame(X)
  
  # fit mediator density
  if (method=='ratio'){
    mx = lm(M ~ ., data=data.frame(A,X,AX=A*X), weights=wts, subset=(D==1&ss))
    res_m = rep(NA, n)
    res_m[D==1&ss] = residuals(mx)^2
    mx1 = predict(mx, newdata=data.frame(A=1,X,AX=X)[!ss,])
    mx0 = predict(mx, newdata=data.frame(A=0,X,AX=0*X)[!ss,])
    sx = lm(log(res_m+0.000001) ~ ., data=data.frame(A,X,AX=A*X), weights=wts, subset=(D==1&ss))
    mediator_logvar = mean(residuals(sx)^2)
    sx1 = exp(predict(sx, newdata=data.frame(A=1,X,AX=X)[!ss,]) + mediator_logvar/2)
    sx0 = exp(predict(sx, newdata=data.frame(A=0,X,AX=0*X)[!ss,]) + mediator_logvar/2)
    mx_1ss = cbind(mx0, mx1)
    sx_1ss = sqrt(cbind(sx0, sx1))
    mx = lm(M ~ ., data=data.frame(A,X,AX=A*X), weights=wts, subset=(D==1&!ss))
    res_m = rep(NA, n)
    res_m[D==1&!ss] = residuals(mx)^2
    mx1 = predict(mx, newdata=data.frame(A=1,X,AX=X)[ss,])
    mx0 = predict(mx, newdata=data.frame(A=0,X,AX=0*X)[ss,])
    sx = lm(log(res_m+0.000001) ~ ., data=data.frame(A,X,AX=A*X), weights=wts, subset=(D==1&!ss))
    mediator_logvar = mean(residuals(sx)^2)
    sx1 = exp(predict(sx, newdata=data.frame(A=1,X,AX=X)[ss,]) + mediator_logvar/2)
    sx0 = exp(predict(sx, newdata=data.frame(A=0,X,AX=0*X)[ss,]) + mediator_logvar/2)
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
  a_seq = rbind(c(0,0,0), c(1,0,0), c(1,1,0), c(1,1,1))
  F_reg = F_eff = se_eff = matrix(0, nrow=4, ncol=length(tt))
  a_names = apply(a_seq, 1, paste0, collapse="")
  rownames(F_reg) = rownames(F_eff) = rownames(se_eff) = a_names
  eif = array(0, dim=c(n,length(tt),4), dimnames=list(NULL,round(tt,4),a_names))
  
  row_cumsum <- function(x) t(apply(x, 1L, cumsum))
  haz0 <- cbind(0, haz10*A + haz00*(1-A))[, seq_len(mt), drop=FALSE]
  observed_haz1 <- haz11*A + haz01*(1-A)
  observed_haz2 <- haz12*A + haz02*(1-A)
  observed_cumhaz1 <- row_cumsum(observed_haz1)
  observed_cumhaz2 <- row_cumsum(observed_haz2)
  censor_cumhaz <- row_cumsum(haz0)
  Yt <- outer(Time, tt, '>=')
  dNt <- outer(Time, tt, '==')
  dN1 <- dNt * (status == 1)
  dN2 <- dNt * (status == 2)
  rm(dNt)

  # Cache quantities independent of the quadrature node.
  effect_cache <- vector('list', nrow(a_seq))
  for (k in seq_len(nrow(a_seq))) {
    a0 <- a_seq[k, 1]; a1 <- a_seq[k, 2]; a2 <- a_seq[k, 3]
    haz1 <- event_haz[[a1 + 1L]][[1]]
    haz2 <- event_haz[[a2 + 1L]][[2]]
    rm01 <- rm02 <- rep(1, n)
    if (method == 'ratio') {
      if (a0 != a1) rm01 <- dnorm(M, mx[,a0+1], sx[,a0+1]) / dnorm(M, mx[,a1+1], sx[,a1+1])
      if (a0 != a2) rm02 <- dnorm(M, mx[,a0+1], sx[,a0+1]) / dnorm(M, mx[,a2+1], sx[,a2+1])
    } else if (method == 'propensity') {
      if (a0 != a1) rm01 <- ps_m[,a0+1] / ps_m[,a1+1]
      if (a0 != a2) rm02 <- ps_m[,a0+1] / ps_m[,a2+1]
    } else {
      if (a0 != a1) rm01 <- (a1 == 1)/r10 + (a1 == 0)*r10
      if (a0 != a2) rm02 <- (a2 == 1)/r10 + (a2 == 0)*r10
    }
    regression_cache <- lapply(seq_along(folds), function(fold) {
      idx <- which(A == a0 & D & folds[[fold]])
      lp1 <- predict(event_fits[[a1 + 1L]][[fold]], df_c[idx, , drop=FALSE],
                     type='lp', reference='zero')
      lp2 <- predict(event_fits[[a2 + 1L]][[fold]], df_c[n + idx, , drop=FALSE],
                     type='lp', reference='zero')
      h1 <- exp(lp1) %*% t(event_lam[[a1 + 1L]][[1]][[fold]])
      h2 <- exp(lp2) %*% t(event_lam[[a2 + 1L]][[2]][[fold]])
      list(idx=idx, test=which(!folds[[fold]]), h1=h1,
           cumhaz1=row_cumsum(h1), cumhaz2=row_cumsum(h2))
    })
    effect_cache[[k]] <- list(
      haz1=haz1, cumhaz1=row_cumsum(haz1), cumhaz2=row_cumsum(haz2),
      weight1=(A == a1)/ps * pmin(99, pmax(0.01, rm01)),
      weight2=(A == a2)/ps * pmin(99, pmax(0.01, rm02)),
      weight0=(A == a0)/ps, regression=regression_cache
    )
  }

  for (b in seq_along(frailty_nodes)) {
    # One latent draw, loaded +1 on cause 1 and cause2_sign on cause 2.
    multiplier1 <- exp(frailty_nodes[b])
    multiplier2 <- exp(cause2_sign * frailty_nodes[b])
    node_weight <- quad_weights[b]
    S <- pmax(0.01, exp(-observed_cumhaz1 * multiplier1 -
                          observed_cumhaz2 * multiplier2 - censor_cumhaz))
    dG1t <- (dN1 - Yt * observed_haz1 * multiplier1) / S
    dG2t <- (dN2 - Yt * observed_haz2 * multiplier2) / S
    G1t <- row_cumsum(dG1t)
    G2t <- row_cumsum(dG2t)
    for (k in seq_len(nrow(a_seq))) {
      cached <- effect_cache[[k]]
      haz1 <- cached$haz1 * multiplier1
      S <- exp(-cached$cumhaz1 * multiplier1 - cached$cumhaz2 * multiplier2)
      term1 <- cached$weight1 * (dG1t - G1t * haz1)
      term2 <- -cached$weight2 * G2t * haz1
      term3 <- cached$weight0 * haz1
      gamma_obs <- matrix(NA_real_, n, mt)
      h_hat <- matrix(0, n, mt)
      for (fold in seq_along(folds)) {
        reg <- cached$regression[[fold]]
        gamma_obs[reg$idx, ] <- exp(-reg$cumhaz1 * multiplier1 -
                                      reg$cumhaz2 * multiplier2) * reg$h1 * multiplier1
        fit_h <- lm(gamma_obs ~ ., data=Xdf, weights=wts, subset=reg$idx)
        h_hat[reg$test, ] <- predict(fit_h, newdata=Xdf[reg$test, , drop=FALSE])
      }
      h_hat[h_hat < 0] <- 0
      F_reg[k, ] <- F_reg[k, ] + colMeans(row_cumsum(h_hat), na.rm=TRUE) * node_weight
      term4 <- (1 - cached$weight0) * h_hat
      # Legacy phi fields remain conditional on the last node; they are not marginal EIFs.
      if (k == 1L) {
        phi0 <- term1 + term2; phi00 <- term1
      }
      if (k == 4L) {
        phi1 <- term1 + term2; phi10 <- term1
      }
      Fx <- row_cumsum(S * (term1 + term2 + term3) + term4)
      Fx[!D, ] <- 0
      eif[, , k] <- eif[, , k] + Fx * node_weight
    }
  }

  # Weighted least squares is linear in its response: integrate before augmentation.
  # aug_time = FALSE drops Time from the augmentation regression, as in fit_eif_cf
  augmentation_data <- if (aug_time) data.frame(A, X, Time, st=as.factor(status)) else
    data.frame(A, X, st=as.factor(status))
  for (k in seq_len(nrow(a_seq))) {
    Fx <- eif[, , k]
    if (any(!D)) {
      EFx <- matrix(0, n, mt)
      for (fold in seq_along(folds)) {
        train <- D & folds[[fold]]
        test <- which(!folds[[fold]])
        fit_phi <- lm(Fx ~ ., data=augmentation_data, weights=wts, subset=train)
        EFx[test, ] <- predict(fit_phi, newdata=augmentation_data[test, , drop=FALSE])
      }
      Fx <- wts * Fx + (1 - wts) * EFx
    }
    if (any(!is.finite(Fx))) stop('Non-finite integrated scores; inspect overlap and model fits.')
    F_effk <- pmin(1, pmax(0, colMeans(Fx)))
    if (monotone) F_effk <- isoreg(F_effk)$yf
    eifk <- sweep(Fx, 2L, F_effk, '-')
    F_eff[k, ] <- F_effk
    se_eff[k, ] <- sqrt(colMeans(eifk^2)) / sqrt(n - 1)
    eif[, , k] <- eifk
  }
  colnames(F_reg) <- colnames(F_eff) <- round(tt, 4)

  list(data=dat, time=tt, phi0=phi0, phi1=phi1, phi00=phi00, phi10=phi10,
       eif=eif, wts=wts, ps=ps,
       lam01=lam01_cf, lam02=lam02_cf, fit01=fit0_cf, fit02=fit0_cf,
       lam11=lam11_cf, lam12=lam12_cf, fit11=fit1_cf, fit12=fit1_cf,
       F_reg=F_reg, F_eff=F_eff, se=se_eff, cf=TRUE, ss=ss,
       frailty_s=s, frailty_sd=abs(s), cause2_sign=cause2_sign,
       frailty_variance=frailty_var,
       quadrature=list(nodes=frailty_nodes, weights=quad_weights),
       phi_frailty_node=tail(frailty_nodes, 1))
}

# Direct and mediation effects on the restricted mean time to first infection of type
# `event` as the frailty effect s varies. Each band is the estimate +/- 1.96 x the
# robust (MAD) bootstrap SD from app1_sens_boot.R (B = 100 replicates of the whole
# procedure, including the no-frailty TE that sets the scale factor); the analytic
# intervals are kept in the results for comparison.
run_sens <- function(res, s_seq = seq(-1, 1, length.out = 11), seed = 2025,
                     boot_file = 'app1_sens_boot.rds') {
  fit <- res$fit; natural <- res$natural; event <- res$event
  sens <- do.call(rbind, lapply(s_seq, function(s) {
    set.seed(seed)
    fit_s <- fit_eif_cf_f(synth$trt, synth$Time, synth$status, synth$M, synth[, Xname],
                          method = 'riesz', event = event, s = s, ss = fit$ss,
                          aug_time = FALSE)
    effects <- est_decompose(fit_s)
    data.frame(s = s, variance = s^2,
               TE = effects$RMSTL_TE, ME = effects$RMSTL_ME,
               DE = effects$RMSTL_DE, CE = effects$RMSTL_CE,
               TE_se = effects$RMSTL_TE_se, ME_se = effects$RMSTL_ME_se,
               DE_se = effects$RMSTL_DE_se, CE_se = effects$RMSTL_CE_se)
  }))
  TE <- natural$RMSTL_TE
  DCE <- natural$RMSTL_DE + natural$RMSTL_CE
  # Retain the original scaling; do not interpret it as the unscaled frailty-specific DE.
  sens$scale_factor <- -TE / sens$TE
  sens$scaled_direct_RMSTG <- sens$DE * sens$scale_factor
  # The scale factor is held fixed, so these intervals carry only the DE(s) sampling error.
  sens$scaled_se <- sens$DE_se * abs(sens$scale_factor)
  sens$scaled_lower <- sens$scaled_direct_RMSTG - qnorm(0.975) * sens$scaled_se
  sens$scaled_upper <- sens$scaled_direct_RMSTG + qnorm(0.975) * sens$scaled_se
  unstable <- !is.finite(sens$scaled_direct_RMSTG) |
    abs(sens$TE) <= sqrt(.Machine$double.eps) * max(1, abs(TE))
  if (any(unstable)) {
    warning('Some TE denominators are numerically near zero; scaled effects are set to NA.')
    sens[unstable, c('scale_factor', 'scaled_direct_RMSTG', 'scaled_se',
                     'scaled_lower', 'scaled_upper')] <- NA_real_
  }

  # The mediation effect gets the same scaling as the direct effect
  sens$scaled_mediation_RMSTG <- sens$ME * sens$scale_factor
  sens[unstable, 'scaled_mediation_RMSTG'] <- NA_real_
  sens_plots(sens, readRDS(boot_file), natural, event)
}

# Bootstrap bands for the scaled direct and mediation effects (replicate b uses its own
# scale factor -TE_b(0) / TE_b(s)) and the two sensitivity figures. Only
# the bands are recomputed here, so new figures can be drawn from saved fits.
sens_plots <- function(sens, boot, natural, event) {
  boot <- boot[boot$event == event, ]
  boot <- merge(boot, setNames(boot[boot$s == 0, c('b', 'TE')], c('b', 'TE0')), by = 'b')
  sf <- -boot$TE0 / boot$TE
  boot$scaled_DE <- boot$DE * sf
  boot$scaled_ME <- boot$ME * sf
  # Normal intervals, estimate +/- 1.96 x bootstrap SD. Percentile intervals are not used:
  # the bootstrap replicates of the asymptomatic mediation and direct effects are shifted
  # by 1-2 days from the estimates, and percentile limits inherit that shift.
  # The SD is the robust 1.4826 x MAD: in a few replicates TE_b(s) is almost exactly zero,
  # so the scale factor reaches 10-120 and a single replicate would dominate the plain SD.
  bq <- do.call(rbind, lapply(split(boot, boot$s), function(x) data.frame(
    s = x$s[1], B = sum(is.finite(x$scaled_DE)),
    boot_sd = mad(x$scaled_DE, na.rm = TRUE), boot_me_sd = mad(x$scaled_ME, na.rm = TRUE))))
  sens <- merge(sens[, setdiff(names(sens), c(setdiff(names(bq), 's'), 'boot_lower',
                                              'boot_upper', 'boot_me_lower',
                                              'boot_me_upper'))], bq, by = 's')
  z <- qnorm(0.975)
  sens$boot_lower <- sens$scaled_direct_RMSTG - z * sens$boot_sd
  sens$boot_upper <- sens$scaled_direct_RMSTG + z * sens$boot_sd
  sens$boot_me_lower <- sens$scaled_mediation_RMSTG - z * sens$boot_me_sd
  sens$boot_me_upper <- sens$scaled_mediation_RMSTG + z * sens$boot_me_sd

  TE <- natural$RMSTL_TE
  DCE <- natural$RMSTL_DE + natural$RMSTL_CE
  ref <- data.frame(label = c('TE (No frailty)', 'DE+CE (No frailty)'), RMSTG = c(-TE, -DCE))
  one_plot <- function(est, lo, hi, effect_name) {
    yl <- range(c(sens[[lo]], sens[[hi]], ref$RMSTG, 0), na.rm = TRUE)
    ggplot(sens, aes(x = s, y = .data[[est]])) +
      geom_ribbon(aes(ymin = .data[[lo]], ymax = .data[[hi]]), alpha = 0.15, fill = 'brown') +
      geom_line(linewidth = 1.2, colour = 'brown') +
      geom_hline(data = ref, aes(yintercept = RMSTG, colour = label, linetype = label),
                 linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = 'dotted') +
      scale_colour_manual(values = c('DE+CE (No frailty)' = 'steelblue',
                                     'TE (No frailty)' = 'black')) +
      scale_linetype_manual(values = c('DE+CE (No frailty)' = 'longdash',
                                       'TE (No frailty)' = 'solid')) +
      labs(x = 'Frailty effect', y = 'Restricted mean time delayed', colour = NULL,
           linetype = NULL,
           title = paste(effect_name, 'effect on', event_names[event], 'first infection')) +
      coord_cartesian(ylim = c(-ceil_to(-yl[1], 5), ceil_to(yl[2], 5))) +
      theme_bw() +
      theme(legend.position = 'bottom')
  }
  list(results = sens,
       plot = one_plot('scaled_direct_RMSTG', 'boot_lower', 'boot_upper', 'Direct'),
       plot_me = one_plot('scaled_mediation_RMSTG', 'boot_me_lower', 'boot_me_upper',
                          'Mediation'))
}

app1_sens <- lapply(app1_res, run_sens)
for (k in 1:2) {
  save_app1(app1_sens[[k]]$plot, paste0('app1_sensitivity', k))
  save_app1(app1_sens[[k]]$plot_me, paste0('app1_sensitivity_me', k))
}
lapply(app1_sens, function(x) round(x$results[, c('s', 'TE', 'DE', 'scaled_direct_RMSTG',
                                                  'boot_lower', 'boot_upper', 'ME',
                                                  'scaled_mediation_RMSTG',
                                                  'boot_me_lower', 'boot_me_upper')], 2))

## ===========================================================================
## AMP applications (HVTN 703/HPTN 081 and HVTN 704/HPTN 085): controlled CIF under VRC01
## ===========================================================================
## A = VRC01 (10 and 30 mg/kg pooled) vs placebo; M = log10 pt80. X = age_bin, age
## above the median among participants with a measured mediator.
## Primary event: event_cat 1-2 (VRC01-sensitive); competing event: event_cat 3.
## PT80 is measured only in the popPK subcohort, so the placebo arm contributes no
## mediator data and the cross-world arms are not identified. The controlled CIF needs
## only the (1,1,1) path, which fit_eif_cf reports via arm_ok; the remaining arms come
## back as NA with a warning.
## Step-like flat stretches are the monotone projection pooling a non-monotone raw
## estimate, not risk that stops accruing. Switching method= to 'propensity' changes
## nothing here: the (1,1,1) arm has a0 = a1 = a2, so the density ratios are fixed at 1.
## Mediator: pt80 = day 1-7 max predicted VRC01 concentration / placebo-case
## geometric-mean IC80 of the same trial, added to the CSV by build_pt80.R. The file's
## pt80_panel (nab_max2 / panel_ic80) is a different definition and is not used.
##
## HVTN 703: 88 treated participants with a measured mediator (11 primary, 32 competing
## events). At day 601 the EIF estimate of the marginal (1,1,1) CIF is 0.0145 against
## 0.0100 for both F_reg and the treated-arm Aalen-Johansen value; read the shape across
## PT80 more than the absolute level. Marker grid: 7 points equally spaced in log10
## PT80 over 15-75, inside the treated range (14.5-94.7). PT80 is bimodal by dose
## (10 mg/kg: 14.5-30.4; 30 mg/kg: 42.6-94.7), so the grid point at 33.5 lies in the gap
## and is supported only by kernel smoothing (h = 0.2 on the log10 scale).
## HVTN 704: 99 treated participants with a measured mediator (17 primary, 38 competing
## events). PT80 is higher than in 703 (10 mg/kg: 40-87; 30 mg/kg: 85-234), and the grid
## spans its 10%-90% quantiles (46-171).
amp <- read.csv('amp_survival_nab_analysis.csv')
stopifnot('pt80' %in% names(amp))

prep_amp <- function(protocol_name) {
  d <- amp |>
    filter(protocol %in% protocol_name, !is.na(event_cat)) |>
    mutate(A      = as.integer(dose != 'Control'),
           M      = log10(pt80),
           Time   = survday,
           status = case_when(event_cat %in% 1:2 ~ 1L,
                              event_cat == 3     ~ 2L,
                              .default           = 0L))
  # median split among measured participants keeps the two age groups balanced; taken
  # within each trial when several are pooled
  d <- d |>
    mutate(age_bin = as.integer(age > median(age[!is.na(M)])), .by = protocol) |>
    mutate(trial704 = as.integer(protocol == 'HVTN 704/HPTN 085'))
  d
}

# Placebo-arm CIF of the primary event, E[F1(t | A = 0, X)] over the whole trial, by the
# one-step (EIF-based) semiparametric estimator. Cox working models in the placebo arm
# for the primary event, the competing event and censoring, each with covariates xvars
# (discrete); P(A = 0 | X) by stratum proportions. xvars = NULL gives the unadjusted
# estimator, which reproduces Aalen-Johansen and its SE exactly.
# EIF: F1(t|X) - psi + (1-A)/pi0(X) [ int_0^t dM1/G(v-|X)
#                      - int_0^t {F1(t|X) - F1(v|X)} / {S(v|X) G(v-|X)} dM ].
eff_cif_control <- function(d, xvars = NULL, cause = 1) {
  n <- nrow(d)
  Xall <- if (is.null(xvars)) matrix(0, n, 1) else as.matrix(d[, xvars, drop = FALSE])
  key <- apply(Xall, 1, paste, collapse = '_')
  pats <- Xall[!duplicated(key), , drop = FALSE]
  pid <- match(key, key[!duplicated(key)])
  ctl <- d$A == 0
  d0 <- d[ctl, ]
  tt <- sort(unique(d0$Time[d0$Time > 0]))
  J <- length(tt)
  dL <- function(ev) {
    df <- data.frame(Time = d0$Time, ev = ev, Xall[ctl, , drop = FALSE])
    f <- if (is.null(xvars)) coxph(Surv(Time, ev) ~ 1, data = df, ties = 'breslow') else
      coxph(Surv(Time, ev) ~ ., data = df, ties = 'breslow')
    bh <- basehaz(f, centered = FALSE)
    dH0 <- diff(c(0, stepfun(bh$time, c(0, bh$hazard))(tt)))
    lp <- if (is.null(xvars)) rep(0, nrow(pats)) else drop(pats %*% coef(f))
    exp(lp) %o% dH0                                    # covariate patterns x time
  }
  dL1 <- dL(d0$status == cause)
  dL2 <- dL(d0$status != cause & d0$status != 0)
  dLc <- dL(d0$status == 0)
  rcum <- function(m, f) matrix(apply(m, 1, f), nrow = nrow(m), byrow = TRUE)
  S <- rcum(1 - dL1 - dL2, cumprod)                    # S(v | x)
  Sm <- cbind(1, S[, -J, drop = FALSE])                # S(v- | x)
  Gm <- cbind(1, rcum(1 - dLc, cumprod)[, -J, drop = FALSE])   # G(v- | x)
  F1 <- rcum(Sm * dL1, cumsum)                         # F1(t | x)
  xi <- pid[ctl]
  Y <- outer(d0$Time, tt, '>=')
  ev_at <- outer(d0$Time, tt, '==')
  dM1 <- ev_at * (d0$status == cause) - Y * dL1[xi, , drop = FALSE]
  dM <- ev_at * (d0$status != 0) - Y * (dL1[xi, , drop = FALSE] + dL2[xi, , drop = FALSE])
  W <- S[xi, , drop = FALSE] * Gm[xi, , drop = FALSE]
  aug0 <- rcum(dM1 / Gm[xi, , drop = FALSE], cumsum) -
    (F1[xi, , drop = FALSE] * rcum(dM / W, cumsum) - rcum(F1[xi, , drop = FALSE] * dM / W, cumsum))
  pi0 <- as.numeric(tapply(ctl, pid, mean))
  aug <- matrix(0, n, J)
  aug[ctl, ] <- aug0 / pi0[xi]
  phi <- F1[pid, , drop = FALSE] + aug
  psi <- colMeans(phi)
  IF <- sweep(phi, 2, psi)
  list(time = tt, psi = psi, se = sqrt(colMeans(IF^2) / n), IF = IF)
}

# Step-function lookup of the placebo-arm estimate at arbitrary times
ctl_at <- function(ctl, times, what = c('psi', 'se')) {
  what <- match.arg(what)
  c(0, ctl[[what]])[findInterval(times, ctl$time) + 1]
}

# Controlled counterfactual CIF of VRC01-sensitive acquisition under VRC01, with the
# marker fixed on m_seq; draws the CIF over time at the lowest, middle and highest grid
# points (panel a) and the end-of-follow-up CIF against the marker (panel b). The
# placebo-arm reference (panel a, and the denominator of relative efficacy) is the
# semiparametric estimator above with the same covariates.
## protocol_name may list several trials for a pooled analysis; x_names then needs the
## trial indicator 'trial704'.
run_amp <- function(protocol_name, m_seq, file_stub = NULL, h = 0.2, monotone = TRUE,
                    seed = 2025, x_names = 'age_bin') {
  # est_con selects the covariate columns of fit$data through the global Xname
  Xname <<- x_names
  d <- prep_amp(protocol_name)
  set.seed(seed)
  fit <- fit_eif_cf(d$A, d$Time, d$status, d$M, d[, Xname, drop = FALSE],
                    method = 'riesz', event = 1, monotone = monotone, aug_time = FALSE)
  con <- est_con(fit, a = 1, m_seq = m_seq, h = h, monotone = monotone)
  tmax <- max(fit$time)
  idx3 <- c(1, ceiling(length(m_seq) / 2), length(m_seq))

  groups <- paste('Predicted ID80 =', round(10^m_seq[idx3]))
  ref_lab <- 'Control'
  group_cols <- setNames(c('black', 'brown', 'darkorange', 'grey40'), c(groups, ref_lab))
  group_lty <- setNames(c(rep('solid', 3), 'dashed'), c(groups, ref_lab))
  # Reference: semiparametric placebo-arm CIF of the primary event over the whole trial
  ctl <- eff_cif_control(d, Xname)
  if (monotone) ctl$psi <- isoreg(ctl$psi)$yf
  df_ref <- data.frame(time = c(0, ctl$time), CIF = c(0, ctl$psi),
                       group = factor(ref_lab, levels = c(groups, ref_lab)))
  time_vec <- c(0, con$time)
  CIF_mat <- do.call(rbind, lapply(idx3, function(i) {
    f <- c(0, con$Fm_eff[i, ])
    if (monotone) isoreg(f)$yf else f
  }))
  SE_mat <- cbind(0, con$se[idx3, , drop = FALSE])
  df_est <- data.frame(
    time  = rep(time_vec, times = length(idx3)),
    CIF   = as.vector(t(CIF_mat)),
    lower = pmax(0, as.vector(t(CIF_mat - 1.96 * SE_mat))),
    upper = pmin(1, as.vector(t(CIF_mat + 1.96 * SE_mat))),
    group = factor(rep(groups, each = length(time_vec)), levels = c(groups, ref_lab))
  )
  p_cif <- ggplot(df_est, aes(x = time, y = CIF, color = group, linetype = group)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.2, colour = NA,
                show.legend = FALSE) +
    geom_line(linewidth = 0.7) +
    geom_step(data = df_ref, linewidth = 0.7) +
    scale_color_manual(values = group_cols, name = 'Marker',
                       guide = guide_legend(position = 'inside')) +
    scale_linetype_manual(values = group_lty, name = 'Marker') +
    scale_fill_manual(values = group_cols, guide = 'none') +
    coord_cartesian(ylim = c(0, ceil_to(c(df_est$upper, df_ref$CIF), 0.01))) +
    labs(x = 'Study time (days)', y = 'Cumulative incidence', title = '(a)') +
    theme_bw() +
    theme(legend.position.inside = c(0.22, 0.78),
          legend.background = element_rect(fill = alpha('white', 0.6)))

  Jm <- ncol(con$Fm_eff)
  df_all <- data.frame(
    marker = 10^m_seq,
    value  = con$Fm_eff[, Jm],
    lower  = pmax(0, con$Fm_eff[, Jm] - 1.96 * con$se[, Jm]),
    upper  = pmin(1, con$Fm_eff[, Jm] + 1.96 * con$se[, Jm])
  )
  # log10 axis labelled on the PT80 scale, as in the COVID CIF-vs-marker plot
  p_curve <- ggplot(df_all, aes(x = marker, y = value)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line(linewidth = 0.6) +
    scale_x_log10(breaks = 10^m_seq, labels = round(10^m_seq, 1)) +
    coord_cartesian(ylim = c(0, ceil_to(df_all$upper, 0.01))) +
    labs(x = expression('log'[10] * ' Predicted ID80'),
         y = paste0('Cumulative incidence at day ', tmax), title = '(b)') +
    theme_bw()

  res <- list(data = d, fit = fit, con = con, control = ctl, tmax = tmax, cif_plot = p_cif,
              curve_plot = p_curve, end_cif = df_all)
  re <- amp_relative_efficacy(res, monotone = monotone)
  res$re <- re$data
  res$re_plot <- re$plot
  re_hi <- amp_relative_efficacy(res, monotone = monotone, which = 'high', panel = '(d)')
  res$re_high <- re_hi$data
  res$re_high_plot <- re_hi$plot
  if (!is.null(file_stub)) {
    ggsave(paste0(file_stub, '_cif1.pdf'), p_cif, width = 400 / 72, height = 275 / 72)
    ggsave(paste0(file_stub, '_curve1.pdf'), p_curve, width = 400 / 72, height = 275 / 72)
    ggsave(paste0(file_stub, '_re1.pdf'), re$plot, width = 400 / 72, height = 275 / 72)
    ggsave(paste0(file_stub, '_re1_high.pdf'), re_hi$plot, width = 400 / 72, height = 275 / 72)
  }
  res
}

# Marker-specific relative efficacy RE(t) = 1 - F^{a=1,m}(t) / F^0(t), with m the middle
# grid point (the middle line of panel a) and F^0 the semiparametric placebo-arm CIF
# (res$control). The delta method on the log ratio, treating the two estimates as
# independent, gives se{log(F1/F0)} = sqrt(se1^2/F1^2 + se0^2/F0^2); the interval is
# transformed back, which
# keeps its upper end below 1. The curve starts once the placebo arm has min_events
# primary events: before that F^0 rests on one or two events and the band is unbounded
# below.
## which = 'middle' or 'high' picks the middle or highest grid point (the middle or
## orange line of panel a); the curve also starts only once F^{a=1,m} is positive.
amp_relative_efficacy <- function(res, monotone = TRUE, min_events = 5,
                                  which = c('middle', 'high'), panel = '(c)') {
  which <- match.arg(which)
  con <- res$con; d <- res$data
  k <- if (which == 'middle') ceiling(length(con$m_seq) / 2) else length(con$m_seq)
  F1 <- con$Fm_eff[k, ]
  if (monotone) F1 <- isoreg(F1)$yf
  n0 <- sapply(con$time, function(t) sum(d$A == 0 & d$status == 1 & d$Time <= t))
  df <- data.frame(time = con$time, F1 = F1, se1 = con$se[k, ],
                   F0 = ctl_at(res$control, con$time, 'psi'),
                   se0 = ctl_at(res$control, con$time, 'se'), n0_events = n0)
  df <- df[df$n0_events >= min_events & df$F1 > 0, ]
  df$se_log <- sqrt(df$se1^2 / df$F1^2 + df$se0^2 / df$F0^2)
  df$RE <- 1 - df$F1 / df$F0
  df$lower <- 1 - (df$F1 / df$F0) * exp(1.96 * df$se_log)
  df$upper <- 1 - (df$F1 / df$F0) * exp(-1.96 * df$se_log)

  ylo <- -ceil_to(-min(df$lower), 0.25)
  ytop <- 1 + 0.1 * (1 - min(ylo, 0))
  # step-shaped band to match geom_step: hold each bound until the next time point
  n <- nrow(df)
  band <- data.frame(time = c(rbind(df$time, c(df$time[-1], res$tmax))),
                     lower = rep(df$lower, each = 2), upper = rep(df$upper, each = 2))
  p <- ggplot(df, aes(x = time, y = RE)) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_ribbon(data = band, aes(x = time, ymin = lower, ymax = upper), inherit.aes = FALSE,
                alpha = 0.2) +
    geom_step(linewidth = 0.7) +
    # marker label in a strip above RE = 1, which the band (upper bound < 1) never reaches
    annotate('text', x = res$tmax, y = ytop, hjust = 1, vjust = 1,
             label = paste('Predicted ID80 =', round(10^con$m_seq[k]))) +
    coord_cartesian(xlim = c(0, res$tmax), ylim = c(min(ylo, 0), ytop)) +
    labs(x = 'Study time (days)', y = 'Relative efficacy', title = panel) +
    theme_bw()
  list(data = df, plot = p)
}

dir.create('app2', showWarnings = FALSE)
amp703 <- run_amp('HVTN 703/HPTN 081', seq(log10(15), log10(75), length.out = 7),
                  file.path('app2', 'app2'))
m704_q <- quantile(log10(amp$pt80[amp$protocol == 'HVTN 704/HPTN 085' &
                                    !is.na(amp$event_cat)]), c(0.1, 0.9), na.rm = TRUE)
amp704 <- run_amp('HVTN 704/HPTN 085', seq(m704_q[1], m704_q[2], length.out = 7),
                  file.path('app2', 'app2_704'))

## Pooled 703 + 704 with a trial indicator. PT80 is standardised by each trial's own
## placebo-case IC80, and the scales overlap only on about 40-95 (703: 14.5-94.7;
## 704: 40.1-234.1), so grid points outside that range are supported by one trial.
mpool_q <- quantile(log10(amp$pt80[!is.na(amp$event_cat)]), c(0.1, 0.9), na.rm = TRUE)
amp_pool <- run_amp(c('HVTN 703/HPTN 081', 'HVTN 704/HPTN 085'),
                    seq(mpool_q[1], mpool_q[2], length.out = 7),
                    file.path('app2', 'app2_pooled'), x_names = c('age_bin', 'trial704'))
lapply(list(`703` = amp703, `704` = amp704, pooled = amp_pool), function(r)
  list(arm_ok = r$fit$arm_ok, end_cif = round(r$end_cif, 4)))
