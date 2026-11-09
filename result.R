library(statmod)
library(dplyr)
library(survival)
library(ggplot2)
library(xtable)

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
    DM = rbinom(n, 1, plogis(1+(status>0)-0.05*Time))
    M[DM==0] = NA
  }
  return(list(X=X,A=A,M=M,Time=Time,status=status))
}

## Mediation analysis with competing events ##
#setwd('/n0')
B = 1000; n = 5000; method='riesz'
tseq = 1:10; mseq = 0:3; h = 0.2
F_reg_B = F_eff_B = se_B = array(NA, dim=c(4,10,B))
Fnat_reg_B = Fnat_eff_B = senat_B = array(NA, dim=c(3,10,B))
Fcon_reg_B = Fcon_eff_B = secon_B = array(NA, dim=c(4,10,B))
Fcon2_reg_B = Fcon2_eff_B = secon2_B = array(NA, dim=c(4,10,B))
h2 = function(tt,x,m) {0.2*rep(1,nrow(x))%*%t(log(1+tt))}
for (b in 1:B){
  cat(b,'')
  tr = try(load(paste0('fit',b,'.Rdata')))
  if ('try-error' %in% class(tr)) next
  tt = fit$time
  tindex = sapply(tseq, function(k) max(which(tt<=k)))
  F_reg_B[,,b] = fit$F_reg[,tindex]
  F_eff_B[,,b] = fit$F_eff[,tindex]
  se_B[,,b] = fit$se[,tindex]
  Fnat_reg_B[1,,b] = natural$ME_reg[tindex]
  Fnat_reg_B[2,,b] = natural$DE_reg[tindex]
  Fnat_reg_B[3,,b] = natural$CE_reg[tindex]
  Fnat_eff_B[1,,b] = natural$ME_eff[tindex]
  Fnat_eff_B[2,,b] = natural$DE_eff[tindex]
  Fnat_eff_B[3,,b] = natural$CE_eff[tindex]
  senat_B[1,,b] = natural$se_ME[tindex]
  senat_B[2,,b] = natural$se_DE[tindex]
  senat_B[3,,b] = natural$se_CE[tindex]
  se_B[,,b] = fit$se[,tindex]
  Fcon_reg_B[,,b] = control$Fm_reg[,tindex]
  Fcon_eff_B[,,b] = control$Fm_eff[,tindex]
  secon_B[,,b] = control$se[,tindex]
  Fcon2_reg_B[,,b] = control2$Fm_reg[,tindex]
  Fcon2_eff_B[,,b] = control2$Fm_eff[,tindex]
  secon2_B[,,b] = control2$se[,tindex]
}

set.seed(999)
F_true = matrix(NA,nrow=4,ncol=10)
dat = gendata(n=100000, a=c(0,0,0))
F_true[1,] = sapply(tseq, function(k) mean(dat$Time<=k&dat$status==1))
dat = gendata(n=100000, a=c(1,0,0))
F_true[2,] = sapply(tseq, function(k) mean(dat$Time<=k&dat$status==1))
dat = gendata(n=100000, a=c(1,1,0))
F_true[3,] = sapply(tseq, function(k) mean(dat$Time<=k&dat$status==1))
dat = gendata(n=100000, a=c(1,1,1))
F_true[4,] = sapply(tseq, function(k) mean(dat$Time<=k&dat$status==1))
Fa_true = array(F_true,dim=c(4,10,B))
Fnat_true = apply(F_true, 2, diff)
Fnata_true = array(Fnat_true, dim=c(3,10,B))
Fcon_true = Fcon2_true = matrix(0,nrow=length(mseq),ncol=length(tseq))
for (m in 0:3){
  for (b in 1:ngrid){
    m_b = m + sqrt(2)*h*b_list[b]
    dat = gendata(n=100000, a=c(0,0,0), m=m_b)
    Fcon_true[m+1,] = Fcon_true[m+1,] +
      sapply(tseq, function(k) mean(dat$Time<=k&dat$status==1)) * w_list[b]
    dat = gendata(n=100000, a=c(0,0,0), m=m_b, h2=h2)
    Fcon2_true[m+1,] = Fcon2_true[m+1,] +
      sapply(tseq, function(k) mean(dat$Time<=k&dat$status==1)) * w_list[b]
  }
}
Fcona_true = array(Fcon_true, dim=c(4,10,B))
Fcon2a_true = array(Fcon2_true, dim=c(4,10,B))


F_eff_B[F_eff_B<0|F_eff_B>1] = NA
Fnat_eff_B[Fnat_eff_B<(-1)|Fnat_eff_B>1] = NA
biasF_reg = apply(F_reg_B - Fa_true, c(1,2), mean, na.rm=TRUE)
biasF_eff = apply(F_eff_B - Fa_true, c(1,2), mean, na.rm=TRUE)
seF_eff = apply(se_B, c(1,2), mean, na.rm=TRUE)
sdF_eff = apply(F_eff_B, c(1,2), sd, na.rm=TRUE)
cpF_eff = apply(F_eff_B-1.96*se_B<=Fa_true & F_eff_B+1.96*se_B>=Fa_true, c(1,2), 
                mean, na.rm=TRUE)
biasFnat_reg = apply(Fnat_reg_B - Fnata_true, c(1,2), mean, na.rm=TRUE)
biasFnat_eff = apply(Fnat_eff_B - Fnata_true, c(1,2), mean, na.rm=TRUE)
seFnat_eff = apply(senat_B, c(1,2), mean, na.rm=TRUE)
sdFnat_eff = apply(Fnat_eff_B, c(1,2), sd, na.rm=TRUE)
cpFnat_eff = apply(Fnat_eff_B-1.96*senat_B<=Fnata_true & 
                     Fnat_eff_B+1.96*senat_B>=Fnata_true, c(1,2), mean, na.rm=TRUE)
Fcon_eff_B[Fcon_eff_B<0|Fcon_eff_B>1] = NA
Fcon2_eff_B[Fcon2_eff_B<0|Fcon2_eff_B>1] = NA
biasFcon_reg = apply(Fcon_reg_B - Fcona_true, c(1,2), mean, na.rm=TRUE)
biasFcon_eff = apply(Fcon_eff_B - Fcona_true, c(1,2), mean, na.rm=TRUE)
seFcon_eff = apply(secon_B, c(1,2), mean, na.rm=TRUE)
sdFcon_eff = apply(Fcon_eff_B, c(1,2), sd, na.rm=TRUE)
cpFcon_eff = apply(Fcon_eff_B-1.96*secon_B<=Fcona_true & 
                     Fcon_eff_B+1.96*secon_B>=Fcona_true, c(1,2), mean, na.rm=TRUE)
biasFcon2_reg = apply(Fcon2_reg_B - Fcon2a_true, c(1,2), mean, na.rm=TRUE)
biasFcon2_eff = apply(Fcon2_eff_B - Fcon2a_true, c(1,2), mean, na.rm=TRUE)
seFcon2_eff = apply(secon2_B, c(1,2), mean, na.rm=TRUE)
sdFcon2_eff = apply(Fcon2_eff_B, c(1,2), sd, na.rm=TRUE)
cpFcon2_eff = apply(Fcon2_eff_B-1.96*secon2_B<=Fcon2a_true & 
                      Fcon2_eff_B+1.96*secon2_B>=Fcon2a_true, c(1,2), mean, na.rm=TRUE)

# Plot: Counterfactual CIF
groups   <- c('(0,0,0)','(1,0,0)','(1,1,0)','(1,1,1)')
time_vec <- c(0, tseq)
df_true <- data.frame(
  time  = rep(time_vec, times = 4),
  CIF   = c(c(0, F_true[1,]), c(0, F_true[2,]),
            c(0, F_true[3,]), c(0, F_true[4,])),
  group = rep(groups, each = length(time_vec)),
  type  = 'True'
)
df_bias <- data.frame(
  time  = rep(time_vec, times = 4),
  CIF   = c(c(0, F_true[1,] + biasF_eff[1,]),
            c(0, F_true[2,] + biasF_eff[2,]),
            c(0, F_true[3,] + biasF_eff[3,]),
            c(0, F_true[4,] + biasF_eff[4,])),
  group = rep(groups, each = length(time_vec)),
  type  = 'Estimated'
)
df_all <- bind_rows(df_true, df_bias) %>%
  mutate(
    group = factor(group, levels = groups),
    type  = factor(type, levels = c('True', 'Estimated'))
  )
p <- ggplot(df_all, aes(x = time, y = CIF, color = group, linetype = type)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(
    values = c('(0,0,0)' = 'black', '(1,0,0)' = 'darkred', 
               '(1,1,0)' = 'darkgreen', '(1,1,1)' = 'darkblue'),
    name = expression('('*a[0]*','*a[1]*','*a[2]*')'),
    guide = guide_legend(position = 'inside')
  ) +
  scale_linetype_manual(
    values = c('True' = 'solid', 'Estimated' = 'twodash'),
    name = NULL,
    guide = guide_legend(position = 'bottom')
  ) +
  coord_cartesian(ylim = c(0, 0.8)) +
  labs(x = 'Time', y = 'CIF', title = '(a) Counterfactual CIF of Event 1') +
  theme_bw() +
  theme(
    legend.position.inside = c(0.15, 0.75),
    legend.background = element_rect(fill = alpha('white', 0.6))
  )
print(p)

# Plot: Decomposed effects
effect_labels <- c('Direct effect', 'Mediation effect', 'Competing effect')
build_effect_df <- function(i, label) {
  true_vals   <- c(0, Fnat_true[i,])
  biased_vals <- c(0, Fnat_true[i,] + biasFnat_eff[i,])
  se_vals     <- c(0, seFnat_eff[i,])
  upper_vals  <- biased_vals + 1.96 * se_vals
  lower_vals  <- biased_vals - 1.96 * se_vals
  data.frame(
    time   = time_vec,
    value  = c(true_vals, biased_vals, upper_vals, lower_vals),
    curve  = rep(c('True', 'Estimated', 'Upper', 'Lower'), each = length(time_vec)),
    effect = label
  )
}
df_all <- bind_rows(lapply(1:3, function(i) build_effect_df(i, effect_labels[i]))) %>%
  mutate(
    curve  = factor(curve, levels = c('True', 'Estimated', 'Upper', 'Lower')),
    effect = factor(effect, levels = effect_labels)
  )
p <- ggplot(df_all, aes(x = time, y = value, linetype = curve)) +
  geom_line(linewidth = 0.6, color = 'black') +
  scale_linetype_manual(
    values = c('True' = 'solid', 'Estimated' = 'twodash',
               'Upper' = 'dashed', 'Lower' = 'dashed'),
    name = NULL
  ) +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  facet_wrap(~ effect, nrow = 1) +
  labs(x = 'Time', y = 'Treatment effect', title = '(b) Decomposed treatment effects on the CIF') +
  theme_bw() +
  theme(
    strip.text = element_text(face = 'bold'),
    legend.position = 'bottom'
  )
print(p)

# Plot: Counterfactual with mediator (4 curves)
group_labels <- paste('m =', mseq)
df_true <- data.frame(
  time  = rep(time_vec, times = 4),
  CIF   = c(c(0, Fcon_true[1,]), c(0, Fcon_true[2,]),
            c(0, Fcon_true[3,]), c(0, Fcon_true[4,])),
  group = rep(group_labels, each = length(time_vec)),
  type  = 'True'
)
df_bias <- data.frame(
  time  = rep(time_vec, times = 4),
  CIF   = c(c(0, Fcon_true[1,] + biasFcon_eff[1,]),
            c(0, Fcon_true[2,] + biasFcon_eff[2,]),
            c(0, Fcon_true[3,] + biasFcon_eff[3,]),
            c(0, Fcon_true[4,] + biasFcon_eff[4,])),
  group = rep(group_labels, each = length(time_vec)),
  type  = 'Estimated'
)
df_all <- bind_rows(df_true, df_bias) %>%
  mutate(
    group = factor(group, levels = group_labels),
    type  = factor(type, levels = c('True', 'Estimated'))
  )
p <- ggplot(df_all, aes(x = time, y = CIF, color = group, linetype = type)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(
    values = setNames(c('black', 'darkred', 'darkgreen', 'darkblue'), group_labels),
    name = 'm',
    guide = guide_legend(position = 'inside')
  ) +
  scale_linetype_manual(
    values = c('True' = 'solid', 'Estimated' = 'twodash'),
    name = NULL,
    guide = guide_legend(position = 'bottom')
  ) +
  coord_cartesian(ylim = c(0, 0.8)) +
  labs(x = 'Time', y = 'CIF', title = '(c) Controlled counterfactual CIF of Event 1') +
  theme_bw() +
  theme(
    legend.position.inside = c(0.15, 0.75),
    legend.background = element_rect(fill = alpha('white', 0.6)),
    legend.text = element_text(size = 8)
  )
print(p)

df_true <- data.frame(
  time  = rep(time_vec, times = 4),
  CIF   = c(c(0, Fcon2_true[1,]), c(0, Fcon2_true[2,]),
            c(0, Fcon2_true[3,]), c(0, Fcon2_true[4,])),
  group = rep(group_labels, each = length(time_vec)),
  type  = 'True'
)
df_bias <- data.frame(
  time  = rep(time_vec, times = 4),
  CIF   = c(c(0, Fcon2_true[1,] + biasFcon2_eff[1,]),
            c(0, Fcon2_true[2,] + biasFcon2_eff[2,]),
            c(0, Fcon2_true[3,] + biasFcon2_eff[3,]),
            c(0, Fcon2_true[4,] + biasFcon2_eff[4,])),
  group = rep(group_labels, each = length(time_vec)),
  type  = 'Estimated'
)
df_all <- bind_rows(df_true, df_bias) %>%
  mutate(
    group = factor(group, levels = group_labels),
    type  = factor(type, levels = c('True', 'Estimated'))
  )
p <- ggplot(df_all, aes(x = time, y = CIF, color = group, linetype = type)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(
    values = setNames(c('black', 'darkred', 'darkgreen', 'darkblue'), group_labels),
    name = 'm',
    guide = guide_legend(position = 'inside')
  ) +
  scale_linetype_manual(
    values = c('True' = 'solid', 'Estimated' = 'twodash'),
    name = NULL,
    guide = guide_legend(position = 'bottom')
  ) +
  coord_cartesian(ylim = c(0, 0.8)) +
  labs(x = 'Time', y = 'CIF', title = '(e) Controlled counterfactual CIF of Event 1') +
  theme_bw() +
  theme(
    legend.position.inside = c(0.15, 0.75),
    legend.background = element_rect(fill = alpha('white', 0.6)),
    legend.text = element_text(size = 8)
  )
print(p)

# Plot: Counterfactual CIF at t=10 with mediator
true_vals   <- Fcon_true[,10]
biased_vals <- Fcon_true[,10] + biasFcon_eff[,10]
upper_vals  <- biased_vals + 1.96 * seFcon_eff[,10]
lower_vals  <- biased_vals - 1.96 * seFcon_eff[,10]
df_all <- data.frame(
  m     = rep(mseq, times = 4),
  CIF   = c(true_vals, biased_vals, upper_vals, lower_vals),
  curve = rep(c('True', 'Estimated', 'Upper', 'Lower'), each = length(mseq))
) %>%
  mutate(curve = factor(curve, levels = c('True', 'Estimated', 'Upper', 'Lower')))
p <- ggplot(df_all, aes(x = m, y = CIF, linetype = curve)) +
  geom_line(linewidth = 0.6, color = 'black') +
  scale_linetype_manual(
    values = c('True' = 'solid', 'Estimated' = 'twodash',
               'Upper' = 'dashed', 'Lower' = 'dashed'),
    name = NULL,
    breaks = c('True', 'Estimated', 'Upper', 'Lower')
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = 'Mediator', y = 'CIF', title = '(d) Controlled counterfactual CIF at t=10 with m') +
  theme_bw() +
  theme(legend.position = 'bottom')
print(p)

true_vals   <- Fcon2_true[,10]
biased_vals <- Fcon2_true[,10] + biasFcon2_eff[,10]
upper_vals  <- biased_vals + 1.96 * seFcon2_eff[,10]
lower_vals  <- biased_vals - 1.96 * seFcon2_eff[,10]
df_all <- data.frame(
  m     = rep(mseq, times = 4),
  CIF   = c(true_vals, biased_vals, upper_vals, lower_vals),
  curve = rep(c('True', 'Estimated', 'Upper', 'Lower'), each = length(mseq))
) %>%
  mutate(curve = factor(curve, levels = c('True', 'Estimated', 'Upper', 'Lower')))
p <- ggplot(df_all, aes(x = m, y = CIF, linetype = curve)) +
  geom_line(linewidth = 0.6, color = 'black') +
  scale_linetype_manual(
    values = c('True' = 'solid', 'Estimated' = 'twodash',
               'Upper' = 'dashed', 'Lower' = 'dashed'),
    name = NULL,
    breaks = c('True', 'Estimated', 'Upper', 'Lower')
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = 'Mediator', y = 'CIF', title = '(f) Controlled counterfactual CIF at t=10 with m') +
  theme_bw() +
  theme(legend.position = 'bottom')
print(p)


result = rbind(biasF_eff, sdF_eff, seF_eff, cpF_eff)
result = result[c(4*1:4-3,4*1:4-2,4*1:4-1,4*1:4),]
rownames(result) = rep(c('Bias','SD','SE','CP'), 4)
xtable(result, digits=3)

result = rbind(biasFnat_eff, sdFnat_eff, seFnat_eff, cpFnat_eff)
result = result[c(3*1:4-2,3*1:4-1,3*1:4),]
rownames(result) = rep(c('Bias','SD','SE','CP'), 3)
xtable(result, digits=3)

result = rbind(biasFcon_eff, sdFcon_eff, seFcon_eff, cpFcon_eff)
result = result[c(4*1:4-3,4*1:4-2,4*1:4-1,4*1:4),]
rownames(result) = rep(c('Bias','SD','SE','CP'), 4)
xtable(result, digits=3)

result = rbind(biasFcon2_eff, sdFcon2_eff, seFcon2_eff, cpFcon2_eff)
result = result[c(4*1:4-3,4*1:4-2,4*1:4-1,4*1:4),]
rownames(result) = rep(c('Bias','SD','SE','CP'), 4)
xtable(result, digits=3)
