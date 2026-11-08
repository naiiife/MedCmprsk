# MedCmprsk
Causal mediation analysis for a time-to-event endpoint in the presence of competing events

### Method
#### Data structure
A: treatment<br>
X: covariates<br>
M: mediator (NA = missing)<br>
Time: failure time<br>
status: event indicator (1 for primary event, 2 for competing event, 0 for censoring)

#### Implementation
##### fit_eif_cf(A, Time, status, M, X, method='riesz', event=1, monotone=TRUE)
To compute the efficient influence function and estimate counterfactual CIF<br>
time: time points when the counterfactual CIF is evaluated<br>
F_reg: regression (plug-in) estimate of counterfactual CIF<br>
F_eff: efficient influence function (EIF)-based estimate of counterfactual CIF<br>
se_eff: pointwise standard error of the estimated EIF-based estimate of counterfactual CIF
##### est_decompose(fit)
To decompose the total effect into mediation, direct, and indirect effects
##### est_con(fit, a=0, m_seq=NULL, h2=NULL, time=NULL, h=0.5, monotone=TRUE)
To estimate the controlled counterfactual CIFs, setting the mediator at m_seq and the competing event cumulative hazard at h2

### Files
##### main.R
Main functions
##### simulation.R
Code for a single run of simulation
