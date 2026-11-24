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
To compute the efficient influence function and estimate counterfactual CIF; the mediator density ratio is estimated by Riesz regression<br>
The main output consists of the following:<br>
time: time points when the counterfactual CIF is evaluated<br>
F_reg: regression (plug-in) estimate of counterfactual CIF<br>
F_eff: efficient influence function (EIF)-based estimate of counterfactual CIF<br>
se_eff: pointwise standard error of the estimated EIF-based estimate of counterfactual CIF
##### est_decompose(fit)
To decompose the total effect into mediation, direct, and indirect effects; the input is an object from fit_eif_cf
##### est_con(fit, a=0, m_seq=NULL, h2=NULL, time=NULL, h=0.5, monotone=TRUE)
To estimate the kernel-smoothed controlled counterfactual CIFs, setting the mediator at m_seq and the competing event cumulative hazard at h2; the bandwidth of the kernel that smooths the controlled CIF is h

### Files
#### Simulation
##### main.R
Main functions
##### simulation.R
Code for a single run of simulation
##### result.R
Code for summarizing parallel replicates
#### Applications (COV2008; AMP)
##### app1_sens_boot.R
Sensitivity analysis for the COV2008 data
##### app2_build_pt80.R
To generate predicted ID80 in the AMP dataset
##### application_analysis.R
Main code for two applications (COV2008, AMP)
