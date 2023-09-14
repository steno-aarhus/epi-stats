* Change this to match your folder name on your disk
cd D:\stovring\SDCA\EpiSpace_EpiStats\EffMod

capture log close
log using effmod_analysis.log, replace
clear

use partic_data_dist


** Analysis with logistic regression

* Effect in full cohort
logit nonpartic b0.distance, or

* Stratified analyses
bysort access: logit nonpartic b0.distance, or

* Model with interaction
logit nonpartic b0.distance##b0.access, or
estimates store model1

qui logit nonpartic b0.distance b0.access, or
* Test for interaction
lrtest model1


** Same analysis but with risk ratio instead of odds ratio

* Effect in full cohort
binreg nonpartic b0.distance, rr

* Stratified analyses
bysort access: binreg nonpartic b0.distance, rr

* Model with interaction (note that we need to use ML estimation to use 
* lrtest (likelihood ratio test))
binreg nonpartic b0.distance##b0.access, rr ml
estimates store model2

qui binreg nonpartic b0.distance b0.access, rr ml
* Test for interaction
lrtest model2

log close
