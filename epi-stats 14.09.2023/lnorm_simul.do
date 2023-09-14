* do D:\stovring\SDCA\EpiSpace_EpiStats\Transformations\lnorm_simul.do

cd D:\stovring\SDCA\EpiSpace_EpiStats\Transformations\

set seed 18956439
capture log close
log using lnorm_simul.log, replace

clear

**********************************************************************
* Creating a simulated dataset to analyze

set obs 342
gen bmi20 = (runiform() - .5) * 16 + 27 - 20
la var bmi20 "BMI - 20"
gen group = runiform() < .5

gen hba1c = exp(.07 * bmi + .2 * (1 - group) + rnormal(.9, 0.2)) 

sa lnorm_data_ex, replace

gen ln_hba1c = ln(hba1c)


**********************************************************************
* Analysis of dataset

regress ln_hba1c b0.group bmi

regress ln_hba1c b0.group bmi, eform("MedRat")

predict pred_lnhba1c, xb
predict resid_lnhba1c, resid

qnorm resid_lnhba1c
graph export qq_lnhba1c.png, replace width(1400)

separate pred_lnhba1c, by(group)
twoway (scatter ln_hba1c bmi20, mlabel(group) mlabp(0) msy(i)) ///
	(line pred_lnhba1c0 pred_lnhba1c1 bmi20, c(L L))
graph export pred_lnhba1c.png, replace width(1400)
	
gen pred_hba1c = exp(pred_lnhba1c)
separate pred_hba1c, by(group)
sort group bmi20
twoway (scatter hba1c bmi20, mlabel(group) mlabp(0) msy(i)) ///
	(line pred_hba1c0 pred_hba1c1 bmi20, c(L L)), ylab(0(2)12)
graph export pred_hba1c.png, replace width(1400)

log close