* do D:\stovring\SDCA\EpiSpace_EpiStats\Analysis_of_Change\bl_change_ancova.do

cd D:\stovring\SDCA\EpiSpace_EpiStats\Analysis_of_Change

capture log close
log using bl_change_ancova.log, replace

clear
set seed 5158887


*****************************************************************************
* Values controlling the setting
scalar dx = - 0.0    // baseline average difference
scalar eff_size = .7 // Effect size at follow-up

mat m1 = (25, 25)     // BMI at baseline and follow-up (control group)
mat covmat = (5, 4.5 \ 4.5, 5) // Covariance BL vs FU (high correlation)


*****************************************************************************
* Creating simulated dataset
drawnorm bl_x1 fu_y1, n(300) means(m1) cov(covmat)

mat m2 = (25 + dx, 25 + dx + eff_size)
drawnorm bl_x2 fu_y2, n(300) means(m2) cov(covmat)

gen id = _n
reshape long bl_x fu_y, i(id) j(group)


*****************************************************************************
* Assessment of baseline balance
ttest bl_x, by(group)


*****************************************************************************
* SACS
gen diff = fu_y - bl_x
ttest diff, by(group)


*****************************************************************************
* ANCOVA
regress fu_y b1.group c.bl_x

log close
