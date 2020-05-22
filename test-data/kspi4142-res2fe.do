cap log close
cap rm kspi4142-res2fe.txt
log using kspi4142-res2fe.txt, text

use kspi4142-res2fe.dta,clear
set seed 4142

// Consecutive integers: useful for generating the dummies below
replace project_n2 = project_n2 - 2
gen long project_id3 = 1
sort project_id2
replace project_id3 = project_id3[_n-1] + 1 - 1*(project_id2[_n] == project_id3[_n-1]) if _n>1

// Data Generation
gen x1=rnormal(0,1)
gen x2=rnormal(0,1)
gen d1=rnormal(0,1)
gen d2=rnormal(0,1)
gen c1=d1[project_n2]
gen c2=d2[project_id3]
drop d1 d2

gen y = x1 + x2 + c1 + c2 + (mod(project_n2,2)+1)*rnormal(0,10)
// Added some heteroskedasticity


// REGHDFE
timer on 1
reghdfe y x1 x2, abs(project_n2 project_id3) vce(robust)
timer off 1
matrix V_reghdfe = e(V)

// RES2FE
timer on 2
twowayset project_n2 project_id3
projvar y x1 x2, p(w_)
reg w_y w_x1 w_x2, robust
dofadj
timer off 2
matrix V_res2fe = e(V)
matrix adj_V_res2fe = V_res2fe * (`r(dofadj)'^2)
matrix list  V_reghdfe 
matrix list adj_V_res2fe
matrix list V_res2fe
timer list
log close


