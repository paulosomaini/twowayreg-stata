clear all
do twowayreg.ado

*** 0) Preliminaries

forvalues lo = 3/3 {
di `lo'
forvalues wo = 2/2 {
di `wo'
foreach vars of numlist 2 10 {

di `vars'

loc long = 10^`lo'
loc wide = 10^`wo'
*loc vars = 2
loc lout = 0.1
loc reps = 1

loc toto = `long'*`wide'
set more off

forvalues rep = 1/`reps' {


*** 1) Generate Data
drop _all
set obs `toto'
** Variables
forvalues var = 1/`vars' {
	gen x`var'= rnormal(0)
	}
** Fixed Effects
* Indicators
gen hhid = floor((_n-1)/`wide')
gen ttid = _n-1-hhid*`wide'
** Drop a fraction of observations;
gen out= uniform()
sort out
drop if _n<`lout'*`toto'
* Effects
gen hhef = rnormal(0)
gen ttef = rnormal(0)
bysort hhid: replace hhef = hhef[1]
gen hid = 1
replace hid = hid[_n-1] + 1*(hhid[_n-1]~=hhid[_n]) if _n>1
bysort ttid: replace ttef = ttef[1]
gen tid = 1
replace tid = tid[_n-1] + 1*(ttid[_n-1]~=ttid[_n]) if _n>1


** Dependent Variable
gen y = hhef + ttef + rnormal(0)
forvalues var = 1/`vars' {
	replace y= y + x`var'
	}

*** 2) Run Our procedure
twowayset hid tid, gen(sample)
projvar y x*, p(w_)
twowayreg w_y w_x* if sample==1

drop w_*
drop sample
twowaysave
}
}
}
}

save Example2.dta,replace 

clear all
do twowayreg.ado
use Example2.dta
*** 2) Run Our procedure
twowayload hid tid, gen(sample)
projvar y x*, p(w_)
twowayreg w_y w_x* if sample==1, robust

drop w_*

