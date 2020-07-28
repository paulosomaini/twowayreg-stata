cd "d:\Users\user\Documents\GitHub\twowayreg-stata"

clear all
do twowayreg.ado
use fortest.dta

twowayset ExperimentID dayHour if lnKWH!=., gen(hola1 hola2)
twowaysave 

projvar lnKWH tier1Post,p(w_)
twowayreg w_lnKWH w_tier1Post,  robust 
ereturn list
projvar tier2Post, p(w_)
twowayreg w_lnKWH w_tier1Post w_tier2Post,  robust abs(hola1 hola2)

ereturn list

clear all
do twowayreg.ado
use fortest.dta
twowayload ExperimentID dayHour if lnKWH!=., gen(hola1 hola2) 
ereturn list
projvar lnKWH tier1Post tier2Post,p(w_)
twowayreg w_lnKWH w_tier1Post w_tier2Post 

twowayregwrap lnKWH tier1Post tier2Post, absorb(ExperimentID dayHour) gen(hola1 hola2) newv(w_) vce(cluster hola1) 
twowayregwrap w_lnKWH w_tier1Post w_tier2Post, absorb(ExperimentID dayHour) noproj vce(cluster hola1) statadof
twowayregwrap w_lnKWH w_tier1Post w_tier2Post,noproj vce(cluster hola1) statadof




