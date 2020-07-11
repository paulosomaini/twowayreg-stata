clear all
do twowayreg.ado


capture program drop twowayregwrap

program define twowayregwrap, eclass sortpreserve
version 14 
syntax varlist(numeric ts fv) [if] [in], [,ABSorb(varlist min=2 max=3) DROP GENerate(name)] [, NEWVars(name) REPLACE] [, ROBUST VCE VCE_2]
gettoken depvar indepvars : varlist
_fv_check_depvar `depvar'
fvexpand `indepvars' 
marksample touse

if ("`ABSorb('varlist')'"=="`absorb('varlist')'" & "`drop'"=="drop"){
    twowayset `absorb', drop
	
	if ("`NEWVars(`name')'"=="`newvars(`name')'" & "`replace'"==""){
		projvar `depvar' `indepvars', p(`NEWVars' `varname')
		local newdepvar w_dist
		local newindepvar w_mean_np_5 w_x*
		
		if("`robust'"=="robust") {
			twowayreg `newdepvar' `newindepvar', robust
			}
		else if("`vce'"=="vce"){
			twowayreg `w_'`depvar' `w_'`indepvars', vce
		}
		else if("`vce_2'"=="vce_2"){
			twowayreg `w_'`depvar' `w_'`indepvars', vce_2
		}
		else{
			twowayreg `w_'`depvar' `w_'`indepvars'
		}

		}
		
		
	else if ("`NEWVars(`name')'"=="" & "`replace'"=="replace"){
		projvar `depvar' `indepvars', replace
		
		if("`robust'"=="robust") {
		twowayreg `depvar' `indepvars', robust
			}
		else if("`vce'"=="vce"){
			twowayreg `depvar' `indepvars', vce
		}
		else if("`vce_2'"=="vce_2"){
			twowayreg `depvar' `indepvars', vce_2
		}
		else{
			twowayreg `depvar' `indepvars'
		}
	}
}

else if("`ABSorb('varlist')'"=="`absorb('varlist')'" & "`drop'"!="drop"){
	twowayset `absorb', gen(`generate')
	
		if ("`NEWVars(`varname')'"=="`newvars(`varname')'" & "`replace'"==""){
		projvar `depvar' `indepvars', p(`NEWVars' w_)
		
		if("`robust'"=="robust") {
			twowayreg `w_depvar' `w_indepvars' if `generate'==1, robust
			}
		else if("`vce'"=="vce"){
			twowayreg `w_'`depvar' `w_'`indepvars' if `generate'==1, vce
		}
		else if("`vce_2'"=="vce_2"){
			twowayreg `w_'`depvar' `w_'`indepvars' if `generate'==1, vce_2
		}
		else{
			twowayreg `w_'`depvar' `w_'`indepvars' if `generate'==1
		}

		}
		
		
	else if ("`NEWVars(`name')'"=="" & "`replace'"=="replace"){
		projvar `depvar' `indepvars', replace
		
		if("`robust'"=="robust") {
		twowayreg `depvar' `indepvars' if `generate'==1 , robust
			}
		else if("`vce'"=="vce"){
			twowayreg `depvar' `indepvars' if `generate'==1, vce
		}
		else if("`vce_2'"=="vce_2"){
			twowayreg `depvar' `indepvars' if `generate'==1, vce_2
		}
		else{
			twowayreg `depvar' `indepvars' if `generate'==1
		}
	}
	
}

   drop twoWaynewid
	drop twoWaynewt


end 