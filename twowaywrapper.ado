capture program drop twowayregwrap

program define twowayregwrap, eclass sortpreserve
version 14 
syntax varlist(numeric ts fv) [if] [in], [,ABSorb(varlist min=2 max=3) DROP GENerate(name)] [, NEWVars(name) REPLACE NOPROJ] [, VCE(name)] [, SAVE rootsave(name) foldersave(string)]
gettoken depvar indepvars : varlist


if ("`ABSorb('varlist')'"=="`absorb('varlist')'" & "`drop'"=="drop"){
	twowayset `absorb', drop
	
	 if ("`NEWVars(`name')'"=="`newvars(`name')'" & "`replace'"=="" & "`noproj'"==""){
		capture confirm variable `newvars'
		if !_rc {
                 di in red "There is at least one variable with the same prefix chosen, please change the prefix or drop the variable"
				}
        else {
              projvar `depvar' `indepvars', p(`newvars')
			  twowayreg `newvars'*, `vce'
			  }
				
	}
		
		
	else if ("`NEWVars(`name')'"=="" & "`replace'"=="replace" & "`noproj'"==""){
		projvar `depvar' `indepvars', replace
		twowayreg `depvar' `indepvars', `vce'
		
	}
	
	else if ("`NEWVars(`name')'"=="" & "`replace'"=="" & "`noproj'"=="noproj"){
		twowayreg `depvar' `indepvars', `vce'
	}
drop twoWaynewid 
drop twoWaynewt
}

else if("`ABSorb('varlist')'"=="`absorb('varlist')'" & "`drop'"!="drop"  & "`noproj'"==""){
	twowayset `absorb', gen(`generate')
	
		if ("`NEWVars(`varname')'"=="`newvars(`varname')'" & "`replace'"==""){
		capture confirm variable `newvars'
		if !_rc {
                 di in red "There is at least one variable with the same prefix chosen, please change the prefix or drop the variable"
				}
        else {
              projvar `depvar' `indepvars', p(`newvars')
			  twowayreg `newvars'* if `generate'==1, `vce'
			  }
	}	
		
	else if ("`NEWVars(`name')'"=="" & "`replace'"=="replace" ){
		projvar `depvar' `indepvars', replace
		twowayreg `depvar' `indepvars' if `generate'==1, `vce'
	}

}	

else if ("`ABSorb('varlist')'"=="`absorb('varlist')'" & "`drop'"!="drop" & "`noproj'"=="noproj"){
		twowayreg `depvar' `indepvars' if `generate'==1, `vce'
}	
	


/*
if ("`SAVE'"=="save"){
    twowaysave
}


else if ("`rootsave(`name')'"=="`rootsave(name)'"){
    twowaysave, root(`rootsave')
}

else if ("`foldersave(`string')'"=="`foldersave(`string')'"){
    twowaysave, folder(`foldersave')
}
*/


end 
