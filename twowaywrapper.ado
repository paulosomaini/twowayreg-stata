capture program drop twowayregwrap

program define twowayregwrap, eclass sortpreserve
version 14 
syntax varlist(numeric ts fv) [if] [in], [,ABSorb(varlist min=2 max=3) NOPROJ] [, NEWVars(name) REPLACE] [, VCE(name)] [, SAVE rootsave(name) foldersave(string)]
gettoken depvar indepvars : varlist


if ("`ABSorb('varlist')'"=="`absorb('varlist')'" & "`noproj'"==""){
	qui{
	tempvar auxiliar1
	mark `auxiliar1' `if' `in'
	markout `auxiliar1' `varlist'
	}
    twowayset `absorb'
	
	 if ("`NEWVars(`name')'"=="`newvars(`name')'" & "`replace'"==""){
		capture confirm variable `newvars'
		if !_rc {
                 di in red "There is at least one variable with the same prefix chosen, please change the prefix or drop the variable"
				}
        else {
              projvar `depvar' `indepvars', p(`newvars')
			  twowayreg `newvars'* if `auxiliar1'==1, `vce'
			  }
				
	}
		
		
	else if ("`NEWVars(`name')'"=="" & "`replace'"=="replace" ){
		projvar `depvar' `indepvars', replace
		twowayreg `depvar' `indepvars' if `auxiliar1'==1, `vce'
		
	}
}
else if ("`noproj'"=="noproj"){
	qui{
	tempvar auxiliar1
	mark `auxiliar1' `if' `in'
	markout `auxiliar1' `varlist'
	}
	twowayreg `depvar' `indepvars' if `auxiliar1'==1, `vce'
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