clear all
do twowayreg.ado


capture program drop twowayregwrap

program define twowayregwrap, eclass sortpreserve
version 14 
syntax varlist(numeric ts fv) [if] [in], [,ABSorb(varlist min=2 max=3) DROP GENerate(name)] [, NEWVars(name) REPLACE] [, VCE(name)] [, SAVE rootsave(name) foldersave(string)]
gettoken depvar indepvars : varlist

if ("`ABSorb('varlist')'"=="`absorb('varlist')'" & "`drop'"=="drop"){
    twowayset `absorb', drop
	
	if ("`NEWVars(`name')'"=="`newvars(`name')'" & "`replace'"==""){
		qui{
			des, varlist
		}
		local myvars=r(varlist)
		projvar `depvar' `indepvars', p(`NEWVars' w_) 
		qui{
			des,varlist
		}
		local myvars2=r(varlist)
		local tokeep : list myvars2-myvars  
		twowayreg `tokeep', `vce'
		}
		
		
	else if ("`NEWVars(`name')'"=="" & "`replace'"=="replace"){
		projvar `depvar' `indepvars', replace
		twowayreg `depvar' `indepvars', `vce'
		
	}
}

else if("`ABSorb('varlist')'"=="`absorb('varlist')'" & "`drop'"!="drop"){
	twowayset `absorb', gen(`generate')
	
		if ("`NEWVars(`varname')'"=="`newvars(`varname')'" & "`replace'"==""){
		qui{
			des, varlist
		}
		local myvars=r(varlist)
		projvar `depvar' `indepvars', p(`NEWVars' w_) 
		qui{
			des,varlist
		}
		local myvars2=r(varlist)
		local tokeep : list myvars2-myvars  
		twowayreg `tokeep' if `generate'==1, `vce'
		}
		
		
	else if ("`NEWVars(`name')'"=="" & "`replace'"=="replace"){
		projvar `depvar' `indepvars', replace
		twowayreg `depvar' `indepvars' if `generate'==1, `vce'
	}
	
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