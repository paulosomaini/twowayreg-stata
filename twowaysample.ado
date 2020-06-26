capture program drop twowaysample
// Creates a variable with non-redundant observations to include in the 2-way fe regression
// Redundant observations cause collinearities between the two sets of fixed effects.


program define twowaysample, sortpreserve
version 11
syntax varlist(min=2 max=3) [if] [in], Generate(name) [Replace]

gettoken iid aux: varlist
gettoken tid w: aux


qui{
if ("`replace'"=="replace") {
	cap drop `generate'
}

if !("`w'"==""){
	replace `w' = . if `w'<=0
}


mark `generate' `if' `in'
markout `generate' `varlist'

tempvar howmany
count if `generate' == 1



while `r(N)' {
	bys `iid': gen `howmany' = _N if `generate'
	replace `generate' = 0 if `howmany' == 1
	drop `howmany'

	bys `tid': gen `howmany' = _N if `generate'
	replace `generate' = 0 if `howmany' == 1
	
	count if `howmany' == 1
	drop `howmany'
}
}

end
