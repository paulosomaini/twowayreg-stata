capture program drop twowayregwrap
capture program drop twowayset
capture mata mata drop sparse()
capture mata mata drop excludemissing()
capture mata mata drop proddiag()
capture mata mata drop diagprod()
capture mata mata drop diagminus()
 capture mata mata drop projDummies()
capture mata mata drop saveMat()
capture mata mata drop readMat()

//Mata programs:


mata:
real matrix sparse(real matrix x)
 {
  real matrix y
  real scalar k
 
  y = J(colmax(x[,1]),colmax(x[,2]),0) 
  for (k=1; k<=rows(x); k++) {
    y[x[k,1],x[k,2]] = y[x[k,1],x[k,2]] + x[k,3]
  }
 
  return(y)
 }
 
  //sparse matrix function ends
  

 // multiplying a diagonal matrix represented by a vector times a matrix.
 // Diag*A multiplies each rows.
 real matrix diagprod(real colvector x, real matrix A)
 {
  real matrix y
  real scalar k
  if(rows(x)<cols(x)) x = x'
 
  y = J(rows(A),cols(A),0)
  for (k=1; k<=rows(x); k++) {
    y[k,] = A[k,] * x[k,1]
  }
 
  return(y)
 }
 
   matrix readMat(string s,string n)
 {
  matrix X
  real scalar fh
fh = fopen(s+"_"+n, "r")
X = fgetmatrix(fh)
fclose(fh)
return(X)
 }
 
 
  
 void saveMat(string s,string n,matrix X)
 {
  real scalar fh
fh = fopen(s + "_" + n, "rw")
fputmatrix(fh, X)
fclose(fh)
 }
 
 
  real matrix proddiag(real matrix A,real colvector x)
 {
  real matrix y
  real scalar k
  if(rows(x)<cols(x)) x = x'
 
  y = J(rows(A),cols(A),0)
  for (k=1; k<=rows(x); k++) {
    y[,k] = A[,k] * x[k,1]
  }
 
  return(y)
 }
 
   real matrix diagminus(real colvector x,real matrix A)
 {
  //real matrix y
  real scalar k
  if(rows(x)<cols(x)) x = x'
 
  //y = -A
  for (k=1; k<=rows(x); k++) {
    A[k,k] = A[k,k] - x[k,1]
  }
 
  return(-A)
 }
 

 void projDummies()
{
real matrix D, DH1, DH, CinvHHDH, AinvDDDH, A, B, C
real colvector DD, HH, invDD, invHH
real scalar N, T
string scalar newid,newt,w,sampleVarName
D=.
newid=st_local("var1")
newt=st_local("var2")
w = st_local("twoway_w")
sampleVarName = st_local("touse_set2")
if (w==""){
D = st_data(.,(newid,newt),sampleVarName)
D = (D,J(rows(D),1,1))
}
else {
D = st_data(.,(newid,newt,w),sampleVarName)
}


DH1=sparse(D)
DD=quadrowsum(DH1)
HH=quadcolsum(DH1)'
HH=HH[1..cols(DH1)-1]
DH=DH1[.,1..cols(DH1)-1]
invDD=DD:^-1 
invHH=HH:^-1

N=colmax(D)[.,1]
T=colmax(D)[.,2]
st_numscalar("e(dimN)",N)
st_numscalar("e(dimT)",T)
st_matrix("e(invDD)",invDD)
st_matrix("e(invHH)",invHH) 

if (N<T)
		{
        
        CinvHHDH=diagprod(invHH,DH')
		A=qrinv(diagminus(DD,CinvHHDH'*DH'))
		//st_matrix("CinvHHDH",CinvHHDH)
        B=-A*CinvHHDH'
		st_matrix("e(CinvHHDH)",CinvHHDH)
		st_matrix("e(A)",A)
		st_matrix("e(B)",B)
		
		
		}
    else
	{
        AinvDDDH=diagprod(invDD,DH)
		C=qrinv(diagminus(HH,AinvDDDH'*DH))
		//st_matrix("AinvDDDH",AinvDDDH)
        B=-AinvDDDH*C
		st_matrix("e(AinvDDDH)",AinvDDDH)
		st_matrix("e(C)",C)
		st_matrix("e(B)",B)

		
    }
 }
 
 end
 
program define nonredundants, eclass sortpreserve
version 11
syntax varlist(min=2 max=2) [if] [in], GENerate(name)

gettoken twoway_id twoway_t: varlist

	*touse_red is created to pass the if and in options of twowayset or twowayload to nonredundants 
	tempvar touse_red
	mark `generate' `if' `in'
	
	tempvar howmany
	count if `generate'== 1
	
	*with this part of the code we are able to discard the redundants observations of analysis 
	while `r(N)' {
			bys `twoway_id': gen `howmany' = _N if `generate'
			replace `generate'= 0 if `howmany' == 1
			drop `howmany'

			bys `twoway_t': gen `howmany' = _N if `generate'
			replace `generate'= 0 if `howmany' == 1
						
			count if `howmany' == 1
			drop `howmany'
			}
end


program define twowayset, eclass sortpreserve
version 11
syntax varlist(min=2 max=3) [if] [in], [GENerate(namelist min=2 max=2) Nogen]
gettoken twoway_id aux: varlist
gettoken twoway_t twoway_w: aux

	qui{
	*if and in options to twowayset
	tempvar touse_set
	mark `touse_set' `if' `in'
	markout `touse_set' `varlist'
	*Discard the observations with negative weights
	if !("`twoway_w'"==""){
	replace `twoway_w' = . if `twoway_w'<=0
	replace `touse_set' = 0 if `twoway_w' == .
	}

	tempvar touse_set2 touse_set3
	nonredundants `twoway_id' `twoway_t' if `touse_set', gen(`touse_set2')

	*touse_set2 is used in ProjDummies() as a marker of nonredundants observations
	*touse_set_3 is created to be used in e(sample)
	gen `touse_set3'=`touse_set2'
	
	ereturn post, esample(`touse_set3')
}

	*option nogen to not generate extra fixed effects, this command is useful if the fixed effects are consecutives
	if ("`nogen'"=="nogen"){
		sort `twoway_id' `twoway_t'
		qui{
			tempvar check1
			gen `check1'=`twoway_id'[_n]-`twoway_id'[_n-1]
			replace `check1'=1 if _n==1
			capture assert `check1'<=1 
			local rc = _rc
		}
		if `rc'{
			di "{err} The fixed effects indices are not consecutive integers. Please, use the option gen to generate consecutive variables." 
			exit `rc'
		}
		
		tempvar var1 var2
		local var1  "`twoway_id'"
		local var2  "`twoway_t'"
	}
	else{
		*if the fixed effects are not consecutive the user has to create new variables to be used to create D matrix
		qui{
		gettoken var1 var2: generate
		egen `var1'= group(`twoway_id') if `touse_set2'==1
		egen `var2'= group(`twoway_t') if `touse_set2'==1
	
		}		
		}
	ereturn local absorb "`var1' `var2' `twoway_w'"
	mata projDummies()
	drop `touse_set'
	
	


end


capture program drop projvar
capture mata mata drop projVar()


mata
void projVar()
{
	real matrix V, varIn, D,aux,delta,tau,varOut,A,B,CinvHHDH,AinvDDDH,C
	real colvector invHH,invDD,Dy,Ty
	real scalar N,T
	string scalar newid, newt, currvar,newvar,sampleVarName,w,linear_index,var1,var2
	currvar = st_local("currvar")
	newvar = st_local("newvar")
	var1 = st_local("var1")
	var2 = st_local("var2")
	N=st_numscalar("e(dimN)")
	T=st_numscalar("e(dimT)")
	
	w=st_local("w")
	sampleVarName = st_local("touse_proj")
	linear_index = st_local("linear_index")
	V = st_data(.,(var1, var2,currvar),sampleVarName)
	varIn=V[.,3]
	
	if (w==""){
	D = st_data(.,(var1, var2),sampleVarName)
	D = (D,J(rows(D),1,1))
	}
	else {
	D = st_data(.,(var1, var2,w),sampleVarName)
	}
	
	V[.,3]=V[.,3]:*D[.,3]
	aux=sparse(V)
	//printf("3")
	Dy=rowsum(aux)
	Dy=Dy
	Ty=colsum(aux)
	Ty=Ty[1,1..cols(aux)-1]'
	B=st_matrix("e(B)")
			

	 if (N<T)
			{
			
			A=st_matrix("e(A)")
			invHH=st_matrix("e(invHH)")
			CinvHHDH=st_matrix("e(CinvHHDH)")
			delta=A*Dy+B*Ty
			tau=B'*(Dy-CinvHHDH'*Ty)+(invHH:*Ty) \0
			}
		else
		{
			C=st_matrix("e(C)")
			invDD=st_matrix("e(invDD)")
			AinvDDDH=st_matrix("e(AinvDDDH)")
			delta=(invDD:*Dy)+B*(Ty-AinvDDDH'*Dy)
			tau=B'*Dy+C*Ty \0 
		}

	varOut=(varIn-delta[V[.,1]]-tau[V[.,2]]):*sqrt(D[.,3])
	st_store(st_data(.,linear_index,sampleVarName), newvar, varOut)

}
end



program define projvar, nclass
version 11
syntax varlist, [Prefix(name)] [REPLACE]
	
	*in e(absorb) there is the fixed effects that we use to generate the new matrix V
	local absorb = "`e(absorb)'"
	gettoken var1 aux : absorb
	gettoken var2 twoway_w : aux
	local w `twoway_w'
	
	*set the vars to be projected
	gettoken depvar indepvars : varlist
    _fv_check_depvar `depvar'
    fvexpand `indepvars' 
	
	
	qui{
		
	tempvar touse_proj touse_check linear_index
	*set the sample that is projected
	gen byte `touse_proj'= e(sample)
	gen byte `touse_check'  =  `touse_proj'
	*touse_check is a marker of missing values of the vars projected
	markout `touse_check'  `varlist'
	*check if there is a missing in some of the variables
	capture assert  `touse_proj' ==  `touse_check' 
	local rc = _rc
	}
	if `rc' {
		  di "{err} Some of the included variables have missing values."
	   exit `rc'
	}
	drop `touse_check'
	
	
	*variable created to store only the observations that are non-missings and in that way the arrays are conformables
	gen `linear_index' = _n	
	
	*we create the new variables, if a variable already exists then it is ignored and not projected
	foreach currvar of varlist `varlist' {
		local newvar="`prefix'`currvar'"
		if ("`replace'" != "") {
		local newvar="`currvar'"
			mata projVar()
		}
		else {
				capture confirm variable `prefix'`currvar'
		if !_rc { 
			}
			else{
				qui gen `newvar'=.
				mata projVar()
			}
		}

	
	}

drop `linear_index'
*save in scalars and arrays the macros.
scalar dimN= e(dimN)
scalar dimT= e(dimT)
matrix invDD=e(invDD)
matrix invHH=e(invHH)
if (dimN<dimT){
	matrix CinvHHDH=e(CinvHHDH)
	matrix A= e(A)
	matrix B=e(B)
}
else {
	matrix AinvDDDH=e(AinvDDDH)
	matrix C= e(C)
	matrix B=e(B)
}


end

capture program drop twowaysave
capture mata mata drop matasave()

mata
void matasave()
{
real matrix D, CinvHHDH, AinvDDDH, A, B, C
real colvector  invDD, invHH, var1,var2
real scalar N, T
string scalar newid,newt, w,sampleVarName, root


root =st_local("using")
newid=st_local("var_1")
st_local("newid",newid)
newt=st_local("var_2")
st_local("newt",newt)
w=st_local("twoway_w")
st_local("w",w)
N=st_numscalar("dimN")
T=st_numscalar("dimT")
invDD=st_matrix("invDD")
invHH=st_matrix("invHH") 
saveMat(root,"twoWayVar1", newid)
saveMat(root,"twoWayVar2", newt)
saveMat(root,"twoWayW",w)
saveMat(root,"twoWayN1", N)
saveMat(root,"twoWayN2", T)
saveMat(root,"twoWayInvDD", invDD)
saveMat(root,"twoWayInvHH", invHH)


 if (N<T)
		{
        CinvHHDH=st_matrix("CinvHHDH")
		A=st_matrix("A")
		B=st_matrix("B")

		saveMat(root,"twoWayCinvHHDH", CinvHHDH)
		saveMat(root,"twoWayA", A)
		saveMat(root,"twoWayB", B)
		
		
		}
    else
	{
        AinvDDDH=st_matrix("AinvDDDH")
		C=st_matrix("C")
		B=st_matrix("B")
		saveMat(root,"twoWayAinvDDDH", AinvDDDH)
		saveMat(root,"twoWayC", C)
		saveMat(root,"twoWayB", B)
		
    }


}

end


program define twowaysave, eclass
version 11
syntax [using/]

local absorb = "`e(absorb)'"
gettoken var1 aux: absorb
gettoken var2 w: aux
*obtain the name of the fixed effects
local var_1 `var1'
local var_2 `var2'
local twoway_w `w'

qui{
tempvar touse_save
gen byte `touse_save'= e(sample)
}

*create the scalars that are store in e()
scalar dimN= e(dimN)
scalar dimT=e(dimT)
matrix invDD=e(invDD)
matrix invHH=e(invHH)

if (dimN<dimT){
	matrix CinvHHDH=e(CinvHHDH)
	matrix A= e(A)
	matrix B=e(B)
}
else {
	matrix AinvDDDH=e(AinvDDDH)
	matrix C= e(C)
	matrix B=e(B)
}

*create the new list of macros
ereturn clear
ereturn post, esample(`touse_save')
mata matasave()
ereturn local absorb "`newid' `newt' `w'"
ereturn scalar dimN= dimN
ereturn scalar dimT= dimT
ereturn matrix invDD= invDD
ereturn matrix invHH= invHH
if (dimN<dimT){
	ereturn matrix CinvHHDH=CinvHHDH
	ereturn matrix A= A
	ereturn matrix B= B
}
else {
	ereturn matrix AinvDDDH=AinvDDDH
	ereturn matrix C= C
	ereturn matrix B= B
}


end

capture program drop twowayload
capture mata mata drop mataload()

mata
void mataload()
{
real matrix D, CinvHHDH, AinvDDDH, A, B, C
real colvector  invDD, invHH
real scalar N, T
string scalar newid,newt, w,sampleVarName, root

root =st_local("using")
w = readMat(root,"twoWayW")
newid= readMat(root,"twoWayVar1")
newt= readMat(root,"twoWayVar2")

st_local("newid",newid)
st_local("newt",newt)
st_local("w",w)

N=readMat(root,"twoWayN1")
T=readMat(root,"twoWayN2")
invDD=readMat(root,"twoWayInvDD")
invHH=readMat(root,"twoWayInvHH")
st_numscalar("e(dimN)",N)
st_numscalar("e(dimT)",T)
st_matrix("e(invDD)",invDD)
st_matrix("e(invHH)",invHH) 

 if (N<T)
		{
        CinvHHDH=readMat(root,"twoWayCinvHHDH")
		A=readMat(root,"twoWayA")
		B=readMat(root,"twoWayB")
		st_matrix("e(CinvHHDH)",CinvHHDH)
		st_matrix("e(A)",A)
		st_matrix("e(B)",B)
	
		}
    else
	{
        AinvDDDH=readMat(root,"twoWayAinvDDDH")
		C=readMat(root,"twoWayC")
		B=readMat(root,"twoWayB")
		st_matrix("e(AinvDDDH)",AinvDDDH)
		st_matrix("e(C)",C)
		st_matrix("e(B)",B)
		
    }


}

end


program define twowayload, eclass
version 11
syntax [using/] [if] [in]
mata mataload()
*tokenize the names of the fixed effects
gettoken twoway_id: newid
gettoken twoway_t: newt
gettoken twoway_w: w
	qui{
	tempvar touse_set
	mark `touse_set' `if' `in'
	*Discard the observations with negative weights
	if !("`twoway_w'"==""){
	replace `twoway_w' = . if `twoway_w'<=0
	replace `touse_set' = 0 if `twoway_w' == .
	}

	tempvar touse_set2
	nonredundants `twoway_id' `twoway_t' if `touse_set', gen(`touse_set2')
}
*save the arrays, scalars obtained in ProjDummies_load()
scalar dimN= e(dimN)
scalar dimT=e(dimT)
matrix invDD=e(invDD)
matrix invHH=e(invHH)

if (dimN<dimT){
	matrix CinvHHDH=e(CinvHHDH)
	matrix A= e(A)
	matrix B=e(B)
}
else {
	matrix AinvDDDH=e(AinvDDDH)
	matrix C= e(C)
	matrix B=e(B)
}
ereturn post, esample(`touse_set2')
ereturn local absorb "`newid' `newt' `w'"
*store the arrays and scalars in e()
ereturn scalar dimN= dimN
ereturn scalar dimT= dimT
ereturn matrix invDD= invDD
ereturn matrix invHH= invHH
if (dimN<dimT){
	ereturn matrix CinvHHDH=CinvHHDH
	ereturn matrix A= A
	ereturn matrix B= B
}
else {
	ereturn matrix AinvDDDH=AinvDDDH
	ereturn matrix C= C
	ereturn matrix B= B
}


end

capture program drop twowayreg 

program define twowayreg, eclass sortpreserve
    version 11
 
    syntax varlist(numeric ts fv),[,ROBUST VCE(namelist) statadof] 
	local absorb = "`e(absorb)'"
    gettoken depvar indepvars : varlist
    _fv_check_depvar `depvar'
    fvexpand `indepvars'

	qui{
	tempvar touse_reg
	gen byte `touse_reg'= e(sample)
	}
	scalar dimN= e(dimN)
	scalar dimT= e(dimT)
	matrix invDD=e(invDD)
	matrix invHH=e(invHH)
	if (dimN<dimT){
		matrix CinvHHDH=e(CinvHHDH)
		matrix A= e(A)
		matrix B=e(B)
	}
	else {
		matrix AinvDDDH=e(AinvDDDH)
		matrix C= e(C)
		matrix B=e(B)
	}
	
	
   if ("`robust'"=="robust"){
   *standard errors  proposed by Arellano (1987) robust to heteroscedasticity and serial correlation    
   	qui{
  	regress `depvar' `indepvars' , noc robust
	scalar vadj = e(df_r)/(e(df_r)- dimN - dimT)
	scalar df_r1= e(df_r) - dimN - dimT
    }
 }
    else if ("`vce(namelist)'"== "`vce(namelist)'" & "`statadof'"== ""){
	*standard errors robust to heteroscedasticity but assumes no correlation within group or serial correlation.
   qui{
  	regress `depvar' `indepvars'  , noc vce(`vce')
	qui{
		scalar df_r= e(N)-e(df_m)-1
	}
	scalar df_r1= e(df_r)
	scalar vadj = df_r/(df_r- dimN - dimT)
	    }
 }

  else if ("`vce(namelist)'"== "`vce(namelist)'" & "`statadof'"== "statadof"){
	*Arellano standard errors with a degree of freedom correction performed by Stata xtreg, fe.
   qui{
  	regress `depvar' `indepvars' , noc vce(`vce')
	scalar N_1=e(N)
	scalar df_m= e(df_m)
	qui{
		scalar df_r= e(N)-e(df_m)-1
	}
	scalar df_r1= e(df_r)
	scalar vadj = (N_1-1)*(df_r/(df_r - 1))/(N_1 - df_m - 1)
	    }

 
}
 
  else{
     *standard errors assuming homoscedasticity and no within  group correlation or serial correlation
  qui{
  	regress `depvar' `indepvars'  , noc
	scalar df_r1= e(df_r)-dimN-dimT
	scalar vadj = e(df_r)/(e(df_r)- dimN - dimT)
	}
  }
	mat b=e(b)
	scalar N_1=e(N)
	scalar R2= e(r2)
	scalar F= e(F)
	scalar df_m= e(df_m)
	scalar rtms= e(rmse)
	matrix V = vadj*e(V)

  *table of regression with standar errors with dof correction
  eret post b V, esample(`touse_reg')
  *macros 
  ereturn local absorb "`absorb'"
  ereturn scalar N_1= N_1
  ereturn scalar R2= R2
  ereturn scalar F= F
  ereturn scalar df_m=df_m
  ereturn scalar df_r1= df_r1
  ereturn scalar rtms= rtms
  ereturn scalar dimN= dimN
  ereturn scalar dimT= dimT
  ereturn matrix invDD= invDD
  ereturn matrix invHH= invHH
  if (dimN<dimT){
	ereturn matrix CinvHHDH=CinvHHDH
	ereturn matrix A= A
	ereturn matrix B= B
 }
else {
	ereturn matrix AinvDDDH=AinvDDDH
	ereturn matrix C= C
	ereturn matrix B= B
}
  *table display
  display _newline "Two-Way Regression" _col(45) "Number of obs" _col(60)"=" _col(65) N_1
  display _col(45) "F(" df_m "," df_r1 ")"  _col(60)"="  _col(65) F
  display _col(45) "R-squared" _col(60)"="  _col(65) R2
  display _col(45) "Root MSE " _col(60)"="  _col(65) rtms
  eret display

  

  
end 
 


program define twowayregwrap, eclass sortpreserve
version 11
syntax varlist(numeric ts fv) [if] [in], [,ABSorb(varlist min=2 max=3) GENerate(namelist) NOGEN NOPROJ] [, NEWVars(name) REPLACE] [, VCE(namelist) statadof]
gettoken depvar indepvars : varlist


	qui{
	*tempvar to use if and in and delete missings observation of analysis
	tempvar touse_wrap
	mark `touse_wrap' `if' `in'
	markout `touse_wrap' `varlist'
	}
if ("`noproj'"=="" & "`nogen'"=="nogen"){
	*make the whole regression without creating new fixed effects
	twowayset `absorb' if `touse_wrap', nogen
	gettoken twoway_id twoway_t : absorb
	if ("`NEWVars(`name')'"=="`newvars(`name')'" & "`replace'"==""){
		capture confirm variable `newvars'
		if !_rc {
                 di "{err} There is at least one variable with the same prefix chosen, please change the prefix or drop the variable"
				}
        else {
              projvar `depvar' `indepvars', p(`newvars')
			  twowayreg `newvars'* , vce(`vce') `statadof'
			  }
				
	}
		
		
	else if ("`NEWVars(`name')'"=="" & "`replace'"=="replace" ){
		projvar `depvar' `indepvars' , replace
		twowayreg `depvar' `indepvars' , vce(`vce') `statadof'
		
	}
}	
if ("`noproj'"=="" & "`nogen'"==""){
	*make the whole regression creating new fixed effects
	twowayset `absorb' if `touse_wrap', gen(`generate')
	gettoken twoway_new_id twoway_new_t : generate

	 if ("`NEWVars(`name')'"=="`newvars(`name')'" & "`replace'"==""){
		capture confirm variable `newvars'
		if !_rc {
                 di "{err} There is at least one variable with the same prefix chosen, please change the prefix or drop the variable"
				}
        else {
              projvar `depvar' `indepvars', p(`newvars')
			  twowayreg `newvars'* , vce(`vce') `statadof'
			  }
				
	}
		
		
	else if ("`NEWVars(`name')'"=="" & "`replace'"=="replace" ){
		projvar `depvar' `indepvars', replace
		twowayreg `depvar' `indepvars' if `touse_wrap', vce(`vce') `statadof'
		
	}
}	

else if ("`noproj'"=="noproj"){
	*option just to make the regression without setting the fixed effects or projecting varlist
	gettoken twoway_id aux: absorb
	gettoken twoway_t w: aux

	
	twowayreg `depvar' `indepvars', vce(`vce') `statadof'
}




end

** End of twowayreg.ado
