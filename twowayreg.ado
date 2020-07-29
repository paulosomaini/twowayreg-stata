capture program drop twowayregwrap
capture program drop twowayprep
capture program drop twowayset
capture mata mata drop dofadj()
capture mata mata drop dofadj_l()
capture mata mata drop sparse()
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
 
  real matrix dofadj(real dof)
 {
  real scalar adj, T, N
  N=st_numscalar("N")
  T=st_numscalar("T")
 adj = sqrt(dof/(dof - N - T + 1))
  return(adj)
 }

 
 real matrix dofadj_l(string root, real dof)
 {
  real scalar adj, T, N
  N=readMat(root,"twoWayN1")
  T=readMat(root,"twoWayN2")
  adj = sqrt(dof/(dof - N - T + 1))
  return(adj)
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
sampleVarName = st_local("touse_set")
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
 
 

program define twowayset, eclass sortpreserve
version 11
syntax varlist(min=2 max=3) [if] [in], [GENerate(namelist min=2 max=2) Nogen]
gettoken twoway_id aux: varlist
gettoken twoway_t twoway_w: aux
//summ `varlist'
// I need to make it robust to non 1,2,3... ids.



	if !("`w'"==""){
	replace `w' = . if `w'<=0
	}
	
	qui{
	tempvar touse_set
	mark `touse_set' `if' `in'


	tempvar howmany
	count if `touse_set'== 1

	while `r(N)' {
			bys `twoway_id': gen `howmany' = _N if `touse_set'
			replace `touse_set'= 0 if `howmany' == 1
			drop `howmany'

			bys `twoway_t': gen `howmany' = _N if `touse_set'
			replace `touse_set'= 0 if `howmany' == 1
						
			count if `howmany' == 1
			drop `howmany'
			}
	tempvar touse_set2
	gen `touse_set2'= `touse_set'
	}	
	
	ereturn post, esample(`touse_set2')
	
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
			di "{err} The fixed effects are not consecutive, please use the option gen to generate consecutive variables." 
			exit `rc'
		} 
		tempvar var1 var2
		gen `var1'= `twoway_id'
		gen `var2'= `twoway_t'

		ereturn local absorb "`twoway_id' `twoway_t'"

	}
	
	else{
		gettoken var1 var2: generate
		egen `var1'= group(`twoway_id')
		egen `var2'= group(`twoway_t')
		ereturn local absorb "`var1' `var2'"
		}
		
	mata projDummies()
	drop `touse_set'
	
	
	//di in gr "Checkpoint 1"
//ret li
//di in gr "Checkpoint 2"
	scalar twoWaynewid= `var1'
	scalar twoWaynewt= `var2'
	scalar twoWayw="`twoway_w'"
	scalar twoWayif="`if'"
	scalar twoWayin="`in'"

//return post r(B), esample(`twoway_sample') 
//obs(`nobs') dof(`dof')


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
	
	w=st_strscalar("twoWayw")
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
	
	//rows(Ty)
    //cols(Ty)
	//rows(Dy)
	//cols(Dy)
			

	 if (N<T)
			{
			
			A=st_matrix("e(A)")
			invHH=st_matrix("e(invHH)")
			CinvHHDH=st_matrix("e(CinvHHDH)")
			//printf("b")
			delta=A*Dy+B*Ty
			tau=B'*(Dy-CinvHHDH'*Ty)+(invHH:*Ty) \0
			}
		else
		{
			//printf("1")
			C=st_matrix("e(C)")
			invDD=st_matrix("e(invDD)")
			AinvDDDH=st_matrix("e(AinvDDDH)")
			delta=(invDD:*Dy)+B*(Ty-AinvDDDH'*Dy)
			tau=B'*Dy+C*Ty \0 
			//printf("c")
		}

	//how to index
	//varout=(var-delta(struc.hhid)-tau(struc.tid')).*sqrt(struc.w);
	varOut=(varIn-delta[V[.,1]]-tau[V[.,2]]):*sqrt(D[.,3])
	//printf("4")
	//st_matrix("DD2",B)
	st_store(st_data(.,linear_index,sampleVarName), newvar, varOut)
	//printf("5")
}
end



program define projvar, nclass
version 11
syntax varlist, [Prefix(name)] [REPLACE]
	
	local absorb = "`e(absorb)'"
	gettoken var1 var2 : absorb
	
	
	gettoken depvar indepvars : varlist
    _fv_check_depvar `depvar'
    fvexpand `indepvars' 
	
	
	qui{
	tempvar touse_proj touse_check linear_index
	gen byte `touse_proj'= e(sample)
	gen byte `touse_check'  =  `touse_proj'
	markout `touse_check'  `varlist'   
	capture assert  `touse_proj' ==  `touse_check' 
	local rc = _rc
	}
	if `rc' {
		  di "{err} Some of the included variables have missing values."
	   exit `rc'
	}
	drop `touse_check'
	
	
	
	gen `linear_index' = _n	
	
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
capture mata mata drop projDummies_save()

capture program drop twowaysave
capture mata mata drop projDummies_save()

mata
void projDummies_save()
{
real matrix D, CinvHHDH, AinvDDDH, A, B, C
real colvector  invDD, invHH, var1,var2
real scalar N, T
string scalar twoWaynewid,twoWaynewt, w,sampleVarName, root


root =st_local("using")
newid=st_local("var_1")
st_local("newid",newid)
newt=st_local("var_2")
st_local("newt",newt)
N=st_numscalar("dimN")
T=st_numscalar("dimT")
invDD=st_matrix("invDD")
invHH=st_matrix("invHH") 
saveMat(root,"twoWayVar1", newid)
saveMat(root,"twoWayVar2", newt)
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
gettoken var1 var2: absorb
local var_1 `var1'
local var_2 `var2'

qui{
tempvar touse_save
gen byte `touse_save'= e(sample)
}
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

ereturn clear
ereturn post, esample(`touse_save')
mata projDummies_save()
ereturn local absorb "`newid' `newt'"
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
capture mata mata drop projDummies_load()

mata
void projDummies_load()
{
real matrix D, CinvHHDH, AinvDDDH, A, B, C
real colvector  invDD, invHH
real scalar N, T
string scalar twoWaynewid,twoWaynewt, w,sampleVarName, root

root =st_local("using")
w = st_local("twoway_w")
newid= readMat(root,"twoWayVar1")
newt= readMat(root,"twoWayVar2")

st_local("newid",newid)
st_local("newt",newt)

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
mata projDummies_load()
gettoken twoway_id: newid
gettoken twoway_t: newt
	qui{
	tempvar touse_set
	mark `touse_set' `if' `in'


	tempvar howmany
	count if `touse_set'== 1

	while `r(N)' {
			bys `twoway_id': gen `howmany' = _N if `touse_set'
			replace `touse_set'= 0 if `howmany' == 1
			drop `howmany'

			bys `twoway_t': gen `howmany' = _N if `touse_set'
			replace `touse_set'= 0 if `howmany' == 1
						
			count if `howmany' == 1
			drop `howmany'
			}
	tempvar touse_set2
	gen `touse_set2'= `touse_set'
	}	
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
ereturn local absorb "`newid' `newt'"
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
scalar twoWayw="`twoway_w'"
scalar twoWayif="`if'"
scalar twoWayin="`in'"

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


  eret post b V, esample(`touse_reg')
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
  
  display _newline "Two-Way Regression" _col(45) "Number of obs" _col(60)"=" _col(65) N_1
  display _col(45) "F(" df_m "," df_r1 ")"  _col(60)"="  _col(65) F
  display _col(45) "R-squared" _col(60)"="  _col(65) R2
  display _col(45) "Root MSE " _col(60)"="  _col(65) rtms
  eret display

  

  
end 
 


capture program drop dofadj
program define dofadj, rclass
version 11
syntax ,[Root(name)]
local dof = `e(N)'-`e(df_m)'-1
mata dofadj(`dof')
mata st_local("dofadj", strofreal(dofadj(`dof')))
return scalar dofadj = `dofadj'
end

capture program drop dofadj_l

program define dofadj_l, rclass
version 11
syntax ,[Root(name)]
local dof = `e(N)'-`e(df_m)'-1
mata dofadj_l("`root'",`dof')
mata st_local("dofadj", strofreal(dofadj_l("`root'",`dof')))
return scalar dofadj = `dofadj'
end

program define twowayprep, eclass sortpreserve
version 11
syntax varlist(numeric ts fv) [if] [in], [,ABSorb(varlist min=2 max=3)] [,NEWVars(name) REPLACE]
gettoken depvar indepvars : varlist


	gettoken twoway_id aux: absorb
	gettoken twoway_t w: aux
	qui{
	tempvar touse_prep
	mark `touse_prep' `if' `in'
	markout `touse_prep' `varlist'

	tempvar howmany
	count if `touse_prep'== 1

	while `r(N)' {
			bys `twoway_id': gen `howmany' = _N if `touse_prep'
			replace `touse_prep'= 0 if `howmany' == 1
			drop `howmany'

			bys `twoway_t': gen `howmany' = _N if `touse_prep'
			replace `touse_prep'= 0 if `howmany' == 1
						
			count if `howmany' == 1
			drop `howmany'
			}
	}	

	twowayset `absorb' if `touse_prep'
	if ("`NEWVars(`name')'"=="`newvars(`name')'" & "`replace'"==""){
	          projvar `depvar' `indepvars' if `touse_prep', p(`newvars')
	}
	else if ("`NEWVars(`name')'"=="" & "`replace'"=="replace" ){
		projvar `depvar' `indepvars' if `touse_prep', replace
	}

drop twoWaynewid twoWaynewt

end 

program define twowayregwrap, eclass sortpreserve
version 11
syntax varlist(numeric ts fv) [if] [in], [,ABSorb(varlist min=2 max=3) GENerate(namelist) NOGEN NOPROJ] [, NEWVars(name) REPLACE] [, VCE(namelist) statadof] [, SAVE rootsave(name) foldersave(string)]
gettoken depvar indepvars : varlist


	qui{
	tempvar touse_wrap
	mark `touse_wrap' `if' `in'
	markout `touse_wrap' `varlist'
	}
if ("`noproj'"=="" & "`nogen'"=="nogen"){
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
	gettoken twoway_id aux: absorb
	gettoken twoway_t w: aux

	
	twowayreg `depvar' `indepvars', vce(`vce') `statadof'
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

** End of twowayreg.ado
