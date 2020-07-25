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
 
  real matrix readMat(string s,string n)
 {
  real matrix X
  real scalar fh
fh = fopen(s+"_"+n, "r")
X = fgetmatrix(fh)
fclose(fh)
return(X)
 }
 
  
 void saveMat(string s,string n,real matrix X)
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
D = st_data(.,("twoWaynewid","twoWaynewt",w),sampleVarName)
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

st_numscalar("e(H)",N)
st_numscalar("e(T)",T)
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
syntax varlist(min=2 max=3) [if] [in], [GENerate(namelist) Nogen]
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

		tempvar var1
		gen `var1'= `twoway_id'
		tempvar var2
		gen `var2'= `twoway_t'
	}
	
	else{
		gettoken var1 var2: generate
		egen `var1'= group(`twoway_id')
		egen `var2'= group(`twoway_t')
		}
		
	/*
	[Gen(names)]- for new names-- if gen!="" {group(name`1',name`2')} else{check consecutive}
	check consecutive numbers.
	*/
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
	real colvector invHH,invDD,Dy,Ty, linear_index
	real scalar N,T
	string scalar newid, newt, currvar,newvar,sampleVarName,w
	currvar = st_local("currvar")
	newvar = st_local("newvar")
	newid=st_local("var1")
	N=st_numscalar("e(H)")
	T=st_numscalar("e(T)")
	w=st_strscalar("twoWayw")
	newt=st_local("var2")
	sampleVarName = st_local("touse_proj")
	V = st_data(.,(newid, newt,currvar),sampleVarName)
	varIn=V[.,3]
	
	if (w==""){
	D = st_data(.,(newid, newt),sampleVarName)
	D = (D,J(rows(D),1,1))
	}
	else {
	D = st_data(.,(newid, newt,w),sampleVarName)
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
	linear_index = st_data(.,("linear_index"),sampleVarName)
	st_store(linear_index, newvar, varOut)
	//printf("5")
}
end



program define projvar, nclass
version 11
syntax varlist, [ABSorb(varlist)] [Prefix(name)] [REPLACE] 
	
	gettoken depvar indepvars : varlist
    _fv_check_depvar `depvar'
    fvexpand `indepvars' 
	
	
	qui{
	tempvar touse_proj
	gen byte `touse_proj'= e(sample)
	tempvar  touse_check 
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
	
	
	
	gen linear_index = _n	
	gettoken var1 var2: absorb
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

drop linear_index
scalar N= e(H)
scalar T= e(T)
matrix invDD=e(invDD)
matrix invHH=e(invHH)
if (N<T){
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

mata
void projDummies_save()
{
real matrix D, CinvHHDH, AinvDDDH, A, B, C
real colvector  invDD, invHH
real scalar N, T
string scalar twoWaynewid,twoWaynewt, w,sampleVarName, root

root =st_local("root")
N=st_numscalar("N")
T=st_numscalar("T")
invDD=st_matrix("invDD")
invHH=st_matrix("invHH") 
saveMat(root,"twoWayN1", N)
saveMat(root,"twoWayN2", T)
saveMat(root,"twoWayInvDD", invDD)
saveMat(root,"twoWayInvHH", invHH)

//st_matrix("twoWayD", D...)
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


program define twowaysave, rclass
version 11
syntax  [if] [in] [, FOLDER(string)] [Root(name)]

if ("`folder(`string')'"=="`folder(`string')'"){
	cd "`folder'"
}

	tempvar twoway_sample
	mark `twoway_sample' `if' `in'
	markout `twoway_sample' `varlist'
	mata projDummies_save()
	scalar twoWayid="`twoway_id'"
	scalar twoWayt="`twoway_t'"
	scalar twoWayw="`twoway_w'"
	scalar twoWayif="`if'"
	scalar twoWayin="`in'"


	
drop twoWaynewid
drop twoWaynewt

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

root =st_local("root")
w = st_local("twoway_w")
sampleVarName = st_local("twoway_sample")
N=readMat(root,"twoWayN1")
T=readMat(root,"twoWayN2")
invDD=readMat(root,"twoWayInvDD")
invHH=readMat(root,"twoWayInvHH")
st_numscalar("e(H)",N)
st_numscalar("e(T)",T)
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


program define twowayload, rclass
version 11
syntax varlist(min=2 max=3) [if] [in],[, FOLDER(string)] [,DROP] [,Generate(name)] [Replace] [Root(name)] 
gettoken twoway_id aux: varlist
gettoken twoway_t twoway_w: aux

if ("`folder(`string')'"=="`folder(`string')'"){
	cd "`folder'"
}
//summ `varlist'
// I need to make it robust to non 1,2,3... ids.


	egen twoWaynewid= group(`twoway_id')
	egen twoWaynewt= group(`twoway_t')
	local newvars twoWaynewid twoWaynewt

	tempvar twoway_sample
	mark `twoway_sample' `if' `in'
	markout `twoway_sample' `newvars'



	mata projDummies_load()
//di in gr "Checkpoint 1"
//ret li
//di in gr "Checkpoint 2"
	scalar twoWayw="`twoway_w'"
	scalar twoWayif="`if'"
	scalar twoWayin="`in'"
	
  gettoken twoWaynewid aux: varlist
  gettoken twoWaynewt w: aux

  if("`drop'"=="drop"){
		if !("`w'"==""){
		replace `w' = . if `w'<=0
	}
	
	qui{
	gen aux=1
	
	tempvar howmany
	count if aux == 1

	while `r(N)' {
	bys `twoWaynewid': gen `howmany' = _N if aux
	replace aux = 0 if `howmany' == 1
	drop `howmany'

	bys `twoWaynewt': gen `howmany' = _N if aux
	replace aux = 0 if `howmany' == 1
		
	count if `howmany' == 1
	drop `howmany'
	drop if aux!=1
	drop aux
	}
	}

	
}	
	
else if ("`drop'"!="drop") {
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
		bys `twoWaynewid': gen `howmany' = _N if `generate'
		replace `generate' = 0 if `howmany' == 1
		drop `howmany'

		bys `twoWaynewt': gen `howmany' = _N if `generate'
		replace `generate' = 0 if `howmany' == 1
		
		count if `howmany' == 1
		drop `howmany'
	}
	}
}





end

capture program drop twowayreg 

program define twowayreg, eclass sortpreserve
    version 14
 
    syntax varlist(numeric ts fv) [if] [in],[,ROBUST VCE(namelist) statadof] 
    gettoken depvar indepvars : varlist
    _fv_check_depvar `depvar'
    fvexpand `indepvars'
	
	gettoken twoway_id aux: absorb
	gettoken twoway_t w: aux
	qui{
	tempvar touse_reg
	mark `touse_reg' `if' `in'
	replace `touse_reg'= `touse_reg' & e(sample)
	}
	scalar N= e(H)
	scalar T= e(T)
	matrix invDD=e(invDD)
	matrix invHH=e(invHH)
	if (N<T){
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
   	qui{
  	regress `depvar' `indepvars' if `touse_reg' , noc robust
	scalar vadj = e(df_r)/(e(df_r)- N - T)
	scalar df_r1= e(df_r) - N - T
    }
 }
    else if ("`vce(namelist)'"== "`vce(namelist)'" & "`statadof'"== ""){

   qui{
  	regress `depvar' `indepvars' if `touse_reg' , noc vce(`vce')
	qui{
		scalar df_r= e(N)-e(df_m)-1
	}
	scalar df_r1= e(df_r)
	scalar vadj = df_r/(df_r- N - T)
	    }
 }

  else if ("`vce(namelist)'"== "`vce(namelist)'" & "`statadof'"== "statadof"){

   qui{
  	regress `depvar' `indepvars' if `touse_reg' , noc vce(`vce')
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
  qui{
  	regress `depvar' `indepvars' if `touse_reg' , noc
	scalar df_r1= e(df_r)-N-T
	scalar vadj = e(df_r)/(e(df_r)- N - T)
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
  ereturn scalar N_1= N_1
  ereturn scalar R2= R2
  ereturn scalar F= F
  ereturn scalar df_m=df_m
  ereturn scalar df_r1= df_r1
  ereturn scalar rtms= rtms
  ereturn scalar H= N
  ereturn scalar T= T
  ereturn matrix invDD= invDD
  ereturn matrix invHH= invHH
  if (N<T){
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
version 14 
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
version 14 
syntax varlist(numeric ts fv) [if] [in], [,ABSorb(varlist min=2 max=3) NOPROJ] [, NEWVars(name) REPLACE] [, VCE(name)] [, SAVE rootsave(name) foldersave(string)]
gettoken depvar indepvars : varlist


	qui{
	tempvar touse_wrap
	mark `touse_wrap' `if' `in'
	markout `touse_wrap' `varlist'
	}

if ("`noproj'"==""){
	twowayset `absorb' if `touse_wrap'
	
	 if ("`NEWVars(`name')'"=="`newvars(`name')'" & "`replace'"==""){
		capture confirm variable `newvars'
		if !_rc {
                 di "{err} There is at least one variable with the same prefix chosen, please change the prefix or drop the variable"
				}
        else {
              projvar `depvar' `indepvars', p(`newvars')
			  twowayreg `newvars'* if `touse_wrap', `vce'
			  }
				
	}
		
		
	else if ("`NEWVars(`name')'"=="" & "`replace'"=="replace" ){
		projvar `depvar' `indepvars' , replace
		twowayreg `depvar' `indepvars' if `touse_wrap', `vce'
		
	}
}	

else if ("`noproj'"=="noproj"){
	gettoken twoway_id aux: absorb
	gettoken twoway_t w: aux
	egen twoWaynewid= group(`twoway_id')
	egen twoWaynewt= group(`twoway_t')
	
	twowayreg `depvar' `indepvars' if `touse_wrap', `vce'
}
drop twoWaynewid twoWaynewt

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