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

newid=st_matrix("twoWaynewid")
newt=st_matrix("twoWaynewt")
w = st_local("twoway_w")
sampleVarName = st_local("twoway_sample")
if (w==""){
D = st_data(.,("twoWaynewid","twoWaynewt"),sampleVarName)
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

 
program define twowayset, rclass
version 11
syntax varlist(min=2 max=3) [if] [in], [,DROP] [,Generate(name)] [Replace] 
gettoken twoway_id aux: varlist
gettoken twoway_t twoway_w: aux
//summ `varlist'
// I need to make it robust to non 1,2,3... ids.


	egen twoWaynewid= group(`twoway_id')
	egen twoWaynewt= group(`twoway_t')
	local newvars twoWaynewid twoWaynewt

	tempvar twoway_sample
	mark `twoway_sample' `if' `in'
	markout `twoway_sample' `newvars'



	mata projDummies()
//di in gr "Checkpoint 1"
//ret li
//di in gr "Checkpoint 2"
	scalar twoWayw="`twoway_w'"
	scalar twoWayif="`if'"
	scalar twoWayin="`in'"

//return post r(B), esample(`twoway_sample') 
//obs(`nobs') dof(`dof')



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


capture program drop projvar
capture mata mata drop projVar()


mata
void projVar()
{
	real matrix V, varIn, D,aux,delta,tau,varOut,A,B,CinvHHDH,AinvDDDH,C
	real colvector invHH,invDD,Dy,Ty
	real scalar N,T
	string scalar newid, newt, currvar,newvar,sampleVarName,w
	currvar = st_local("currvar")
	newvar = st_local("newvar")
	newid=st_local("newid")
	N=st_numscalar("e(H)")
	T=st_numscalar("e(T)")
	w=st_strscalar("twoWayw")
	newt=st_local("newt")
	sampleVarName = st_local("twoway_sample")
	V = st_data(.,("twoWaynewid","twoWaynewt",currvar),sampleVarName)
	varIn=V[.,3]
	
	if (w==""){
	D = st_data(.,("twoWaynewid","twoWaynewt"),sampleVarName)
	D = (D,J(rows(D),1,1))
	}
	else {
	D = st_data(.,("twoWaynewid","twoWaynewt",w),sampleVarName)
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
	st_store(., newvar, varOut)
	//printf("5")
}
end



program define projvar, nclass
version 11
syntax varlist, [Prefix(name)] [REPLACE] 

	tempvar twoway_sample
	loc tif=twoWayif
	loc tin=twoWayin
	mark `twoway_sample' `tif' `tin'
	markout `twoway_sample' `varlist'


	foreach currvar of varlist `varlist' {
		local newvar="`prefix'`currvar'"
		if ("`replace'" != "") {
		local newvar="`currvar'"
		}
		else {
		qui gen `newvar'=.
		}
	mata projVar()
	
	}

qui{
	ereturn list
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
syntax varlist(min=2 max=3) [if] [in], [,DROP] [,Generate(name)] [Replace] 
gettoken twoway_id aux: varlist
gettoken twoway_t twoway_w: aux
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
 
    syntax varlist(numeric ts fv) [if] [in], [,ROBUST VCE VCE_2] 
    gettoken depvar indepvars : varlist
    _fv_check_depvar `depvar'
    fvexpand `indepvars' 
	marksample touse
 
   if ("`robust'"=="robust"){
   	qui{
  	regress `depvar' `indepvars' if `touse', noc robust
	scalar vadj = e(df_r)/(e(df_r)- N - T)
	scalar df_r1= e(df_r) - N - T
    }
 }
 
  else if ("`vce_2'"== "vce_2"){

   qui{
  	regress `depvar' `indepvars' if `touse', noc vce(cluster twoWaynewt)
	scalar N_1=e(N)
	scalar df_m= e(df_m)
	qui{
		scalar df_r= e(N)-e(df_m)-1
	}
	scalar df_r1= e(df_r)
	scalar vadj = (N_1-1)*(df_r/(df_r - 1))/(N_1 - df_m - 1)
	    }

 }

   else if ("`vce'"== "vce"){

   qui{
  	regress `depvar' `indepvars' if `touse', noc vce(cluster twoWaynewt)
	qui{
		scalar df_r= e(N)-e(df_m)-1
	}
	scalar df_r1= e(df_r)
	scalar vadj = df_r/(df_r- N - T)
	    }
 }

 
  else{
  qui{
  	regress `depvar' `indepvars' if `touse', noc
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

  eret post b V
  ereturn scalar N_1= N_1
  ereturn scalar R2= R2
  ereturn scalar F= F
  ereturn scalar df_m=df_m
  ereturn scalar df_r1= df_r1
  ereturn scalar rtms= rtms
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

