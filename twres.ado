capture program drop twres 
capture mata mata drop sparse()
capture mata mata drop excludemissing()
capture mata mata drop proddiag()
capture mata mata drop diagprod()
capture mata mata drop diagminus()
capture mata mata drop projVar()
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
 

void projVar()
{
	real matrix V, varIn, D,aux,delta,tau,varOut,A,B,CinvHHDH,AinvDDDH,C
	real colvector invHH,invDD,Dy,Ty
	real scalar N,T, correction_rank, save_to_e
	string scalar newid, newt, currvar,newvar,sampleVarName,w,linear_index,var1,var2, root
	root=st_local("using")
	save_to_e=st_numscalar("save_to_e")
	currvar = st_local("currvar")
	newvar = st_local("newvar")
	var1 = st_local("var1")
	var2 = st_local("var2")
	
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
	Dy=rowsum(aux)
	Dy=Dy
	Ty=colsum(aux)
	Ty=Ty[1,1..cols(aux)-1]'
	N=st_numscalar("e(dimN)")
	T=st_numscalar("e(dimT)")
	correction_rank=st_numscalar("e(rank_adj)")
	//load the matrices from eresults
	if (save_to_e>0){
		
		B=st_matrix("e(B)")
	}
	else{
		//load the matrices from using option
		N=readMat(root,"twoWayN1")
		T=readMat(root,"twoWayN2")
		B=readMat(root,"twoWayB")
	}

	 if (N<T)
			{
			if (save_to_e>0){
				A=st_matrix("e(A)")
				invHH=st_matrix("e(invHH)")
				CinvHHDH=st_matrix("e(CinvHHDH)")
			}
			else{
				A=readMat(root,"twoWayA")
				invHH=readMat(root,"twoWayinvHH")
				CinvHHDH=readMat(root,"twoWayCinvHHDH")
			}
			delta=A*Dy+B*Ty
			tau=B'*(Dy-CinvHHDH'*Ty)+(invHH:*Ty) \0
			}
		else
		{
			if (save_to_e>0){
				C=st_matrix("e(C)")
				invDD=st_matrix("e(invDD)")
				AinvDDDH=st_matrix("e(AinvDDDH)")
			}
			else{
				C=readMat(root,"twoWayC")
				invDD=readMat(root,"twoWayinvDD")
				AinvDDDH=readMat(root,"twoWayAinvDDDH")
			}
			delta=(invDD:*Dy)+B*(Ty-AinvDDDH'*Dy)
			tau=B'*Dy+C*Ty \0 
		}

	varOut=(varIn-delta[V[.,1]]-tau[V[.,2]]):*sqrt(D[.,3])
	st_store(st_data(.,linear_index,sampleVarName), newvar, varOut)

}

 end
 


program define twres, nclass
version 11
syntax varlist [using/], [Prefix(name)] [REPLACE]


foreach currvar of varlist `varlist'{ 
	if ("`replace'"=="") {
		capture confirm variable `prefix'`currvar'
		local rc = !_rc
		if !_rc {
				di "{err} The variable `prefix'`currvar' already exists."
				exit !_rc
			}
		}
}
	
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
	
	*check if e() has not been rewritten
	capt confirm scalar e(dimN)
	if _rc { 
	   di "{err} The e() has been rewritten, please run twset again."
	   exit
	}
	
	*check that there is using in the command or not.
	capt assert inlist( "`using/'", "")
	if !_rc {    
		scalar save_to_e=1
			*save in scalars and arrays the macros.
			scalar dimN= e(dimN)
			scalar dimT= e(dimT)
			scalar rank_adj=e(rank_adj)
			matrix invDD=e(invDD)
			matrix invHH=e(invHH)
			if (dimN<dimT){
				matrix CinvHHDH=e(CinvHHDH)
				matrix A= e(A)
				matrix B=e(B)
				}
			else{
				matrix AinvDDDH=e(AinvDDDH)
				matrix C= e(C)
				matrix B=e(B)
				}
		
		}
	else{
		scalar save_to_e=0
		}
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


end
