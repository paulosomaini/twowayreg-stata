capture program drop twowayset
capture program drop twowaysample
capture mata mata drop sparse()
capture mata mata drop dofadj()
capture mata mata drop proddiag()
capture mata mata drop diagprod()
capture mata mata drop diagminus()
capture mata mata drop projDummies()
capture mata mata drop projDummies_save()
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
 
 
void projDummies_save()
{
real matrix D, DH1, DH, CinvHHDH, AinvDDDH, A, B, C
real colvector DD, HH, invDD, invHH
real scalar N, T
string scalar twoWaynewid,twoWaynewt, w,sampleVarName, root
D=.
//printf("Hola Paulo, todo functiona hasta aqui.")

w = st_local("twoway_w")
root =st_local("root")
sampleVarName = st_local("twoway_sample")
if (w==""){
D = st_data(.,("twoWaynewid","twoWaynewt"),sampleVarName)
D = (D,J(rows(D),1,1))
}
else {
D = st_data(.,("twoWaynewid","twoWaynewt",w),sampleVarName)
}
//printf(sampleVarName)
//printf("Incluso aca\n")
//D[1..10,]
//printf("y aca")

DH1=sparse(D)
//printf("Wohoo")
DD=quadrowsum(DH1)
HH=quadcolsum(DH1)'
HH=HH[1..cols(DH1)-1]



DH=DH1[.,1..cols(DH1)-1]

 
invDD=DD:^-1 
invHH=HH:^-1

N=colmax(D)[.,1]
T=colmax(D)[.,2]
saveMat(root,"twoWayN1", N)
saveMat(root,"twoWayN2", T)
saveMat(root,"twoWayinvDD", invDD)
saveMat(root,"twoWayinvHH", invHH)
//st_matrix("twoWayD", D...)
 if (N<T)
		{
        
        CinvHHDH=diagprod(invHH,DH')
		A=qrinv(diagminus(DD,CinvHHDH'*DH'))
		//st_matrix("CinvHHDH",CinvHHDH)
        B=-A*CinvHHDH'
		saveMat(root,"twoWayCinvHHDH", CinvHHDH)
		saveMat(root,"twoWayA", A)
		saveMat(root,"twoWayB", B)
		
		
		}
    else
	{
        AinvDDDH=diagprod(invDD,DH)
		C=qrinv(diagminus(HH,AinvDDDH'*DH))
		//st_matrix("AinvDDDH",AinvDDDH)
        B=-AinvDDDH*C
		saveMat(root,"twoWayAinvDDDH", AinvDDDH)
		saveMat(root,"twoWayC", C)
		saveMat(root,"twoWayB", B)
		
    }
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
syntax varlist(min=2 max=3) [if] [in], [, SAVE]
gettoken twoway_id aux: varlist
gettoken twoway_t twoway_w: aux
//summ `varlist'
// I need to make it robust to non 1,2,3... ids.

if ("`save'"=="save") {
	if ("`root'" == "") {
	local root="last"
	}
	egen twoWaynewid= group(`twoway_id')
	egen twoWaynewt= group(`twoway_t')
//di in gr "`twoway_id'"
//di in gr "`twoway_t'"

	tempvar twoway_sample
	mark `twoway_sample' `if' `in'
	markout `twoway_sample' `varlist'
	mata projDummies_save()
//di in gr "Checkpoint 1"
//ret li
//di in gr "Checkpoint 2"
	scalar twoWayid="`twoway_id'"
	scalar twoWayt="`twoway_t'"
	scalar twoWayw="`twoway_w'"
	scalar twoWayif="`if'"
	scalar twoWayin="`in'"
}

else{

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
}

//return post r(B), esample(`twoway_sample') 
//obs(`nobs') dof(`dof')



end


program define twowaysample, sortpreserve
version 11
syntax varlist(min=2 max=3) [if] [in], Generate(name) [Replace] [, SAVE]

if ("`save'"=="save"){
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


}

else{
gettoken twoWaynewid aux: varlist
gettoken twoWaynewt w: aux


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
capture mata mata drop projVar_save()
capture mata mata drop projVar()

mata
void projVar_save()
{
	real matrix V, varIn, D,aux,delta,tau,varOut,A,B,CinvHHDH,AinvDDDH,C
	real colvector invHH,invDD,Dy,Ty
	real scalar N,T
	string scalar id, t, currvar,newvar,sampleVarName,w, root
	currvar = st_local("currvar")
	newvar = st_local("newvar")
	id=st_strscalar("twoWayid")
	root =st_local("root")
	N=readMat(root,"twoWayN1")
	T=readMat(root,"twoWayN2")
	//D=readMat(root,"twoWayD")
	w=st_strscalar("twoWayw")
	t=st_strscalar("twoWayt")
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
	B=readMat(root,"twoWayB")
	
	//rows(Ty)
    //cols(Ty)
	//rows(Dy)
	//cols(Dy)
			

	 if (N<T)
			{
			
			A=readMat(root,"twoWayA")
			invHH=readMat(root,"twoWayinvHH")
			CinvHHDH=readMat(root,"twoWayCinvHHDH")
			//printf("b")
			delta=A*Dy+B*Ty
			tau=B'*(Dy-CinvHHDH'*Ty)+(invHH:*Ty) \0
			}
		else
		{
			//printf("1")
			C=readMat(root,"twoWayC")
			invDD=readMat(root,"twoWayinvDD")
			AinvDDDH=readMat(root,"twoWayAinvDDDH")
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
	//D=readMat(root,"twoWayD")
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
syntax varlist, [Prefix(name)] [REPLACE] [, SAVE]

if ("`save'"=="save"){
	tempvar twoway_sample
loc tif=twoWayif
loc tin=twoWayin
mark `twoway_sample' `tif' `tin'
markout `twoway_sample' `varlist'
if ("`prefix'" == "") {
	local prefix="proj_"
	}
if ("`root'" == "") {
	local root="last"
	}

foreach currvar of varlist `varlist' {
	local newvar="`prefix'`currvar'"
	if ("`replace'" != "") {
	local newvar="`currvar'"
	}
	else {
	qui gen `newvar'=.
	}
	//di "`currvar'"
	//di "`newvar'"
	mata projVar_save()
}
drop twoWaynewid
drop twoWaynewt
}

else{
	tempvar twoway_sample
loc tif=twoWayif
loc tin=twoWayin
mark `twoway_sample' `tif' `tin'
markout `twoway_sample' `varlist'
//mata mata describe
//summ `varlist'
//summ `twoway_sample'
// I need to make it robust to non 1,2,3... ids.


foreach currvar of varlist `varlist' {
	local newvar="`prefix'`currvar'"
	if ("`replace'" != "") {
	local newvar="`currvar'"
	}
	else {
	qui gen `newvar'=.
	}
	//di "`currvar'"
	//di "`newvar'"
	mata projVar()
}
drop twoWaynewid
drop twoWaynewt

}
end