# An Algorithm to Estimate the Two-Way Fixed Effect Model

The `twowayreg` package is an algorithm to efficiently estimate a two-way fixed effects model in STATA.  This algorithm is adapted from the [Matlab code](https://soma.people.stanford.edu/research) written for [this paper](https://soma.people.stanford.edu/sites/g/files/sbiybj4171/f/jem-2014-0008.pdf) by Paulo Somaini and Frank A Wolak.  

## Abstract
We present an algorithm to estimate the two-way fixed effect linear model. The algorithm relies on
the Frisch-Waugh-Lovell theorem and applies to ordinary least squares (OLS), two-stage least squares (TSLS)
and generalized method of moments (GMM) estimators. The coefficients of interest are computed using the
residuals from the projection of all variables on the two sets of fixed effects. Our algorithm has three desirable
features. First, it manages memory and computational resources efficiently which speeds up the computation of the estimates. Second, it allows the researcher to estimate multiple specifications using the same set
of fixed effects at a very low computational cost. Third, the asymptotic variance of the parameters of interest
can be consistently estimated using standard routines on the residualized data. 

## How do I use it?

There are two ways you can interface with the algorithm.


### Basic

The most basic way is using the command `twowaywrapper` to run ordinary least squares. 
This command runs the whole algorithm. The syntaxis is as follows:
```
.varlist(numeric ts fv) [if] [in] [using/], ABSorb(varlist min=2 max=3) [GENerate(namelist) NOPROJ] [, NEWVars(name) REPLACE] [, VCE(namelist)]
```

The varlist contains the dependent followed by the independent variables.
Adding `using` followed by a path will save a set of matrices in that path. Otherwise, the matrices will be saved in `eresults`. STATA may fail if it tries to store a matrix that exceeds `matsize` in ereturn. Adding `using` avoids that limitation.
The two levels of fixed effects are specified in the `absorb` option. The absorb option takes at least two variables with group identifiers. The optional third input can be a variable of analytic weights. 
If the group identifiers are not consecutive after dropping redundants and missing observations, then you have the option `generate()` to create new group identifiers.

In the second list of options you choose how to store the residualized variables. 
There are two options: 

1. Create new variables The residualized variables will be stored as new variables that will be named with a prefix specified in `newvars'.
2. Or replace variables with their residualized version.

With NOPROJ you can run the regression without creating any arrays or projecting new variables.

The vce() option allows you to choose how to compute standard errors:
```
.robust: standard errors robust to heteroscedasticity but assumes no correlation within group or serial correlation.
.vce(cluster hhid):standard errors  proposed by Arellano (1987) robust to heteroscedasticity and serial correlation. 
.vce(cluster hhid) statadof:Arellano standard errors with a degree of freedom correction performed by Stata xtreg, fe
.If vce() is omitted: standard errors assuming homoscedasticity and no within  group correlation or serial correlation.
```

Examples:
```
.twowayregwrap y x*, absorb(hhid tid w) newv(w_) robust
.twowayregwrap w_y w_x1, noproj vce(cluster hhid)
```
The first line creates new variables that are the residualized variables. The second line uses the variables already projected to save time.

### Advanced

Sometimes the same set of group identifiers have to be used to run several specifications.
In these cases, one can save computation time by reusing some of the previous computations.
The second and more advanced way to interact with the algorithm is to exploit this feature by using the functions  
`twowayset` `projvar` and `twowayreg`.

`twowayset` uses the group identifiers and weights to create a set of matrices that can be used for multiple specifications. You have the option to save these arrays adding `using` and a path. Otherwise, the matrices will be saved in `eresults`. STATA may fail if it tries to store a matrix that exceeds `matsize` in ereturn. Adding `using` avoids that limitation.
If the group identifiers are not consecutive after dropping redundants and missing observations, then you have the option `generate()` to create new group identifiers.

Examples:
```
.twowayset hhid tid, gen(newids newts)
.twowayset hhid tid
.twowayset hhid tid w, gen(newids newts)
.twowayset hhid tid w
.twowayset hhid tid using "../folder/x", gen(newids newts)
.twowayset hhid tid using "../folder/x"

```
The first two lines assign each observation a weight of one (default. The second pair of lines assign to each observation the weight given by variable w. The last two lines save the arrays in "folder" with the prefix "x"
The user can use if and in in this command. This is extremely important if there is a variable that has missing values.

`projvar` performs the second step of the algorithm. It is followed by a list of variables. These variables are projected onto the set of dummies, and the residual of this projection is saved as additional variables with a new prefix. Alternatively, you may choose to replace the existing variables. For example:

```
.projvar y x1 x2, p(w_)
.projvar y x1 x2, replace
.projvar y x1 x2 using "../folder/x", p(w_)
```
If you used `using` to save the matrices in `twowayset`, use the same path to recover them and use them in `projvar`. Otherwise, `projvar` will use the matrices stored in ereturn.
If a missing value exist in one of the variables projected there will be an error indacating that the variable cannot be residualized.

twowayreg performs the last step of the algorithm, making the regression of the projected variables. Here you can calculate the standard errors as in twowaywrapper. For example:
```
.twowayreg w_y w_x*, vce(cluster hhid)

```

## Whats in this repo?

* twowayreg.ado
* example.do


## Set up and Installation

### One-time use
Download this repository and extract it. From the directory where this README.MD file is located type ` do twowayreg.ado` in the R console.

### Manual Installation

TBD

### From SSC

TBD

## Appendix 

The variance formula will be:

![first eq](https://latex.codecogs.com/gif.latex?%5Cwidehat%7B%5Cmathcal%7BV%7D%7D%3Dq_%7Bc%7D%20%5Cwidehat%7B%5Cmathbf%7BV%7D%7D%5Cleft%28%5Csum_%7Bj%3D1%7D%5E%7BN%7D%20w_%7Bj%7D%20%5Cmathbf%7Bu%7D_%7Bj%7D%5E%7B%5Cprime%7D%20%5Cmathbf%7Bu%7D_%7Bj%7D%5Cright%29%20%5Cwidehat%7B%5Cmathbf%7BV%7D%7D)

where:

![second eq]https://latex.codecogs.com/gif.latex?%5Cwidehat%7B%5Cmathbf%7BV%7D%7D , which is typically a conventionally calculated variance matrix

![third eq]https://latex.codecogs.com/gif.latex?u_j , j=1,...,N a row vector of scores; and

qc, a constant degree of freedom correction necessary to calculate correctly the variance-covarinace matrix.

wj are the user-specified weights, equal to 1 if no weights are specified.

The degree of freedom correction:

![fourth eq]https://latex.codecogs.com/gif.latex?q_%7B%5Cmathrm%7Bc%7D%7D%3D%5Cfrac%7BN-1%7D%7BN-k-d-t%7D%20%5Cfrac%7BM%7D%7BM-1%7D

Where k is the number of controls without dummies, t is the number of time fixed effects without the redundants, d is the number of group fixed effects without the redundants, N is the number of observations(each observation is equal to 1 if no weights are specified), and M is the number of cluster if vce(cluster clustvar) is specified.   

If vce(cluster clusvar) is not specified then N=M and 
![fith eq]https://latex.codecogs.com/gif.latex?q_%7B%5Cmathrm%7Bc%7D%7D%3DN%20/%28N-k-d-t%29

## References

Arellano, M. (1987), Computing Robust Standard Errors for Within-Groups Estimators, Oxford Bulletin of Economics and
Statistics, 49, issue 4, p. 431–434. 

Somaini, P. and F.A. Wolak, (2016), An Algorithm to Estimate the Two-Way Fixed Effects Model, Journal of Econometric Methods, 5, issue 1, p. 143-152.
