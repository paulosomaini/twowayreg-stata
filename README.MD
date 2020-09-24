# An Algorithm to Estimate the Two-Way Fixed Effect Model

The `twfem` package is an algorithm to efficiently estimate a two-way fixed effects model in STATA.  This algorithm is adapted from the [Matlab code](https://soma.people.stanford.edu/research) written for [this paper](https://soma.people.stanford.edu/sites/g/files/sbiybj4171/f/jem-2014-0008.pdf) by Paulo Somaini and Frank A Wolak.  

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

The most basic way is using the command `twfem` to run ordinary least squares. 
This command runs the whole algorithm. The syntaxis is as follows:
```
.command varlist(numeric ts fv) [if] [in] [using/], ABSorb(varlist min=2 max=3) [GENerate(namelist) NOPROJ] [, NEWVars(name) REPLACE] [, VCE(namelist)]
```

Command is an estimation command. For example you could choose to run: reg, ivreg, sureg, etc. 
The varlist contains the dependent followed by the independent variables. The syntaxis of the varlist replicates the syntaxis of the original command.
The code creates and store some matrices that can be used later. Adding `using` followed by a path will save a set of matrices in that path. Otherwise, the matrices will be saved in `eresults`. STATA may throw an error if it tries to store a matrix that exceeds `matsize` in ereturn. Adding `using` avoids that limitation.
The two levels of fixed effects are specified in the `absorb` option. The absorb option takes at least two variables with group identifiers. The optional third input can be a variable of analytic weights. 
If the group identifiers are not consecutive after dropping redundants and missing observations, then you have the option `generate()` to create new group identifiers.

In the second list of options you choose how to store the residualized variables. 
There are two options: 

1. Create new variables: the residualized variables will be stored as new variables that will be named with a prefix specified in `newvars'.
2. Replace variables with their residualized version.

If you have already residualized all the variables that will be included in the specification, you may include them directly and add the option NOPROJ.
This option assumes that all the included variables have been residualized. Residualizing the variables is the most time-consuming step. This option is useful
to reuse variables that have been already residualized.

The vce() option allows you to choose how to compute standard errors:
```
.vce(robust): standard errors robust to heteroscedasticity but assumes no correlation within group or serial correlation.
.vce(cluster hhid):standard errors  proposed by Arellano (1987) robust to heteroscedasticity and serial correlation. 
.If vce() is omitted: standard errors assuming homoscedasticity and no within  group correlation or serial correlation.
```

Examples:
```
.twfem reg y x1 x2 x3 x4, absorb(hhid tid w) newv(w_) vce(robust)
.twfem reg w_y w_x1, noproj vce(cluster hhid)
.twfem ivregress 2sls w_y w_x1 (w_x2= w_x3), noproj vce(robust)
```
The first line stores the residualized variables as new variables using the prefix `w_`. The second line reuses these variables to save time. The third line uses them again with a different command. 

### Advanced

Sometimes the same set of group identifiers have to be used to run several specifications.
In these cases, one can save computation time by reusing some of the previous computations.
The second and more advanced way to interact with the algorithm is to exploit this feature by using the functions  
`twset` `twres` and `twest`.

`twset` uses the group identifiers and weights to create a set of matrices that are used to residualize variables. You have the option to save these arrays adding `using` and a path. Otherwise, the matrices will be saved in `eresults`. STATA may throw an error if it tries to store a matrix that exceeds `matsize` in ereturn. Adding `using` avoids that limitation.
If the group identifiers are not consecutive integers (after dropping redundants and missing observations), then you have the option `generate()` to create new group identifiers.
The code may throw an error if the group identifiers are not consecutive integers and the option `generate()` is missing.

Examples:
```
.twset hhid tid, gen(newids newts)
.twset hhid tid
.twset hhid tid w, gen(newids newts)
.twset hhid tid w
.twset hhid tid using "../folder/x", gen(newids newts)
.twset hhid tid using "../folder/x"

```
The first two lines assign each observation a weight of one (default). The second pair of lines assign to each observation the weight given by variable w. The last two lines save the arrays in "folder" with the prefix "x"
You can use if and in in this command. This is extremely important if there is a variable that has missing values.

`twres` performs the second step of the algorithm. It is followed by a list of variables. These variables are projected onto the set of dummies, and the residual of this projection is saved as additional variables with a new prefix. Alternatively, you may choose to replace the existing variables. For example:

```
.twres y x1 x2, p(w_)
.twres y x1 x2, replace
.twres y x1 x2 using "../folder/x", p(w_)
```
If you used `using` to save the matrices in `twset`, use the same path to recover them and use them in `twres`. Otherwise, `twres` will use the matrices stored in ereturn.
If a missing value exist in one of the variables projected there will be an error indacating that the variable cannot be residualized.

twest performs the last step of the algorithm, making the regression of the projected variables. You can select the estimation command and the formula for the standard errors as in twfem. For example:
```
.twest reg w_y w_x*, vce(cluster hhid)

```

If you already run `twset` you can save the matrices in a folder to not re-calculate them. 
The comand `twsave` will complete this task and the command `twload` load the matrices saved previously. Both of them can be followed by `using` to select the path where the matrices will be saved. If using option is ommited the arrays will we stored in the current directory.
```
.twsave using "../folder/x"
....
.twload using "../folder/x"
```
In this example the matrices will be saved in "folder" with the prefix "x". Now you may close STATA, open it again, and use `twload` to load the arrays.

## What's in this repo?

* twset.ado
* twset.sthlp
* twres.ado
* twres.sthlp
* twest.ado
* twest.sthlp
* twsave.ado
* twsave.sthlp
* twload.ado
* twload.sthlp
* twfem.ado
* twfem.sthlp
* example.do


## Set up and Installation

### One-time use
Download this repository and extract it. From the directory where this README.MD file is located type ` do twest.ado` in the STATA console.

### Manual Installation

To install the commands without ssc install; download the commands and put them in PERSONAL or PLUS folder (to search where are this folders run "sysdir" in STATA).	

### From SSC

TBD

## Appendix 

The variance formula will be:

  ![first eq](https://latex.codecogs.com/gif.latex?%5Cwidehat%7B%5Cmathcal%7BV%7D%7D%3Dq_%7Bc%7D%20%5Cwidehat%7B%5Cmathbf%7BV%7D%7D)

where:

![second eq](https://latex.codecogs.com/gif.latex?%5Cwidehat%7B%5Cmathbf%7BV%7D%7D) , which is variance-covariance matrix calculated by the command choosen in twest or twfem.

`qc`, a degree of freedom correction necessary to calculate correctly the variance-covarinace matrix.

The degree of freedom correction:

   ![third eq](https://latex.codecogs.com/gif.latex?q_%7B%5Cmathrm%7Bc%7D%7D%3DN-k/%28N-k-d-t%29)

Where k is the number of controls without dummies, t is the number of time fixed effects without the redundants, d is the number of group fixed effects without the redundants, N is the number of observations(each observation is equal to 1 if no weights are specified).

Based on [this paper](http://cameron.econ.ucdavis.edu/research/Cameron_Miller_JHR_2015_February.pdf) by A. Colin Cameron, and Douglas L. Miller if a fixed effect variable is also a cluster variable (or is nested within a cluster variable), we do not consider it when computing the absorbed degrees of freedom (the degree of freedom correction will be smaller). This is done to replicate the within estimation which leads to the correct degree of freedom correction. 

## References

Arellano, M. (1987), Computing Robust Standard Errors for Within-Groups Estimators, Oxford Bulletin of Economics and
Statistics, 49, issue 4, p. 431–434. 

Cameron, A. C., & Miller, D. L. (2015). A practitioner’s guide to cluster-robust inference. Journal of human resources, 50(2), 317-372.

Somaini, P. and F.A. Wolak, (2016), An Algorithm to Estimate the Two-Way Fixed Effects Model, Journal of Econometric Methods, 5, issue 1, p. 143-152.
