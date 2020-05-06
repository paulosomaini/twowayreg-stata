# An Algorithm to Estimate the Two-Way Fixed Effect Model

This repo contains the ado file needed to run the Two-way fixed effect model in Stata. Run the ado file first, or place it in the ado-files folder in Stata.

twowayset defines the two variables with the household id and the time id. For example:
```
.twowayset hhid tid
.twowayset hhid tid w
```
The first line uses assigns each observation the same weight. The second line assigns to each observation the weight given by variable w.
The routine twowayset computes the first step of the algorithm and stores the results in the working directory. To store the results in a different folder or to avoid overwriting previous results use the option root. To store the results in files with a prefix abs in the folder C:\Results, write:
```
.twowayset hhid tid w, root(“C:\Results\abc”)
```
projvar performs the second step of the algorithm. It is followed by a list of variables. These variables are projected onto the set of dummies, and the residual of this projection is saved as additional variables. Alternatively, the user may choose to replace the existing variables. For example:
```
.projvar y x1 x2, p(w_)
```
This command generates variables w_y, w_x1 and w_x2. Each variable contains the residual of the projection of the original variable on the set of dummies. 
```
.projvar y x1 x2, replace
```
This command replaces y, x1 and x2 with the residual of the projection of the original variable on the set of dummies. To use a different set of first-step results or to use the results stored in a different folder use the option root, e.g:
```
.projvar y x1 x2, replace root(“C:\Results\abc”)
```
Once variables have been projected onto the set of dummies and their residuals are stored in the memory, the estimation can be performed using the standard STATA commands:
```
.reg w_y w_x1, noc
.reg w_y w_x1 w_x2, noc 
.reg w_y w_x1 w_x2, noc robust
```
Because these variables are residuals, we do not need to estimate a intercept in these functions. The resulting standard errors can be used. The resulting standard errors will be different than those calculated using the MATLAB version of this algorithm (or using the full set of dummies) because of small sample degrees of freedom corrections. 
Minimum working example: see Example.do.
