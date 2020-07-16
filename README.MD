# An Algorithm to Estimate the Two-Way Fixed Effect Model

This repo contains the ado file needed to run the Two-way fixed effect model in Stata. Run the ado file first or place it in the ado-files folder in Stata.

twowaywrapper is a command that runs the whole algorithm with many options for the user, the syntaxis is the following:
```
.varlist(numeric ts fv) [if] [in], [,ABSorb(varlist min=2 max=3) DROP GENerate(name)] [, NEWVars(name) REPLACE NOPROJ] [, VCE(name)]
```
In absorb option the user can define the two variables with the household id and the time id. From this option the user can choose two ways to deal with the redundant observations: drop them (typing drop after absorb) or generate a new variable that indicates which observations are not redundant (with gen() option).
In absorb option the user must also specify the weight variable. 

In the second list of options the user chooses how to perform the second step of the algorithm. The user has three options: create new variables that represent the residual projection of the variables in the list(in parenthesis the user choose the prefix of this new variables), replace variables with their projected counterpart or don't poject the variables if they were already projected.

With vce() option the user chooses how to calculate the standar errors:
```
.vce(robust): standard errors robust to heteroscedasticity but assumes no correlation within group or serial correlation.
.vce(vce):standard errors  proposed by Arellano (1987) robust to heteroscedasticity and serial correlation. 
.vce(vce_2):Arellano standard errors with a degree of freedom correction performed by Stata xtreg, fe
.If vce() is omitted: standard errors assuming homoscedasticity and no within  group correlation or serial correlation.
```

Here are some examples of the command:
```
.twowayregwrap y x*, absorb(hhid tid w) drop newv(w_) vce(robust)
.twowayregwrap w_y w_x1, absorb(hhid tid w) drop noproj vce(vce)
In the first line the user chooses to drop the redundant variables and creates new variables that are the projected variables, whereas in the second line the user takes advantage of the variables already projected to make a faster regression.
```

If the user wants to take more advantage of the algorithm it can be used separately.

twowayset defines the two variables with the household id and the time id. The user can choose to generate a new variable that refer to non-redundants variables or drop the redundant variables. For example:
```
.twowayset hhid tid, gen(sample)
.twowayset hhid tid, drop
.twowayset hhid tid w, gen(sample)
.twowayset hhid tid w, drop

```
The first two lines assign each observation the same weight. The second two lines assign to each observation the weight given by variable w.

projvar performs the second step of the algorithm. It is followed by a list of variables. These variables are projected onto the set of dummies, and the residual of this projection is saved as additional variables. Alternatively, the user may choose to replace the existing variables. For example:
```
.projvar y x1 x2, p(w_)
.projvar y x1 x2, replace
```

twowayreg performs the last step of the algorithm, making the regression of the projected variables. Here the user can calculate the standard errors as in twowaywrapper. For example:
```
.twowayreg w_dist w_mean_np w_x* if sample==1, vce
where the user has already choosen to generate a variable that references to non-redundants variables.
```
