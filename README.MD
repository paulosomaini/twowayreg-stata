# An Algorithm to Estimate the Two-Way Fixed Effect Model

This repo contains the ado file needed to run the Two-way fixed effect model in Stata. Run the ado file first or place it in the ado-files folder in Stata.

twowaywrapper is a command that runs the whole algorithm with many options for the user, the syntaxis is the following:
```
.varlist(numeric ts fv) [if] [in], [,ABSorb(varlist min=2 max=3) GENerate(namelist) NOGEN NOPROJ] [, NEWVars(name) REPLACE] [, VCE(namelist)]
```
In absorb option the user can define the two variables with the household id and the time id. 
In absorb option the user must also specify the weight variable. 
If the fixed effects variables are consecutive after dropping the redundants and missing values then the user can use the option nogen, but if they are not then the user will have to use the option generate() to create new fixed effects that are consecutives.
With NOPROJ the user can make the regression without creating any arrays or projecting new variables.

In the second list of options the user chooses how to perform the second step of the algorithm. The user has two options: create new variables that represent the residual projection of the variables in the list(in parenthesis the user choose the prefix of this new variables), or replace variables with their projected counterpart.

With vce() option the user chooses how to calculate the standar errors:
```
.robust: standard errors robust to heteroscedasticity but assumes no correlation within group or serial correlation.
.vce(cluster hhid):standard errors  proposed by Arellano (1987) robust to heteroscedasticity and serial correlation. 
.vce(cluster hhid) statadof:Arellano standard errors with a degree of freedom correction performed by Stata xtreg, fe
.If vce() is omitted: standard errors assuming homoscedasticity and no within  group correlation or serial correlation.
```

Here are some examples of the command:
```
.twowayregwrap y x*, absorb(hhid tid w) nogen newv(w_) robust
.twowayregwrap w_y w_x1, noproj vce(cluster hhid)
In the first line the user creates new variables that are the projected variables, whereas in the second line the user takes advantage of the variables already projected to make a faster regression.
```

If the user wants to take more advantage of the algorithm it can be used separately.

twowayset defines the two with the id and the time id and with them the program can create the required. With this command the user can also weight the observations. Depending on the database the user will have to generate new variables of groups or not. For example:
```
.twowayset hhid tid, gen(newids newts)
.twowayset hhid tid, nogen
.twowayset hhid tid w, gen(newids newts)
.twowayset hhid tid w, nogen

```
The first two lines assign each observation the same weight. The last two lines assign to each observation the weight given by variable w.
The user can use if and in in this command. This is extremely important if there is a variable that has missing values.

projvar performs the second step of the algorithm. It is followed by a list of variables. These variables are projected onto the set of dummies, and the residual of this projection is saved as additional variables. Alternatively, the user may choose to replace the existing variables. For example:

```
.projvar y x1 x2, p(w_)
.projvar y x1 x2, replace
```
If a missing value exist in one of the variables projected there will be an error indacting that the missing value should be deleted to project the matrix or the user will have to run twowayset again with an if that exclude missing values.

twowayreg performs the last step of the algorithm, making the regression of the projected variables. Here the user can calculate the standard errors as in twowaywrapper. For example:
```
.twowayreg w_y w_x*, vce(cluster hhid)

```
