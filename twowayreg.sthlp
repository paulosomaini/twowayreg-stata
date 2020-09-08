{smcl}
{* *! version 1.0  6 Sep 2020}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "Install command2" "ssc install command2"}{...}
{vieweralsosee "Help command2 (if installed)" "help command2"}{...}
{viewerjumpto "Syntax" "twest##syntax"}{...}
{viewerjumpto "Description" "twest##description"}{...}
{viewerjumpto "Options" "twest##options"}{...}
{viewerjumpto "Remarks" "twest##remarks"}{...}
{viewerjumpto "Examples" "twest##examples"}{...}
{title:Title}
{phang}
{bf:twfe} {hline 2} Algorithm to efficiently estimate a two-way fixed effects model based on Somaini and Wolak(2016).

{marker syntax}{...}
{title:Syntax}
{p 8 17 2}
{cmdab:twfe} command varlist
[{help if}]
[{help in}]
[{help using}]
[{cmd:,}
{it:options}]

{pstd}
Command contains the type of regression. For example you could choose to run: reg, ivreg or sureg. The varlist contains the dependent followed by the independent variables. The syntaxis of the varlist replicates the syntaxis of the original cmd. {p_end}

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt {help using}}} followed by a path will save a set of matrices in that path. If using option is omitted the matrices will be stored in eresults {p_end}

{syntab:Required in the first regression}
{synopt:{opt abs:orb(absvars weight)}} categorical variables that indetify the fixed effects to be absorbed and the weight {p_end}
{synopt:{opt gen:erate(newvars)}} to create new group identifiers {p_end}
{synopt:{opt newv:ars(name)}} to create residualized variables with the prefix in "newvars" {p_end}
{synopt:{opt replace}} replace the variables for their residualized version {p_end}

{syntab:Required for the next regressions}
{synopt:{opt noproj}} Run the regression without creating any arrays or projecting new variables{p_end}

{syntab:SE/Robust}
{p2coldent:+ {opt vce}{cmd:(}{help vcetype}{cmd:)}}{it:vcetype}
may be {opt un:adjusted} (default), {opt robust} or {opt cluster} {clustervar}{p_end}



{synoptline}
{p2colreset}{...}
{p 4 6 2}

{marker description}{...}
{title:Description}
{pstd}
{cmd:twfe} is an algorithm to estimate the two-way fixed effect linear model. The algorithm relies on the Frisch-Waugh-Lovell theorem and applies to ordinary least squares (OLS), two-stage least squares (TSLS) and GMM estimators. {p_end}


{marker options}{...}
{title:Options}

{marker opt_model}{...}
{dlgtab: Full Description}

{phang}
{help using} adding using followed by a path will save a set of matrices in that path. Otherwise, the matrices will be saved in eresults. STATA may fail if it tries to store a matrix that exceeds matsize in ereturn. Adding using avoids that limitation.  {p_end}

{phang}
{opth abs:orb(absvars weight)}the two levels of fixed effects are specified in this option, it takes at least two variables with group identifiers. The optional third input can be a variable of analytic weights. {p_end}

{phang}
{opth gen:erate(newvars)} if the group identifiers are not consecutive after dropping redundants and missing observations, then you have this option to create new group identifiers. This option will create new variables in the database. {p_end}

{phang}
{opth newv:ars(name)} create new variables, the residualized variables will be stored as new variables that will be named with a prefix specified in `newvars'. This option will create new variables in the database.  {p_end}

{phang}
{opt replace} replace the variables for their residualized version. This option will re-write the database. {p_end}

{phang}
{opt noproj} run the regression without creating any arrays or projecting new variables. This option wil take the matrices created and stored in eresults or in the path selected by using option. {p_end}


{marker examples}{...}
{title:Examples}
{pstd}twfe reg y x1 x2 x3 x4, absorb(hhid tid w) newv(w_) vce(robust)  {p_end}
{pstd}twfe reg w_y w_x1, noproj vce(cluster hhid)  {p_end}
{pstd}twfe ivregress 2sls w_y w_x1 (w_x2= w_x3), noproj vce(robust)  {p_end}


{title:Possible Pitfalls and Common Mistakes}
{p2col 8 12 12 2: 1.} The algorithm will not work if the fixed effects after the cleaning of the redundants and missing observations are not consecutives and the generate option is not used. {p_end}
{p2col 8 12 12 2: 2.} The algorithm will not work if the user tries to create with "generate" or with "newvars" new variables with a name that already exits. {p_end}
{p2col 8 12 12 2: 3.} vce option is not allowed if the command selected is "sureg". {p_end}
{p2col 8 12 12 2: 4.} It will be an error if the user has in the database a variable with a name of a command such as reg (commonly used for variables of region).{p_end}
{p2col 8 12 12 2: 5.} Factor-variable and time-series operators not allowed.{p_end}

{synoptline}
{p2colreset}{...}
{p 4 6 2}

{title:Advanced}
{pstd}
Sometimes the same set of group identifiers have to be used to run several specifications. In these cases, one can save computation time by reusing some of the previous computations. The second and more advanced way to interact with the algorithm is  to  exploit this feature by using the functions "twset"(uses the group identifiers and weights to create a set of matrices that can be used for multiple specifications), "twres"(performs the second step of the algorithm) and "twest"(performs the last step of the algorithm, making the regression of the projected variables).

{pstd}
Also if "twset" is already ran the matrixes can be saved in a folder to not re-calculate them. The comand "twsave" will complete this task and the command "twload" load the matrixes saved previously. With this two commands the user can close STATA and still does not have to calculate them again.


{marker syntax}{...}
{title:Syntax}
{p 8 17 2}
{cmdab:twset} FixedEffects 
[{help weight}]
[{help if}]
[{help in}]
[{help using}]
[{cmd:,}
{it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt {help using}}} followed by a path will save a set of matrices in that path. If using option is omitted the matrices will be stored in eresults {p_end}
{synopt:{opt gen:erate(newvars)}} to create new group identifiers {p_end}

{marker examples}{...}
{title:Examples}
{pstd}twset hhid tid, gen(newids newts)  {p_end}
{pstd}twset hhid tid  {p_end}
{pstd}twset hhid tid w, gen(newids newts)  {p_end}
{pstd}twset hhid tid w  {p_end}
{pstd}twset hhid tid using "../folder/x", gen(newids newts)  {p_end}
{pstd}twset hhid tid using "../folder/x"  {p_end}

{title:Possible Pitfalls and Common Mistakes}
{p2col 8 12 12 2: 1.} The algorithm will not work if the fixed effects after the cleaning of the redundants and missing observations are not consecutives and the {opth gen:erate}  option is not used. {p_end}
{p2col 8 12 12 2: 2.} In this command the user has to use the if option to discard the missing observations of the variables that will be used in the next steps.{p_end}
{p2col 8 12 12 2: 3.} generate option: The algorithm wil break if the user tries to create a new variable with a name that already exists.{p_end}


{marker syntax}{...}
{title:Syntax}
{p 8 17 2}
{cmdab:twres} varlist
[{help using}]
[{cmd:,}
{it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt {help using}}} followed by a path will use a set of matrices saved in that path. If using option is omitted the command will use the arrays stored in eresults {p_end}
{synopt:{opt p:refix(name)}} creates new variables that are the residual projection of the originals {p_end}
{synopt:{opt replace}} replace the existing variables by their residual projection {p_end}

{title:Possible Pitfalls and Common Mistakes}
{p2col 8 12 12 2: 1.} newvars option:It will be an error if the use tries to create a variable that already exists.{p_end}
{p2col 8 12 12 2: 2.} If using option is omitted: it will an error if the user runs twset or twest and then another command since projvar needs the eresults from twset or twest.{p_end}
{p2col 8 12 12 2: 3.} If some of the variables has missing values the command will break.{p_end}
{p2col 8 12 12 2: 4.} Factor-variable and time-series operators not allowed.{p_end}


{marker examples}{...}
{title:Examples}
{pstd}twres y x1 x2, p(w_)  {p_end}
{pstd}twres y x1 x2, replace  {p_end}
{pstd}twres y x1 x2 using "../folder/x", p(w_)  {p_end}


{marker syntax}{...}
{title:Syntax}
{p 8 17 2}
{cmdab:twest} command varlist
[{cmd:,}
{it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:SE/Robust}
{p2coldent:+ {opt vce}{cmd:(}{help vcetype}{cmd:)}}{it:vcetype}
may be {opt un:adjusted} (default), {opt robust} or {opt cluster} {clustervar}{p_end}

{title:Possible Pitfalls and Common Mistakes}
{p2col 8 12 12 2: 1.}If using option is omitted: it will an error if the user runs twres and then another command since projvar needs the eresults from twres.{p_end}
{p2col 8 12 12 2: 2.} It will be an error if the use tries to use vce option will sureg command.{p_end}
{p2col 8 12 12 2: 3.} Factor-variable and time-series operators not allowed.{p_end}

{marker examples}{...}
{title:Examples}
{pstd}twest reg w_y w_x*, vce(cluster hhid)  {p_end}
{pstd}twest ivregress w_y w_x1 w_x2 (w_x3= w_x4 w_x5), vce(robust)  {p_end}

{marker syntax}{...}
{title:Syntax}
{p 8 17 2}
{cmdab:twsave} 
{cmd:}
[{help using}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt {help using}}} followed by a path will save a set of matrices in that path. If using option is omitted the matrices will be stored in the current directory {p_end}

{title:Possible Pitfalls and Common Mistakes}
{p2col 8 12 12 2: 1.} This command only work if "twset" has been ran previously and the e() has not been re-written. 

{marker examples}{...}
{title:Examples}
{pstd}twsave using "../folder/x" {p_end}
{pstd} In this example the matrices will be saved in "folder" with prefix "x".


{marker syntax}{...}
{title:Syntax}
{p 8 17 2}
{cmdab:twload} 
{cmd:}
[{help using}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt {help using}}} followed by a path will load the set of matrices that has been saved in that path. If using option is omitted command will search for matrices in the current directory. {p_end}

{title:Possible Pitfalls and Common Mistakes}
{p2col 8 12 12 2: 1.} This command will work only if there are matrices saved in a folder. 

{marker examples}{...}
{title:Examples}
{pstd}twload using "../folder/x" {p_end}
{pstd} In this example the matrices will be load from "folder" with prefix "x".

{synoptline}
{p2colreset}{...}
{p 4 6 2}



{title:Stored results}
{pstd}
{it:Note: In ereturn will be all the eresults of the original command selected plus the matrices and scalars needed for the algorithm.} 

{synoptset 15 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:e(dof_adj)}} Degree of freedom adjustment {p_end}
{synopt:{cmd:e(dimN)}} number of first fixed effect without the redundants observations{p_end}
{synopt:{cmd:e(dimT)}} number of second fixed effect without the redundants observations{p_end}
{synopt:{cmd:e(nested_adj)}} Adjustment in the degree of freedom if one of the fixed effect or both of them are nested within the cluster var {p_end}
{synopt:{cmd:e(rank_adj)}} Adjustment based on the rank of {cmd:e(A)} if {cmd:e(dimN)}<{cmd:e(dimT)}, otherwise is based on the rank of {cmd:e(C)}{p_end}

{synoptset 15 tabbed}{...}
{syntab:Macros}
{synopt:{cmd:e(absorbed)}} name of the fixed effects absorbed and the weight variable {p_end}

{synoptset 15 tabbed}{...}
{syntab:Matrices}
{synopt:{cmd:e(b)}} coefficient vector{p_end}
{synopt:{cmd:e(V)}} variance-covariance matrix adjusted of the estimators {p_end}


{p 8 8 1}
If "using" is omitted{p_end}
{p 8 8 1}{cmd:e(invDD)} is a vector of (dimN)x1 taking the diagonal elements of the inverse of the matrix D'D{p_end}
{p 8 8 1}{cmd:e(invHH)} is a vector of (dimT-1)x1 taking the diagonal elements of the inverse of the matrix H'H {p_end}
{p 8 8 1}{cmd:e(B)} Matrix B needed for the first and second step of the algorithm {p_end}


{p 12 8 1}If {cmd:e(dimN)}<{cmd:e(dimT)}{p_end}
{p 12 8 1}{cmd:e(A)} Matrix A needed for the first and second step of the algorithm{p_end}
{p 12 8 1}{cmd:e(CinvHHDH)} Matrix CinvHHDH needed for the first and second step of the algorithm{p_end}

{p 12 8 1}
If {cmd:e(dimN)}>={cmd:e(dimT)}{p_end}
{p 12 8 1}{cmd:e(C)} Matrix C needed for the first and second step of the algorithm{p_end}
{p 12 8 1}{cmd:e(AinvDDDH)} Matrix AinvDDDH needed for the first step of the algorithm{p_end}


{synoptset 15 tabbed}{...}
{syntab:Functions}
{synopt:{cmd:e(sample)}} marks estimation sample{p_end}
{p2colreset}{...}


{title:References}
{phang}
Somaini, P. and F.A. Wolak, (2016), An Algorithm to Estimate the Two-Way Fixed Effects Model, Journal of Econometric Methods, 5, issue 1, p. 143-152.

{title: Aditional References}
{phang}
Arellano, M. (1987), Computing Robust Standard Errors for Within-Groups Estimators, Oxford Bulletin of Economics and Statistics, 49, issue 4, p. 431–434.

{phang}
Cameron, A. C., & Miller, D. L. (2015). A practitioner’s guide to cluster-robust inference. Journal of human resources, 50(2), 317-372.



