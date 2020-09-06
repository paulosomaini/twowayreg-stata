{smcl}
{* *! version 1.0  6 Sep 2020}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "Install command2" "ssc install command2"}{...}
{vieweralsosee "Help command2 (if installed)" "help command2"}{...}
{viewerjumpto "Syntax" "twowayreg##syntax"}{...}
{viewerjumpto "Description" "twowayreg##description"}{...}
{viewerjumpto "Options" "twowayreg##options"}{...}
{viewerjumpto "Remarks" "twowayreg##remarks"}{...}
{viewerjumpto "Examples" "twowayreg##examples"}{...}
{title:Title}
{phang}
{bf:twowayregwrap} {hline 2} Algorithm to efficiently estimate a two-way fixed effects model based on Somaini and Wolak(2016).

{marker syntax}{...}
{title:Syntax}
{p 8 17 2}
{cmdab:twowayregwrap} command varlist
[{help if}]
[{help in}]
[{help using}]
[{cmd:,}
{it:options}]

{marker description}{...}
{title:Varlist}
{pstd}
The syntaxis of the varlist replicates the syntaxis of the original command. 

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
{cmd:twowayregwrap} is an algorithm to estimate the two-way fixed effect linear model. The algorithm relies on the Frisch-Waugh-Lovell theorem and applies to ordinary least squares (OLS), two-stage least squares (TSLS) and GMM estimators. {p_end}


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
{pstd}twowayregwrap reg y x1 x2 x3 x4, absorb(hhid tid w) newv(w_) vce(robust)  {p_end}
{pstd}twowayregwrap reg w_y w_x1, noproj vce(cluster hhid)  {p_end}
{pstd}twowayregwrap ivregress 2sls w_y w_x1 (w_x2= w_x3), noproj vce(robust)  {p_end}

{synoptline}
{p2colreset}{...}
{p 4 6 2}

{title:Advanced}
{pstd}
Sometimes the same set of group identifiers have to be used to run several specifications. In these cases, one can save computation time by reusing some of the previous computations. The second and more advanced way to interact with the algorithm is to exploit this feature by using the functions "twowayset" (uses the group identifiers and weights to create a set of matrices that can be used for multiple specifications), "projvar"(performs the second step of the algorithm) and "twowayreg"(performs the last step of the algorithm, making the regression of the projected variables).{p_end}

{marker syntax}{...}
{title:Syntax}
{p 8 17 2}
{cmdab:twowayset} FixedEffects 
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
{pstd}twowayset hhid tid, gen(newids newts)  {p_end}
{pstd}twowayset hhid tid  {p_end}
{pstd}twowayset hhid tid w, gen(newids newts)  {p_end}
{pstd}twowayset hhid tid w  {p_end}
{pstd}twowayset hhid tid using "../folder/x", gen(newids newts)  {p_end}
{pstd}twowayset hhid tid using "../folder/x"  {p_end}


{marker syntax}{...}
{title:Syntax}
{p 8 17 2}
{cmdab:projvar} varlist
[{help using}]
[{cmd:,}
{it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt {help using}}} followed by a path will use a set of matrices saved in that path. If using option is omitted the command will use the arrays stored in eresults {p_end}
{synopt:{opt p:refix(name)}} creates new variables that are the residual projection of the originals {p_end}
{synopt:{opt replace}} replace the existing variables by their residual projection {p_end}


{marker examples}{...}
{title:Examples}
{pstd}projvar y x1 x2, p(w_)  {p_end}
{pstd}projvar y x1 x2, replace  {p_end}
{pstd}projvar y x1 x2 using "../folder/x", p(w_)  {p_end}


{marker syntax}{...}
{title:Syntax}
{p 8 17 2}
{cmdab:twowayreg} command varlist
[{cmd:,}
{it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:SE/Robust}
{p2coldent:+ {opt vce}{cmd:(}{help vcetype}{cmd:)}}{it:vcetype}
may be {opt un:adjusted} (default), {opt robust} or {opt cluster} {clustervar}{p_end}


{marker examples}{...}
{title:Examples}
{pstd}twowayreg reg w_y w_x*, vce(cluster hhid)  {p_end}
{pstd}twowayreg ivregress w_y w_x1 w_x2 (w_x3= w_x4 w_x5), vce(robust)  {p_end}

{synoptline}
{p2colreset}{...}
{p 4 6 2}

{title:Stored results}
In eresults results will be all the eresults of the original command selected plus the matrixes and scalars needed for the regressions.  {p_end}

If using option is omitted:
{synoptset 15 tabbed}{...}


{title:Author}
{pstd}Paulo Somaini{break}
Stanford GSB{break}
Email: {browse "mailto:soma@stanford.edu":soma@stanford.edu}
{p_end}

{title:References}
{phang}
Somaini, P. and F.A. Wolak, (2016), An Algorithm to Estimate the Two-Way Fixed Effects Model, Journal of Econometric Methods, 5, issue 1, p. 143-152.

{title: Aditional References}
{phang}
Arellano, M. (1987), Computing Robust Standard Errors for Within-Groups Estimators, Oxford Bulletin of Economics and Statistics, 49, issue 4, p. 431–434.

{phang}
Cameron, A. C., & Miller, D. L. (2015). A practitioner’s guide to cluster-robust inference. Journal of human resources, 50(2), 317-372.



