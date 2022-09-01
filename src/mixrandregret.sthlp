{smcl}
{* *! version 1.0.0  21apr2022}{...}
{title:Title}

{p 5 15 2}{hi:mixrandregret} {hline 2} Mixed Random Regret Minimization Model


{title:Syntax}

{p 8 15 2}
{cmd:mixrandregret}
{depvar}
[{indepvars}] {ifin} {weight} {cmd:,}
{cmdab:gr:oup(}{varname}{cmd:)}
{cmdab:rand:(}{varlist}{cmd:)}
{cmdab:alt:ernatives(}{varname}{cmd:)}
[{opth base:alternative(#)}
 {opt nocons:tant}
 {opth id(varname)}
 {opth cl:uster(varname)}
 {opt r:obust}
 {opth ln(#)}
 {opth nrep(#)}
 {opth burn(#)}
 {opth l:evel(#)}
 {it:maximize_options}]
 
 
{p 8 15 2}
{cmd:mixrpred}
{newvar} {ifin}
[{cmd:,} 
 {opt proba}
 {opt xb}
 {opth nrep(#)} 
 {opth burn(#)}]


{p 8 15 2}
{cmd:mixrbeta}
{varlist} {ifin}
{cmd:,} {opth sav:ing(filename)} 
[{opt plot}
 {opth nrep(#)}
 {opth burn(#)}
 {opt replace}]
 
{phang}
{opt depvar} equal to 1 identifies the chosen alternatives,
whereas a 0 indicates the alternatives that were not chosen.
There can be only one chosen alternative for each case.

{phang}
{opt fweight}s, {opt iweight}s, and {opt pweight}s are allowed (see {help weight}), 
but they are interpreted to apply to decision-makers, not to individual observations.{p_end}


{title:Description}

{pstd}
{cmd:mixrandregret} utilizes the mixed random regret minimization model described in {help mixrandregret##hensher2016:Hensher et al. (2016)}, which is a mixed version of the classic random regret minimization model introduced in {help mixrandregret##chorus2010:Chorus. C. (2010)}. {cmd:mixrandregret} extends the {cmd:randregret} ({help mixrandregret##gutierrez2021:Gutiérrez-Vargas et al, 2021}) and allows to specify normal and log-normally distributed taste parameters inside the regret function. The command uses maximum simulated likelihood for estimation ({help mixrandregret##train2003:Train. K., 2003}). 

{pstd}
{cmd:mixrpred} can be used following {cmd:mixrandregret} to obtain the predicted probabilities or the predicted systematic regret. 

{pstd}
{cmd:mixrbeta} can be used following {cmd:mixrandregret} to calculate
individual-level parameters corresponding to the variables in the specified
{it:varname} using the method proposed by {help mixrandregret##train2003:Train. K. (2003)} chap. 11.  The individual-level parameters are stored in a data
file specified by the user.


{title:Options}

{marker classicoptions}{...}
{synoptset 23 tabbed}{...}
{synopthdr :mixrandregret}
{synoptline}
{syntab:Model}
{p2coldent :* {opth gr:oup(varname)}} is required and specifies a numeric identifier variable ({it:varname}) for the choice occasions.{p_end}
{p2coldent :* {opth rand:(varlist)}} is required and specifies the independent variables whose coefficients are random. The random coefficients can be specified to be normally or lognormally distributed (see the ln() option). The variables immediately following the dependent variable in the syntax are specified to have fixed coefficients.{p_end}
{p2coldent :* {opth alt:ernatives(varname)}} use {it:varname} to identify the
        alternatives available for each case.{p_end}
{synopt:{opth base:alternative(#)}} sets base Alternative Specific Constants (ASC).{p_end}
{synopt:{opt nocons:ant}} suppress the alternative specific constants.{p_end}
{synopt:{opth id(varname)}} specifies a numeric identifier variable for the decision
makers. This option should be specified only when each individual performs
several choices; i.e., the dataset is a panel.{p_end}
{synopt :{opt r:obust}, {opt cl:uster}} see {help estimation options}. The cluster variable must be numeric.{p_end}
{synopt:{opth ln(#)}} specifies that the last {it:#} variables in {opt rand()} have lognormally rather than normally distributed coefficients. The default is {cmd:ln(0)}.{p_end}

{syntab:Simulation}
{synopt :{opth nrep(#)}} specifies the number of Halton draws used for the simulation. The default is {cmd:nrep(50)}.{p_end}
{synopt :{opth burn(#)}} specifies the number of initial sequence elements to drop when creating the Halton sequences. The default is {cmd:burn(15)}. Specifying this option helps reduce the correlation between the sequences in each dimension.{p_end}

{syntab:Reporting}
{synopt :{opth l:evel(#)}} set confidence level; default is {cmd:level(95)}.{p_end}

{syntab:Maximization}
{synopt:{it:maximize_options}}
{opt dif:ficult},
{opt tech:nique(algorithm_spec)}, 
{opt iter:ate(#)}, {opt tr:ace}, {opt grad:ient}, 
{opt showstep}, {opt hess:ian}, {opt tol:erance(#)}, 
{opt ltol:erance(#)} {opt gtol:erance(#)}, {opt nrtol:erance(#)}, 
{opt from(init_specs)}; see {helpb maximize}.{p_end}
{synoptline}
{p2colreset}{...}

{marker geneoptions}{...}
{synoptset 23 tabbed}{...}
{synopthdr :mixrpred}
{synoptline}
{syntab:Model}
{synopt : {opt proba}} calculate probability of a positive outcome; the default{p_end}
{synopt : {opt xb}} calculate linear prediction of the systematic regret{p_end}
{synopt :{opth nrep(#)}} specifies the number of Halton draws used for the simulation. The default is {cmd:nrep(50)}.{p_end}
{synopt :{opth burn(#)}} specifies the number of initial sequence elements to drop when creating the Halton sequences. The default is {cmd:burn(15)}. Specifying this option helps reduce the correlation between the sequences in each dimension.{p_end}
{synoptline}
{p2colreset}{...}

{marker muoptions}{...}
{synoptset 23 tabbed}{...}
{synopthdr :mixrbeta}
{synoptline}
{syntab:Model}
{p2coldent :* {opth saving(filename)}} saves individual-level parameters to {it:filename}{p_end}
{synopt:{opt plot}} save conditional distribution graphs (histogram and Kdensity) for betas.{p_end}
{synopt :{opth nrep(#)}} specifies the number of Halton draws used for the simulation. The default is {cmd:nrep(50)}.{p_end}
{synopt :{opth burn(#)}} specifies the number of initial sequence elements to drop when creating the Halton sequences. The default is {cmd:burn(15)}. Specifying this option helps reduce the correlation between the sequences in each dimension.{p_end}
{synopt :{opt replace}} overwrites {it:filename}.{p_end}
{synoptline}
{p2colreset}{...}


{title:Examples}

{pstd}
Consider the following toy example that contains the data for two individuals who makes two choices. On each choice occasion, the individual has three alternatives. {cmd:choice} is the dependent variable, and {cmd:x_rnd} and {cmd:x_fix} are the independent variables or alternative attributes:

{cmd}
        id_ind id_cs choice_set alt  x_rnd   x_fix  choice
          1      1       1       1   -4.17   -0.83    0
          1      1       1       2    1.39   -0.87    0
          1      1       1       3    2.63   -0.37    1
          1      2       2       1   -0.45    1.12    1
          1      2       2       2    2.98    1.07    0
          1      2       2       3    1.37    0.75    0
          2      3       1       1    2.88   -0.25    0
          2      3       1       2   -3.12   -0.52    1
          2      3       1       3   -1.44    0.39    0
          2      4       2       1   -0.56    0.51    0
          2      4       2       2   -1.03    0.02    0
          2      4       2       3    0.99   -2.22    1{txt}
      
{pstd}
A mixed random regret model where {cmd:x_rnd} has a normally distributed coefficient and {cmd:x_fix} has a fixed coefficient, and ASCs are suppressed can be specified as follows:

{phang2}{cmd:. mixrandregret choice x_fix, group(id_cs) id(id_ind) rand(x_rnd) alternatives(alt) nocons}{p_end}

{pstd}
A model where {cmd:x_rnd} has a lognormally distributed coefficient can be specified as follows:

{phang2}{cmd:. mixrandregret choice x_fix, group(id_cs) id(id_ind) rand(x_rnd) alternatives(alt) ln(1) nocons}{p_end}

{pstd}
Use alternative 1 as base alternative to calculate ASC can be specified as follows:

{phang2}{cmd:. mixrandregret choice x_fix, group(id_cs) id(id_ind) rand(x_rnd) alternatives(alt) ln(1) basealternative(1)}{p_end}

{pstd}
Warning: initial values not provided. We highly advise the users to use starting values from {cmd:randregret} estimates for the means of the distributions of the random coefficients.

{pstd}
An example of setting starting values (2 parameters from fixed effects, 2 parameters from random effects - mean and SD, and 2 parameters from ASC):

{phang2}{cmd:. mixrandregret choice x_fix, group(id_cs) id(id_ind) rand(x_rnd) alternatives(alt) ln(1) basealternative(1) from(1 1 1 1 1 1)}{p_end}

{pstd}
Predict choice probabilities and save values in a variable named {cmd:p}:

{phang2}{cmd:. mixrpred p, proba}{p_end}

{pstd}
Linear prediction of the systematic regret and save values in a variable named {cmd:regret}:

{phang2}{cmd:. mixrpred regret, xb}{p_end}

{pstd}
Obtain the individual-level parameters, generate the plot, and save in the file named {cmd:ind_rnd}:

{phang2}{cmd:. mixrbeta x_rnd, plot saving(ind_rnd) replace}{p_end}


{marker author}{...}
{title:Authors}

{pstd}
Ziyue Zhu{break}
Faculty of Science{break}
KU Leuven{break} 
Leuven, Belgium{break}
ziyue.zhu@student.kuleuven.be

{pstd}Álvaro A. Gutiérrez Vargas{break}
Faculty of Economics and Business{break}
KU Leuven{break} 
Research Centre for Operations Research and Statistics (ORSTAT){break}
Leuven, Belgium{break}
alvaro.gutierrezvargas@kuleuven.be

{pstd}
Martina Vandebroek{break}
Faculty of Economics and Business{break}
KU Leuven{break} 
Research Centre for Operations Research and Statistics (ORSTAT){break}
Leuven, Belgium{break}
martina.vandebroek@kuleuven.be>


{title:References}

{marker hensher2016}{...}
{phang} Hensher, David A and Greene, William H and Ho, Chinh Q. 2016
{browse "https://ascelibrary.org/doi/abs/10.1061/%28ASCE%29TE.1943-5436.0000827":Random regret minimization and random utility maximization in the presence of preference heterogeneity: an empirical contrast}.
{it:Journal of Transportation Engineering} Volume 142 Number 4: pp. 04016009.

{marker chorus2010}{...}
{phang}Chorus. C. 2010.
{browse "https://ojs-lib2.tudelft.nl/ejtir/article/view/2881":A New Model of Random Regret Minimization}.
{it:European Journal of Transport and Infrastructure Research} 10: pp. 181-196.

{marker gutierrez2021}{...}
{phang}Gutiérrez-Vargas. Á., Meulders. M., and Vandebroek. M. 2021.
{browse "https://www.stata-journal.com/article.html?article=st0649":randregret: A command for fitting random regret minimization models using Stata}.
{it:The Stata Journal} Volume 21 Number 3: pp. 626-658.

{marker train2003}{...}
{phang}Train. K. 2003.
{browse "https://eml.berkeley.edu/books/choice2.html":Discrete Choice Methods with Simulation}.
{it:Cambridge University Press}.

