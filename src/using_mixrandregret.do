*! 1.1.1 17 Jan 2022
* Author: Ziyue ZHU (ziyue.zhu@student.kuleuven.be) 

/*==================================*/
/*====      Setting up        ======*/
/*==================================*/
clear all
/*You need to adjust this and change it for the folder where you will be developing the command.*/
global route = "/Users/zhuziyue/Documents/Master (KU Leuven)/Thesis/wrapping_command_102"
cd "$route"
/*This will search for the adofile of the command.*/
findfile mixrandregret.ado



set seed 777
local J  3   /*Number of alternatives*/
local N  20 /*Number of Individuals*/
local t  5  /*Number of choice sets per individual*/

/*===================================================*/
/*===  Define the inputs for mata evaluator    ======*/
/*===================================================*/

/*Star defining the inputs for the mata evaluator defined latter.*/
local rhs "x1 x2"   // right hand side variables: defined by hand
local lhs "choice"  // left hand side variables: defined by hand
local group "id_cs" // variable with choice sets
local id "id_ind"   // individual ID

local nrep 2       // # of halton draws per random variable
local kfix 0        // # of fixed variable (0 in this script)
local krnd 2        // # of random variables (both 2 variables, x1 x2)
local burn 15       // # of burning draws in halton 

/*Set total number of observations on the data*/
set obs  `=`N'*`t'*`J''
* 20 individuals choose from 25 choice sets, and each choice set includes 3 alternatives
*ssc install seq 

/* id_ind = id for each individual */
seq id_ind, b(`=`J'*`t'')

/* id_cs_of_ind_n = id for each choice situation for each individual n */
seq id_cs_of_ind_n,  f(1) t(`=`J'') b(`=`J'')

/* id_cs = id choice situations without considering the individual */
seq id_cs,  f(1) t(`=`J'*`N'*`t'') b(`=`J'')
sort id_cs

/* alternative = size of the choice set from where the individual choose*/
seq alternative, t(`J')


local s = 3
/*two generic attributes (random)*/
gen x1 =  runiform(-`s',`s')
gen x2 =  runiform(-`s',`s')


/*create a fix attribute */
gen x_fix = rnormal(0,1)


/* parameters for distribution */
local mu_1  = 1
local mu_2  = 2

local sigma_1  = 0.5
local sigma_2  = 0.5


mata: b1 = `mu_1' :+ sort(J(`=`J'*`t'', 1, rnormal(`N',1,0,1)),1) :* `sigma_1'
mata: b2 = `mu_2' :+ sort(J(`=`J'*`t'', 1, rnormal(`N',1,0,1)),1) :* `sigma_2'

mata: b = b1 , b2




global individuals = "id_ind"
global choice_sets = "id_cs"
sort id_cs


/*===================================================*/
/*====             DATA Generation             ======*/
/* Two (independent) normally distributed parameters */
/*===================================================*/

mata:
/* program a function to calculate observed regret */
function RRM_log_sim_data(real matrix x_n, 
						  real rowvector betas)
{
  real scalar i, j
  real matrix regret_n 
  
  regret_n = J(rows(x_n), cols(x_n), 0)
  for(i=1; i <= rows(x_n); ++i) { 
	for(j=1; j <= rows(x_n); ++j) { 
		if (i!=j) { 
			r_i = ln(1 :+ exp( betas :* ( x_n[j , . ] :-  x_n[i, . ]))) 				
			regret_n[i, . ] = regret_n[i, . ] :+ r_i 		
			} 
		}	  
	}
return(regret_n)
}
		
// Generates a view of all attributes (M) x*
st_view(X = ., ., "x1 x2")

// Generates a view id choice situations
st_view(panvar_choice_sets = ., ., "id_cs")
		
// set up panel information where each panel unit refers to a choice set 
task_n = panelsetup(panvar_choice_sets, 1)

// # of choice sets 
npanels_choice_sets = panelstats(task_n)[1]
task_n
npanels_choice_sets

for(t=1; t<=npanels_choice_sets; t++) {


	// read in data rows pertaining to subject n & store in a matrix suffixed _n
	x_n = panelsubmatrix(X, t, task_n) 
	b_n = panelsubmatrix(b, t, task_n) 
	
	
	/*Observed (deterministic) Regret*/
	R_i = rowsum(RRM_log_sim_data(x_n, b_n[1,])) /*take first row of block of parameters*/

    /* Type 1 error */
	epsilon  = -1*log(-log(runiform(   rows(R_i),1,0,1)))
    /*Random regret */
	RR_i = R_i  :+ epsilon

    EXP_RR_i = exp(-RR_i)


	P_i   = EXP_RR_i :/ quadcolsum(EXP_RR_i, 1)

	choice_i =  (P_i	 :==max(P_i))
	// collect all choices of choice situations (t) of individual (n)
	  if (t==1)  choice = choice_i
	  else       choice = choice   \ choice_i

}
// Creates a new Stata variable called "choice"    
 idx = st_addvar("float", "choice")
// save the content of Y on "choice" Stata variable
 st_store(., idx, choice)

end


*** Certification lines ***

* The group variable must be numeric
gen group_string = "a"
rcof "noisily mixrandregret choice x_fix, group(group_string) id(id_ind) rand(x1 x2) alt(alter)" == 498

* The alternative variable must be numeric
gen alter_string = "b"
rcof "noisily mixrandregret choice x_fix, group(id_cs) id(id_ind) rand(x1 x2) alt(alter_string)" == 498

* The id variable must be numeric
gen id_string = "c"
rcof "noisily mixrandregret choice x_fix, group(id_cs) id(id_string) rand(x1 x2) alt(alternative)" == 498

* The cluster variable must be numeric
gen cluster_str = "d"
rcof "noisily mixrandregret choice x_fix, group(id_cs) id(id_ind) rand(x1 x2) alt(alternative) cluster(cluster_str)" == 498

* prevent user from using same variable as random and fix 
rcof "noisily mixrandregret choice x_fix x2 , group(id_cs) id(id_ind) rand(x1 x2) alt(alternative)" == 498

* multicollinearity
gen x_coll = x_fix
rcof "noisily mixrandregret choice x_fix x_coll, group(id_cs) id(id_ind) rand(x1 x2) alt(alternative)" == 498
*rcof "noisily mixrandregret choice x_fix x_coll, group(id_cs) id(id_ind) rand(x1 x2) alt(alternative) constr(x_coll = 1)" == 498


* independent variables vary within groups
gen x_fix_f = 1 if id_cs > 10
replace x_fix_f = 0 if id_cs <= 10

rcof "noisily mixrandregret choice x_fix_f, group(id_cs) id(id_ind) rand(x1 x2) alt(alternative)" == 459

* dependent variable only takes values 0-1
gen choice_f = 2
rcof "noisily mixrandregret choice_f x_fix, group(id_cs) id(id_ind) rand(x1 x2) alt(alternative)" == 450

@

*Calling the mixrandregret command *

mixrandregret choice x_fix, group(id_cs) id(id_ind) rand(x1 x2) alt(alternative) ln(1) iter(1)
