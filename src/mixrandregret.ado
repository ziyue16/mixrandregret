*! mixrandregret 1.1.0 4Mar2022
*! [aut & dev] Ziyue Zhu  &  Álvaro A. Gutiérrez-Vargas

*  1.1.0:  	-mixrandregret- allows robust and cluster SEs and weights

program mixrandregret
	version 11
	if replay() {
		if (`"`e(cmd)'"' != "mixrandregret") error 301
		Replay `0'
	}
	else	Estimate `0'
end


program Estimate, eclass sortpreserve
	syntax varlist [if] [in] 		///
		[fweight pweight iweight/],	///
		GRoup(varname) 				///
		ALTernatives(varname)       ///
		RAND(string) [				/// random variable list
		BASEalternative(string)	    /// set base alternative specific constants (ASC) 
		NOCONStant                  /// suppress the ASC
		ID(varname) 				///
		CLuster(varname)			///
		Robust						///
		LN(integer 0) 				/// specifies that the last # variables in rand() have lognormally coefficients
		CORR						/// random coefficients are correlated
		NREP(integer 50)			/// the number of Halton draws used for the simulation
		BURN(integer 15)			/// the number of initial sequence elements to drop
		COLL						/// overrides checks for multicollinearity		
		Level(integer `c(level)')	/// set confidence level; default is level(95)
		USERdraws					///		
		FRom(string)				/// maximize options
        TRace						///
		GRADient					///
		HESSian						///
		SHOWSTEP					///
		ITERate(passthru)			///
		TOLerance(passthru)			///
		LTOLerance(passthru)		///
		GTOLerance(passthru)		///
		NRTOLerance(passthru)		///
		CONSTraints(passthru)		///
		TECHnique(passthru)			///
		DIFficult					///
	]
	
		local mlopts `trace' `gradient' `hessian' `showstep' `iterate' `tolerance' ///
	    `ltolerance' `gtolerance' `nrtolerance' `technique' `difficult'
		
	** globals to mata **
	global group_mata =  "`group'"
	global alternatives_mata =  "`alternatives'"
	global cluster_mata =  "`cluster'"
	
	capture mata: mata drop RRM_log() mixRRM_gf0()
	findfile "mixRRM_gf0.mata"
	do "`r(fn)'"
	
	/*=======================================================================*/	
	/*===================       CHECK BASICS        =========================*/	
	/*=======================================================================*/	
    
	** Don't need technique checks **
	** Check that group, id, cluster and alternative variables are numeric **
	capture confirm numeric var `group'
		if _rc != 0 {
			di in r "The group variable must be numeric"
			exit 498
		}
	
	capture confirm numeric var `alternatives'
	if _rc != 0 {
	    di in r "The alternative variable must be numeric"
		exit 498
	}
		
	if ("`id'" != "") {
		capture confirm numeric var `id'
		if _rc != 0 {
			di in r "The id variable must be numeric"
			exit 498
		}
	}
		
	if ("`cluster'" != "") {
		capture confirm numeric var `cluster'
		if _rc != 0 {
			di in r "The cluster variable must be numeric"
			exit 498
		}
	}

	
	** Mark the estimation sample **
	marksample touse
	markout `touse' `group' `alternatives' `rand' `id' `cluster'

	** Create locals for later: 
	gettoken lhs fixed : varlist
	
	** Create the rhs using the fixed variables + the random ones
	local rhs `fixed' `rand'
		* lhs  : variable with alternatives selected by individuals
		* fixed: fixed variables 
		* rand : random variables
		* rhs  : fixed + rand variables (this is used later to define the equations)

	** Prevent user from using same variable as random and fix **
	local k1 : word count `fixed'
	local k2 : word count `rand' 
	forvalues i = 1(1)`k1' {
		forvalues j = 1(1)`k2' {
			local w1 : word `i' of `fixed' 
			local w2 : word `j' of `rand'
			if ("`w1'" == "`w2'") {
				di in red "The variable `w1' is specified to have both fixed and random coefficients"
				exit 498
			} 	
		}
	}

	/* robust, wgt, coll and constraints setting, safe to omit? */
	
	** Check for multicollinearity **
	local rhs `fixed' `rand'
	if ("`coll'" == "") {
		qui _rmcoll `rhs' 
		if ("`r(varlist)'" != "`rhs'" & "`constraints'" != "") {
			di in gr "Some variables are collinear - make sure this is intended, i.e. because you are"
			di in gr "estimating an error-components model with the necessary constraints imposed"
		}
		if ("`r(varlist)'" != "`rhs'" & "`constraints'" == "") {
			di in red "Some variables are collinear - check your model specification"
			exit 498
		}
	}
	
	/* We don't estimate random regret model to get inital values in mixed version. ///
	   Users are advised to estimate on their own to make comparsion */
	   
	** Drop missing data **
	preserve
	qui keep if `touse'
	
	** Check that the independent variables vary within groups **
	sort `group'
	foreach var of varlist `rhs' {
		capture by `group': assert `var'==`var'[1]
		if (_rc == 0) {
			di in red "Variable `var' has no within-group variance"
			exit 459		
		}
	}
	
	** Check that the dependent variable only takes values 0-1 **
	capture assert `lhs' == 0 | `lhs' == 1
	if (_rc != 0) {
		di in red "The dependent variable must be a 0-1 variable indicating which alternatives are chosen"
		exit 450		
	}
	
	** Check that each group has only one chosen alternative **
	tempvar chonum
	sort `group'
	qui by `group': egen `chonum' = sum(`lhs')
	capture assert `chonum' == 1
	if (_rc != 0) {
		di in red "At least one group has either more than one chosen alternative or none"
		exit 498		
	}
	
	*** ASC Checks (to be added) ***		
		
	** Generate individual id **
	if ("`id'" != "") {
		tempvar nchoice pid
		sort `group'
		by `group': gen `nchoice' = cond(_n==_N,1,0)
		sort `id'
		by `id': egen `pid' = sum(`nchoice')		
		qui duplicates report `id'
		mata: mixr_np = st_numscalar("r(unique_value)")
		mata: mixr_T = st_data(., st_local("pid"))
	}
	else {
		qui duplicates report `group'
		mata: mixr_np = st_numscalar("r(unique_value)")
		mata: mixr_T = J(st_nobs(),1,1)
	}

	** Generate choice occasion id **
	tempvar csid
	sort `group'
	by `group': egen `csid' = sum(1)
	
	** Sort data **
	sort `id' `group'
	
	
	/*=======================================================================*/	
	/*=================== MATA MATRICES AND SCALARS =========================*/	
	/*=======================================================================*/	

	** Set Mata matrices and scalars to be used in optimisation routine **
	local kfix: word count `fixed'
	local krnd: word count `rand'

	mata: mixr_X = st_data(., tokens(st_local("rhs")))
	mata: mixr_Y = st_data(., st_local("lhs"))
	mata: mixr_CSID = st_data(., st_local("csid"))
	mata: mixr_IND = st_data(., st_local("id"))


	mata: mixr_nrep = strtoreal(st_local("nrep"))
	mata: mixr_kfix = strtoreal(st_local("kfix"))
	mata: mixr_krnd = strtoreal(st_local("krnd"))
	mata: mixr_krln = strtoreal(st_local("ln"))
	mata: mixr_burn = strtoreal(st_local("burn"))
	
    if ("`userdraws'" != "") mata: mixr_user = 1	
	else mata: mixr_user = 0

	
	/*=======================================================================*/	
	/*===================       LOCALS FOR ML       =========================*/	
	/*=======================================================================*/	

	** Create macro to define equations for optimisation routine **
	local mean (Mean: `lhs' = `rhs', noconst)
	if ("`corr'" == "") {
		mata: mixr_corr = 0
		local sd (SD: `rand', noconst)
		local max `mean' `sd' 
	}

	** not setting up starting values **
	
	
	/*=======================================================================*/	
	/*===================             ML            =========================*/	
	/*=======================================================================*/		

	** Run optimisation routine **
	ml model gf0  mixRRM_gf0() `max', maximize `mlopts' 

    * To be returned as e() **
    ereturn local title "Mixed random regret model"
    ereturn local indepvars `rhs'
    ereturn local depvar `lhs'
    ereturn local group `group'
    
	ereturn scalar kfix = `kfix'
    ereturn scalar krnd = `krnd'
	ereturn scalar krln = `ln'
	ereturn scalar nrep = `nrep'
	ereturn scalar burn = `burn'
	
		
	if ("`id'" != "") ereturn local id `id'
	if ("`robust'" != "") ereturn local vcetype Robust
	
	if ("`userdraws'" != "") ereturn scalar userdraws = 1
	else ereturn scalar userdraws = 0
	
	ereturn local cmd "mixrandregret"
	
	Header
	Replay, level(`level')

end


program Header
     loc df_m = e(df_m)
	 di ""
     di in gr e(title)
     di ""
     di in gr `"Case ID variable: `=abbrev("${group_mata}",24)'"' _col(48) /// Case ID
     "Number of cases" _col(67) "= " in ye %10.0g e(n_cases)
    
	di in gr `"Alternative variable: `=abbrev("${alternatives_mata}",24)'"' _col(48)
    //  "Number of obs" _col(67) "= " in ye %10.0g e(N)       
end


program Replay
	syntax [, Level(integer `c(level)') CORR]
	ml display , level(`level')

end


