*! mixrbeta 1.1.0 11Apr2022
*! [aut & dev] Ziyue Zhu & Álvaro A. Gutiérrez-Vargas

program define mixrbeta
	version 9.2

	syntax varlist [if] [in], ///
	SAVing(string)            ///
	[PLOT                     ///
	NREP(integer 50)          ///
	BURN(integer 15)          ///
	REPLACE]
	
	** Pre-checks **
	if "`e(cmd)'"=="" exit 301 

	local cmd `e(cmd)'

	** Mark the prediction sample **
	marksample touse, novarlist
	markout `touse' `e(depvar)' `e(indepvars)' `e(group)' `e(id)' `e(basealternative)' `e(alternatives)'

	** Mark groups with no chosen alternatives due to missing data **
	tempvar cho
	sort `e(group)'
	qui by `e(group)': egen `cho' = sum(`e(depvar)'*`touse')
	qui replace `cho' = . if `cho' == 0
	markout `touse' `cho'
	
	** Drop data not in prediction sample **
	preserve
	qui keep if `touse'

	** Generate individual id **
	if ("`e(id)'" != "") {
		tempvar nchoice pid
		sort `e(group)'
		by `e(group)': gen `nchoice' = cond(_n == _N, 1, 0)
		sort `e(id)'
		by `e(id)': egen `pid' = sum(`nchoice')
		qui duplicates report `e(id)'
		local np = r(unique_value)
		mata: mixrbeta_np = st_numscalar("r(unique_value)")
		mata: mixrbeta_T = st_data(., st_local("pid"))
	}
	else {
		tempvar id
		clonevar `id' = `e(group)'
		qui duplicates report `e(group)'
		local np = r(unique_value)
		mata: mixrbeta_np = st_numscalar("r(unique_value)")
		mata: mixrbeta_T = J(st_nobs(), 1, 1)
	}

	** Generate dummy for last obs for each decision-maker **
	tempvar last
	if ("`e(id)'" != "") {
		by `e(id)': gen `last' = cond(_n == _N, 1, 0)
	}
	else {
		by `e(group)': gen `last' = cond(_n == _N, 1, 0)
	}

	** Generate choice occasion id **
	tempvar csid
	sort `e(group)'
	by `e(group)': egen `csid' = sum(1)
	qui duplicates report `e(group)'
	
	** Sort data **
	sort `e(id)' `e(group)'

	** Set Mata matrices to be used in prediction routine **
	local rhs `e(indepvars)'
	local lhs `e(depvar)'
	if ("`e(id)'" != "") local id `e(id)'
	else local id `e(group)'
	mata: mixrbeta_X = st_data(., tokens(st_local("rhs")))
	mata: mixrbeta_Y = st_data(., st_local("lhs"))
	mata: mixrbeta_CSID = st_data(., st_local("csid"))
	mata: mixrbeta_ID = st_data(., st_local("id"), st_local("last"))
	
	** Create dataset containing beta estimates **
	drop _all
	qui set obs `np'
	qui gen double `id' = .
	foreach var of local rhs {
		qui gen double `var' = .
	}
	mata: st_store(., st_local("id"), mixrbeta_ID)
	mata: st_view(mixrbeta_PB=., ., tokens(st_local("rhs")))
	
	** Restore data **
	restore
	preserve
	
	** Parsing parameter vector **
	tempname b_all b_beta ASC_beta
	matrix `b_all' = e(b)
	if "${cons_demanded}" =="NO" {
		matrix `b_beta' = `b_all'
	}
	else if "${cons_demanded}" =="YES" {
		qui tab `e(alternatives)'
		local n_altern = r(r)
		matrix `b_beta'   = `b_all'[1, 1..(`e(rank)'-`n_altern')+1]
		matrix `ASC_beta' = `b_all'[1, `e(rank)'-(`n_altern'-2)..`e(rank)']
	}
		
	** Generate ASC if needed **
	if "`e(ASC)'"=="YES" {
		// generate ASC as tempvars
		qui levelsof `e(alternatives)', local(levels_altern)
		tempvar ASC_
		foreach i of local levels_altern {
			tempvar ASC_`i'
			qui gen `ASC_'`i' = (`e(alternatives)' == `i')
		}
		qui tab `e(alternatives)'
		local n_altern = r(r)
		if "`e(basealternative)'" != "" {
			drop `ASC_'`e(basealternative)'
		}
		else { // drop the alternative with lower number
			qui sum `e(alternatives)' , meanonly
			local min_alt = r(min)
			qui drop `ASC_'`min_alt'
		}
		// Mata allocation of the ASC variables
		mata: st_view(ASC = ., ., "`ASC_'*")
	}
	else {
		// dummy ASC (=0 for conformability)
		mata: ASC = J(1, 1, 0)
		// coefficients ASC
		matrix `ASC_beta' =J(1, 1, 0)
	}
		
	mata: b_beta= st_matrix("`b_beta'")
	mata: b_all= st_matrix("`b_all'")
	mata: ASC_beta= st_matrix("`ASC_beta'")
	
	mata: mixr_beta("`b_beta'")
	
	keep `id' `varlist'

	if ("`replace'" != "") save "`saving'", replace
	else save "`saving'"
	
	if ("`plot'" != "") {
		foreach var in `varlist' {
			graph twoway histogram `var' || kdensity `var', title("Conditional Distribution for Coefficient of `var'") 
			graph save `var', replace
		}
	}
	
	** Restore data **
	restore
end

