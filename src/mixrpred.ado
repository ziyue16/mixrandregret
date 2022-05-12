*! mixrpred 1.1.0 8Apr2022
*! [aut & dev] Ziyue Zhu & Álvaro A. Gutiérrez-Vargas

program define mixrpred, eclass
		version 9.2
	
		syntax newvarname [if] [in], ///
		[PROBA                       ///
		NREP(integer 50)             ///
		BURN(integer 15)             ///
		XB]

		** Pre-checks **
		if "`e(cmd)'" == "" exit 301
		
		if ("`proba'" == "" & "`xb'" == "") {
			loc proba = "proba"
			di "(option proba assumed; probability of success given one success within group)"
		}
		
		** Mark the prediction sample **
		marksample touse, novarlist
		markout `touse' `e(indepvars)' `e(group)' `e(id)' `e(alternatives)'
		
		** Generate variables used to sort data **
		tempvar sorder altid
		gen `sorder' = _n
		sort `touse' `e(id)' `e(group)'
		by `touse' `e(id)' `e(group)': gen `altid' = _n 
		
		** Drop data not in prediction sample **
		preserve
		qui keep if `touse'
		
		** Generate individual id **
		if ("`e(id)'" != "") {
			tempvar nchoice pid
			sort `e(group)'
			by `e(group)': gen `nchoice' = cond(_n==_N, 1, 0)
			sort `e(id)'
			by `e(id)': egen `pid' = sum(`nchoice')
			qui duplicates report `e(id)'
			mata: mixrpre_np = st_numscalar("r(unique_value)")
			mata: mixrpre_T = st_data(., st_local("pid"))
		}
		else {
			qui duplicates report `e(group)'
			mata: mixrpre_np = st_numscalar("r(unique_value)")
			mata: mixrpre_T = J(st_nobs(), 1, 1)
		}

		** Generate choice occacion id **
		tempvar csid
		sort `e(group)'
		by `e(group)': egen `csid' = sum(1)
		qui duplicates report `e(group)'

		** Sort data **
		sort `e(id)' `e(group)' `altid'

		** Set Mata matrices to be used in prediction routine **
		local rhs `e(indepvars)'
		mata: mixrpre_X = st_data(., tokens(st_local("rhs")))
		mata: mixrpre_CSID = st_data(., ("`csid'"))
	
		** Restore data **
		restore
		
		** panvar **
		mata: st_view(panvar = ., ., "`e(group)'")
		
		** X's **
		loc cmd = "`e(cmd)'" 
		gettoken mixrandregret rhs: cmdline, parse(",") // cmdline
		gettoken y x_options: rhs, parse(" ,")          // dependent name
		gettoken covars options: x_options, parse(",")  // covars names

		** Mata allocation of the regressors **
		mata: st_view(X = ., ., "`covars'")
			
		** Parsing parameter vector **
		tempname b_all b_hat ASC_hat
		matrix `b_all' = e(b)
		if "${cons_demanded}" == "NO" {
			matrix `b_hat' = `b_all'
		}
		else if "${cons_demanded}" == "YES" {
			qui tab `e(alternatives)'
			local n_altern = r(r) 
			matrix `b_hat' = `b_all'[1, 1..(`e(rank)'-`n_altern')+1]
			matrix `ASC_hat' = `b_all'[1, `e(rank)'-(`n_altern'-2)..`e(rank)']
		}
			
		** Generate ASC if needed **
		if "`e(ASC)'" == "YES" {
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
				qui sum `e(alternatives)', meanonly
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
			matrix `ASC_hat' =J(1, 1, 0)
		}
		
		mata: b_hat = st_matrix("`b_hat'")
		mata: b_all = st_matrix("`b_all'")
		mata: ASC_hat = st_matrix("`ASC_hat'")

		** Predicted Probability Computations **
		mata: prediction = pbb_pred(X, ASC, panvar)
		
		mata: st_store(., st_addvar("float", "`varlist'"), prediction)
end
