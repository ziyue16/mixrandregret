*! mixrpred 1.1.0 4Mar2022
*! [aut & dev] Ziyue Zhu  &  Álvaro A. Gutiérrez-Vargas

program define mixrpred, eclass
        version 12
	
        syntax newvarname [if] [in],  	///
		GRoup(varname) 			 		///
		ALTernatives(varname)			///
		ID(varname) 				    ///
		[PROBA							///
		NREP(integer 50)                ///
		BURN(integer 15)                ///
		XB]
        
		** Pre-checks **
		if "`e(cmd)'"=="" exit 301 
		
		if ("`proba'"=="" & "`xb'"=="") {
			loc proba = "proba"
			di "(option proba assumed; probability of success given one success within group)"
		}  
		
		** Mark the prediction sample **
		marksample touse , novarlist
		markout `touse' `group' `id' `e(basealternative)' `alternatives'
		
		** Generate variables used to sort data **
	    tempvar sorder altid
	    gen `sorder' = _n
	    sort `touse' `id' `group'
	    by `touse' `id' `group': gen `altid' = _n 
		
		** Drop data not in prediction sample **
	    preserve
	    qui keep if `touse'
		
		** Generate individual id **
	    if ("`id'" != "") {
		   tempvar nchoice pid
		   sort `group'
		   by `group': gen `nchoice' = cond(_n==_N,1,0)
		   sort `id'
		   by `id': egen `pid' = sum(`nchoice')		
		   qui duplicates report `id'
		   mata: mixrpre_np = st_numscalar("r(unique_value)")
		   mata: mixrpre_T = st_data(., st_local("pid"))
	     }
	     else {
		   qui duplicates report `group'
		   mata: mixrpre_np = st_numscalar("r(unique_value)")
		   mata: mixrpre_T = J(st_nobs(),1,1)
	}
        
		** Generate choice occacion id **
	    tempvar csid
	    sort `group'
	    by `group': egen `csid' = sum(1)
	    qui duplicates report `group'
	    local nobs = r(unique_value)

	    ** Sort data **
	    sort `id' `group' `altid'

		** Set Mata matrices to be used in prediction routine **
	    local rhs `e(indepvars)'
	    mata: mixrpre_X = st_data(., tokens(st_local("rhs")))
	    mata: mixrpre_CSID = st_data(., ("`csid'"))
	    local totobs = _N	 
	
     	** Restore data **
	    restore
		
		** panvar **
		mata: st_view(panvar = ., ., "`group'") 
		
		** X's **
		loc cmdline =  "`e(cmdline)'" 
		gettoken mixrandregret rhs: cmdline, parse(" ,")  // cmdline
		gettoken y x_options:  rhs, parse(" ,") 	      // dependent name
		gettoken covars options: x_options, parse(",")    // covars names
			
		** Mata allocation of the regressors **
		mata: st_view(X = ., ., "`covars'") 	
			
		** Parsing parameter vector **
		tempname b_all b_hat ASC_hat	
		matrix `b_all' = e(b)
		if "${cons_demanded}" =="NO"{ 
				matrix `b_hat' = `b_all' 	
		} 
		else if "${cons_demanded}" =="YES"{
				qui tab `alternatives' 
				local n_altern = r(r) 
					matrix `b_hat' 	 = `b_all'[1,1..(`e(rank)'-`n_altern')+1] 		 
					matrix `ASC_hat' = `b_all'[1,`e(rank)'-(`n_altern'-2)..`e(rank)']
		}
			
		** Generate ASC if needed **
		if "`e(ASC)'"=="YES" {
			// generate ASC as tempvars
			qui levelsof `alternatives', local(levels_altern)
			tempvar ASC_
			foreach i of local levels_altern {
				tempvar ASC_`i'
				qui gen `ASC_'`i' = (`alternatives' == `i')
			}				
			qui tab `alternatives' 
			local n_altern = r(r) 
			if "`e(basealternative)'"!=""{
				drop `ASC_'`e(basealternative)'
			}
			else{ // drop the alternative with lower number
				qui sum  `alternatives' , meanonly 
				local min_alt = r(min)
				qui drop `ASC_'`min_alt'
			}
				// Mata allocation of the ASC variables
			mata: st_view(ASC = ., ., "`ASC_'*") 		
			}
			else{ 
				// dummy ASC (=0 for conformability)
				mata: ASC = J(1,1,0) 	
				// coefficients ASC
				matrix `ASC_hat' =J(1,1,0)
			}				
		
		mata: b_hat= st_matrix("`b_hat'")		
		mata: b_all= st_matrix("`b_all'")
		mata: ASC_hat= st_matrix("`ASC_hat'")
		
		// Predicted Probability Computations
		mata: prediction= pbb_pred(X, ASC, panvar) 

		qui gen double	`varlist' = .	
		mata: empty_view = .
		mata: st_view(empty_view, ., "`varlist'")
		mata: empty_view[.,.] =prediction[.,.]
		qui replace `varlist' =. if e(sample)  != 1  
		qui replace `varlist' =. if   `touse' !=1 
}			
end
