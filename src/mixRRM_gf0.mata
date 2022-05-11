*! mixRRM_gf0 1.1.0 10May2022
*! [aut & dev] Ziyue Zhu & Álvaro A. Gutiérrez-Vargas


mata: mata set matastrict off 

/*=========================================================*/
/*=================== RRM FUNCTION ========================*/
/*=========================================================*/

mata:
real matrix RRM_log(real matrix    x_n,
                    real rowvector b,
                    real scalar    mu,
                    real scalar    gamma)
{
	real scalar i, j, r_m
	real matrix regret_n
	
	regret_n = J(rows(x_n), cols(x_n), 0)
	for(i = 1; i <= rows(x_n); ++i) {
		for(j = 1; j <= rows(x_n); ++j) {
			if (i != j) {
				r_i = ln(gamma :+ exp((b :/ mu) :* (x_n[j, .] :- x_n[i, .])))
				regret_n[i, . ] = regret_n[i, . ] :+ mu*r_i
			}
		}
	}
	return(regret_n)
}
end

/*=============================================================*/
/*=================== MIXED RRM FUNCTION ======================*/
/*=============================================================*/

version 12
mata:
void mixRRM_gf0(transmorphic scalar MM,
                real scalar todo, 
                real rowvector b, 
                real colvector lnfj, /* <- a colvector with all the individual loglikelihood constributions */
                S, 
                H)
 {
 
 	/*--------------------------*/
	/*--- External Variables ---*/
	/*--------------------------*/
	
	external mixr_T     // observations in the (long) data  [N x J x T]
	external mixr_CSID  // choice sets ID
	external mixr_IND   // individual ID
	external mixr_nrep  
	external mixr_np
	external mixr_kfix
	external mixr_krnd
	external mixr_krln
	external mixr_burn
	external mixr_corr
	external mixr_user

	nrep = mixr_nrep
	np   = mixr_np
	kfix = mixr_kfix
	krnd = mixr_krnd
	krln = mixr_krln
	burn = mixr_burn
	corr = mixr_corr
	user = mixr_user
	ID_IND = mixr_IND 

	/*-----------------------------*/
	/*--- Variables Declaration ---*/
	/*-----------------------------*/
	
	real matrix subject     // general info about individuals 
	real matrix N_subject   // intermediate step to extract info about num panels 

	/* Loop scalar variables */
	real scalar n           // for individual
	real scalar t           // for choice set
	real scalar i           // alternative i 
	real scalar j           // alternative j 
	real scalar m           // index 
	real scalar n_rows      // scalar for num rows
	real scalar n_cols      // scalar for num columns
	
	/* Matrix explanatory and explained variable */
	real matrix Y           // explained variable
	real matrix X           // covariates
	
	/* Sub matrices per individual */
	real matrix x_n         // sub matrix covariates individual n 
	real matrix y_n         // sub matrix   choice   individual n 
	
	cons_demanded = st_global("cons_demanded")

	/*--------------------------*/
	/*--- Variables Creation ---*/
	/*--------------------------*/
	
	mixr_Y = moptimize_util_depvar(MM, 1)       // recovering the explained variable
	mixr_X = moptimize_init_eq_indepvars(MM, 1) // recovering the explanatory variable
	
	if (cons_demanded=="YES"){
		ASC = moptimize_init_eq_indepvars(MM, 3) // 3rd equation estimates ASC
		id_ASC_eq = moptimize_util_eq_indices(MM, 3)
		b_ASC = b[|id_ASC_eq|]
	}  

	/*------------------------------------------------*/
	/*--- Parse Parameter vector (Random vs Fixed) ---*/
	/*------------------------------------------------*/
	
	B = b' // the program works with parameters (b) in column vector format (B = b') 

	kall = kfix + krnd // total Number of Parameters (kall) = fixed (kfix) + random (krand)
	if (kfix > 0) {
		MFIX = B[|1, 1\kfix, 1|]
		MFIX = MFIX :* J(kfix, nrep, 1)
	}
	
	/* Recover mean and SD */
	MRND = B[|(kfix+1), 1\kall, 1|] // recovering the mean of the random parameters (MRND, "M" for \mu)
	SRND = diag(B[|(kall+1), 1\(kfix+2*krnd), 1|]) // independet distributions only take the diagonal
	
	/*------------------------------*/
	/*--- Compute log-likelihood ---*/
	/*------------------------------*/
	
	/* Set up panel information */
	subject = panelsetup(ID_IND, 1) // input is [N x 1] vector; output is [N_subject x 2] matrix
	N_subject = panelstats(subject)[1] // # of panel units (subjects), identified by number of rows in "panel"; given this 1 we recover the Num. of panels
	
	st_view(panvar = ., ., st_global("group_mata")) // panel information for ASC
	subject_ASC = panelsetup(panvar, 1)
	npanels = panelstats(subject_ASC)[1]
	
	moptimize_init_by(MM, ID_IND) // log-likelihood is computed at subject level
	
	if (user == 1) external mixr_USERDRAWS
	
	/* Object P */
	P = J(np, 1, 0) // will store individual contributions to the LL of each individual
	
	m = 1
	/* Loop over individuals (npanels) */
	for(n = 1; n <= N_subject; n++) {
		
		if (user == 1) {
			ERR = invnormal(mixr_USERDRAWS[|1, (1+nrep*(n-1))\krnd, (nrep*n)|])
		}
		else {/* Regular (non-scrambled) Halton integration */
			ERR = invnormal(halton(nrep, krnd, (1+burn+nrep*(n-1)))')
		}
		
		/* Modify parameters to construct the [beta = mu + sigma * draw] structure. */
		if (kfix > 0) BETA = MFIX \ (MRND :+ (SRND*ERR)) 
		else BETA = MRND :+ (SRND*ERR) // get all random variables
		 
		/* Log-normal distribution */
		if (krln > 0) {
			if ((kall-krln) > 0) {
				BETA = BETA[|1, 1\(kall-krln), nrep|]\exp(BETA[|(kall-krln+1), 1\kall, nrep|])
			}
			else {
				BETA = exp(BETA)
			}
		}

		/* Object R */
		R = J(1, nrep, 1) // will store the product of sequences of choices made by the same individual
		
		/* Looping over choice situations (t) of individual (n) */
		t = 1 
		nc = mixr_T[m, 1] // number of choice sets
		
		for(t = 1; t <= nc; t++) {

			/* Parse the matrix X and Y at the level of choice situation (t) */
			YMAT = mixr_Y[|m, 1\(m+mixr_CSID[m,1]-1), cols(mixr_Y)|]
			XMAT = mixr_X[|m, 1\(m+mixr_CSID[m,1]-1), cols(mixr_X)|]
			
			/* Shape of block individual n */
			n_rows = rows(XMAT) // number of faced alternatives
			n_cols = cols(XMAT) // number of covariates 
			
			ER_rep = J(n_rows, nrep, 0) // store the regret of each draws
			
			if (cons_demanded == "YES") {
				asc = panelsubmatrix(ASC, 1, subject_ASC)
				ASC_prod = asc*b_ASC'
			}
			
			for(rr = 1; rr <= nrep; rr++) {
				regret_draw_r = RRM_log(XMAT, BETA[.,rr]', 1, 1) // regret
				
				if (cons_demanded == "YES") { 
					ER_rep[ ., rr] =rowsum(regret_draw_r :+ ASC_prod) // with ASC condition
				}
				else if (cons_demanded == "NO") {
					ER_rep[ ., rr] =rowsum(regret_draw_r) 
				}
			}
           
			ER = exp(-ER_rep) // the negative of the regret
			ER = (ER :/ colsum(ER, 1)) 
			
			R = R :* colsum(YMAT :* ER, 1) // product of repeated choices from the same individual
			
			m = m + mixr_CSID[m, 1] // move the index for the next choice situation
		}

		P[n, 1] = mean(R', 1) // average all the probabilities over all the draws
	}
		
	/* The vector containing the LL of every individual is simply the log of the probability. */
	lnfj = ln(P) 
}
end



/*=====================================================================*/
/*================== MIXED RRM PREDICT FUNCTION =======================*/
/*=====================================================================*/

mata:
real matrix pbb_pred(real matrix X, real matrix ASC, real colvector panvar)
{

	/*--------------------------*/
	/*--- External Variables ---*/
	/*--------------------------*/
	
	external mixrpre_X
	external mixrpre_T
	external mixrpre_CSID
	external mixrpre_np
	external b_hat
	external b_all
	external ASC_hat 
	external ASC
	external mixr_IND
	
	np = mixrpre_np
	ID_IND = mixr_IND 
	command = st_local("cmd")
	nrep = strtoreal(st_local("nrep"))
	kfix = st_numscalar("e(kfix)")
	krnd = st_numscalar("e(krnd)")
	krln = st_numscalar("e(krln)")
	burn = strtoreal(st_local("burn"))
	user = st_numscalar("e(userdraws)")
	cons_demanded = st_global("cons_demanded")
	proba = st_local("proba") 

	/*------------------------------*/
	/*--- Parse Parameter vector ---*/
	/*------------------------------*/
	
	B_hat = b_all' 

	kall = kfix + krnd
	
	if (kfix > 0) {
		MFIX_hat = B_hat[|1, 1\kfix, 1|]
		MFIX_hat = MFIX_hat :* J(kfix, nrep, 1)
	}
	
	MRND_hat = B_hat[|(kfix+1), 1\kall, 1|] 
	SRND_hat = diag(B_hat[|(kall+1), 1\(kfix+2*krnd), 1|]) 
	
	/*--------------------------*/
	/*--- Compute Prediction ---*/
	/*--------------------------*/
	
	/* Set up panel information */
	subject = panelsetup(ID_IND, 1)
	N_subject = panelstats(subject)[1] 

	/* Panel information for ASC */
	subject_ASC = panelsetup(panvar, 1)
	npanels = panelstats(subject_ASC)[1]
	
	if (user == 1) external mixr_USERDRAWS

	m = 1
	for(n = 1; n <= N_subject; n++) {
	  
		if (user == 1) {
			ERR = invnormal(mixr_USERDRAWS[|1, (1+nrep*(n-1))\krnd, (nrep*n)|])
		}
		else {
			ERR = invnormal(halton(nrep, krnd, (1+burn+nrep*(n-1)))')
		}
		
		if (kfix > 0) BETA_hat = MFIX_hat \ (MRND_hat :+ (SRND_hat*ERR))
		else BETA_hat = MRND_hat :+ (SRND_hat*ERR)
		
		if (krln > 0) {
			if ((kall-krln) > 0) {
				BETA_hat = BETA_hat[|1, 1\(kall-krln), nrep|] \ exp(BETA_hat[|(kall-krln+1), 1\kall, nrep|])
			}
			else {
				BETA_hat = exp(BETA_hat)
			}
		}
		
		t = 1 
		nc = mixrpre_T[m, 1]
		
		for(t = 1; t <= nc; t++) {

			XMAT = mixrpre_X[|m, 1\(m+mixrpre_CSID[m, 1]-1), cols(mixrpre_X)|]

			n_rows = rows(XMAT)
			n_cols = cols(XMAT)
			
			if (cons_demanded == "YES") {
				asc = panelsubmatrix(ASC, 1, subject_ASC)
				ASC_prod = asc*ASC_hat'
			}
			
			ER_rep_hat = J(n_rows, nrep, 0)
			for(rr = 1; rr <= nrep; rr++) {
				regret_draw_hat = RRM_log(XMAT, BETA_hat[.,rr]', 1, 1)

				if (cons_demanded == "YES") {
					ER_rep_hat[., rr] = rowsum(regret_draw_hat :+ ASC_prod)
				}
				else if (cons_demanded == "NO") {
					ER_rep_hat[., rr] =rowsum(regret_draw_hat)
				}
			}
			
			ER_hat = exp(-ER_rep_hat)
			p_hat_i = (ER_hat :/ colsum(ER_hat, 1))
			
			N_cs = mixrpre_CSID[1]
			/* Object output */
			p_out_cs = J(N_cs, 1, 0)
			regret_out_cs = J(N_cs, 1, 0)
			
			for(aa = 1; aa <= N_cs; aa++) {
				p_out_cs[aa, 1] = mean(p_hat_i[aa,.]', 1)
				regret_out_cs[aa, 1] = mean(ER_rep_hat[aa,.]', 1)
			}
			
			m = m + mixrpre_CSID[m, 1]

			/* Probability prediction: proba option [default] */
			if (proba != "") {
				if (t == 1) p_hat_cs = p_out_cs
				else	    p_hat_cs = p_hat_cs \ p_out_cs
			}
			else {
				/* Linear prediction: xb option */
				if (t == 1) regret_hat_cs = regret_out_cs
				else	    regret_hat_cs = regret_hat_cs \ regret_out_cs
			}
		}
		/* Append individuals */
		if (proba != "") {
			if (n == 1) p_hat = p_hat_cs
			else	    p_hat = p_hat \ p_hat_cs
		}
		else {
			if (n == 1) regret_hat = regret_hat_cs
			else	    regret_hat = regret_hat \ regret_hat_cs
		}
	}
	/* Output prediction */
	if (proba != "") return(p_hat)
	else return(regret_hat)
}
end



/*=====================================================================*/	
/*================== INDIVIDUAL LEVEL PARAMETERS ======================*/	
/*=====================================================================*/

version 9.2
mata:
function mixr_beta(string scalar B_s)
{

	/*--------------------------*/
	/*--- External Variables ---*/
	/*--------------------------*/
	
	external mixrbeta_X
	external mixrbeta_Y
	external mixrbeta_T
	external mixrbeta_CSID
	external mixrbeta_PB
	external mixrbeta_np
	external mixrbeta_ID
	external b_beta
	external b_all
	external ASC_beta
	external ASC

	np = mixrbeta_np
	ID_IND = mixrbeta_ID 
	command = st_local("cmd")
	nrep = strtoreal(st_local("nrep"))
	kfix = st_numscalar("e(kfix)")
	krnd = st_numscalar("e(krnd)")
	krln = st_numscalar("e(krln)")
	burn = strtoreal(st_local("burn"))
	user = st_numscalar("e(userdraws)")
	cons_demanded = st_global("cons_demanded")

	/*------------------------------*/
	/*--- Parse Parameter vector ---*/
	/*------------------------------*/
	
	B_beta = b_all'

	kall = kfix + krnd
	
	if (kfix > 0) {
		MFIX_beta = B_beta[|1, 1\kfix, 1|]
		MFIX_beta = MFIX_beta :* J(kfix, nrep, 1)
	}
	
	MRND_beta = B_beta[|(kfix+1), 1\kall, 1|]
	SRND_beta = diag(B_beta[|(kall+1), 1\(kfix+2*krnd), 1|])
	
	/*--------------------------*/
	/*--- Compute Prediction ---*/
	/*--------------------------*/
	
	/* Set up panel information */
	subject = panelsetup(ID_IND, 1)
	N_subject = panelstats(subject)[1]

	/* Panel information for ASC */
	st_view(panvar = ., ., st_global("group_mata")) // panel information for ASC
	subject_ASC = panelsetup(panvar, 1)
	npanels = panelstats(subject_ASC)[1]
			
	if (user == 1) external mixr_USERDRAWS
	
	/* Object P */
	P = J(np, 1, 0)
	
	m = 1
	/* Loop over individuals (npanels) */
	for(n = 1; n <= N_subject; n++) {
		
		if (user == 1) {
			ERR = invnormal(mixr_USERDRAWS[|1, (1+nrep*(n-1))\krnd, (nrep*n)|])
		}
		else {
			/* Regular (non-scrambled) Halton integration */
			ERR = invnormal(halton(nrep, krnd, (1+burn+nrep*(n-1)))')
		}
		
		
		/* Modify parameters to construct the [beta = mu + sigma * draw] structure. */
		if (kfix > 0) BETA_new = MFIX_beta \ (MRND_beta :+ (SRND_beta*ERR)) 
		else BETA_new = MRND_beta :+ (SRND_beta*ERR) // get all random variables
		
		/* Log-normal distribution */
		if (krln > 0) {
			if ((kall-krln) > 0) {
				BETA_new = BETA_new[|1, 1\(kall-krln), nrep|] \ exp(BETA_new[|(kall-krln+1), 1\kall, nrep|])
			}
			else {
				BETA_new = exp(BETA_new)
			}
		}

		/* Object R */
		R = J(1, nrep, 1)
		
		/* Looping over choice situations (t) of individual (n) */
		t = 1
		nc = mixrbeta_T[m, 1] // number of choice sets
		
		for(t = 1; t <= nc; t++) {
		
			/* Parse the matrix X and Y at the level of choice situation (t) */
			YMAT = mixrbeta_Y[|m, 1\(m+mixrbeta_CSID[m,1]-1), cols(mixrbeta_Y)|]
			XMAT = mixrbeta_X[|m, 1\(m+mixrbeta_CSID[m,1]-1), cols(mixrbeta_X)|]
			
			/* Shape of block individual n */
			n_rows=rows(XMAT) 
			n_cols=cols(XMAT)
			
			if (cons_demanded == "YES") {
				asc = panelsubmatrix(ASC, 1, subject_ASC)
				ASC_prod = asc*ASC_beta'
			}
			
			ER_rep_beta = J(n_rows, nrep, 0) 
			for(rr = 1; rr <= nrep; rr++) {
				regret_draw_beta = RRM_log(XMAT, BETA_new[., rr]', 1, 1)

				if (cons_demanded == "YES") {
					ER_rep_beta[., rr] =rowsum(regret_draw_beta :+ ASC_prod)
				}
				else if (cons_demanded == "NO") {
					ER_rep_beta[ .,rr] =rowsum(regret_draw_beta)
				}
			}
			
			ER_beta = exp(-ER_rep_beta)
			ER_beta = (ER_beta :/ colsum(ER_beta, 1))
			
			R = R :* colsum(YMAT :* ER_beta, 1)
			m = m + mixrbeta_CSID[m, 1]
		}
		P[n, 1] = mean(R', 1)
		mixrbeta_PB[n, .] = mean((R :* BETA_new)', 1) / P[n, 1]
	}
}
end

exit

