*! mixRRM_gf0 0.0.1 12Feb2022
*! authors Ziyue Zhu  &  Álvaro A. Gutiérrez-Vargas


mata: mata set matastrict off 

/*--------------------*/
/*--- RRM Function ---*/
/*--------------------*/

mata:
real matrix RRM_log(real matrix 	x_n, 
					real rowvector 	b, 
					real scalar 	mu, 
					real scalar 	gamma)	
{
	real scalar i, j ,r_m
	real matrix regret_n 
	
	regret_n = J(rows(x_n), cols(x_n), 0)
	for(i=1; i <= rows(x_n); ++i) { 
		for(j=1; j <= rows(x_n); ++j) { 
			if (i!=j) { 
				r_i = ln(gamma :+ exp( (b :/ mu) :* ( x_n[j , . ] :-  x_n[i, . ]))) 				
				regret_n[i, . ] = regret_n[i, . ] :+ mu*r_i 		
				} 
			}	  
		}
	return(regret_n)
}
end


version 12
mata:
void mixRRM_gf0(transmorphic scalar MM,
                real scalar todo, 
                real rowvector b, 
                real colvector  lnfj, /*<- lnfj is a colvector with all the individual loglikelihood constributions*/
                S, 
                H)
 {
 
 	/*---------------------------*/
	/*--- external variables  ---*/
	/*---------------------------*/
	
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
	ID_IND = mixr_IND  // Variable that stores Individuals 


	/*-----------------------------*/
	/*--- variables declaration ---*/
	/*-----------------------------*/
	
	real matrix subject     // general info information about individuals 
	real matrix N_subject   // intermediate step to extract into about num panels 

	/*	Loop scalar variables */
	real scalar n // for individual
	real scalar t // for choice set
	real scalar i // alternative i 
	real scalar j // alternative j 
	real scalar m // alternative j 
	real scalar l // for nrep
	real scalar n_rows // scalar for num rows
	real scalar n_cols // scalar for num columns
	
	/*	Matrix explanatory and explained variable */
	real matrix Y // explained variable
	real matrix X // covariates
	
	/*	Sub matrices pero individual */
	real matrix r_i // r_i matrix with regret 
	real matrix x_n // sub matrix covariates individual n 
	real matrix y_n // sub matrix   choice   individual n 
	
	//cons_demanded = st_global("cons_demanded")

	/*--------------------------*/
	/*--- variables creation ---*/
	/*--------------------------*/
	
	mixr_Y = moptimize_util_depvar(MM, 1) 	// recovering the explained variable
	mixr_X = moptimize_init_eq_indepvars(MM, 1)	// recovering the explanatory variable
	
	if (cons_demanded=="YES"){
		ASC = moptimize_init_eq_indepvars(MM,2) 
		id_ASC_eq = moptimize_util_eq_indices(MM,2)
		b_ASC = b[|id_ASC_eq|]	
	}  

	/*------------------------------------------------*/
	/*--- Parse Parameter vector (Random vs Fixed) ---*/
	/*------------------------------------------------*/
	
	/* the program works with parameters (b) in column vector format (B = b') */
	B = b' 
	
	/* Total Number of Parameters (kall) = fixed (kfix) + random (krand) */
	kall = kfix + krnd
	
	if (kfix > 0) {
		MFIX = B[|1,1\kfix,1|]
		MFIX = MFIX :* J(kfix,nrep,1)	
	}

	/* Recovering the mean of the random parameters parameters (MRAND, "M" for \mu ) */
	MRND = B[|(kfix+1),1\kall,1|] 
	
	/* if the user only want independet distributions then only take the diagonal */
    SRND = diag(B[|(kall+1),1\(kfix+2*krnd),1|])   // no correlation for now
	
	/*------------------------------*/
	/*--- Compute log-likelihood ---*/
	/*------------------------------*/
	
	/* set up panel information: input is [N x 1] vector "id", and output is [N_subject x 2] matrix "subject" */
	subject = panelsetup(ID_IND,1)
	
	/* # of panel units (subjects), identified by number of rows in "panel" */
	N_subject = panelstats(subject)[1] // given this 1 we recover the Num. of panels
		
	/* tell Stata that log-likelihood is computed at subject level */
	moptimize_init_by(MM,ID_IND)	
	
	if (user == 1) external mixr_USERDRAWS
	
	/* Object P: will store individual contributions to the LL of each individual */
	P = J(np,1,0)
	
	m = 1                           
	/* loop over individuals (npanels) */
	for(n=1; n <= N_subject; n++) {
		
		if (user == 1) {
			ERR = invnormal(mixr_USERDRAWS[|1,(1+nrep*(n-1))\krnd,(nrep*n)|])
		}
		else {
            /* Regular (non-scrambled) Halton integration */
			ERR = invnormal(halton(nrep,krnd,(1+burn+nrep*(n-1)))')
		}
		
		/* Modify parameters to construct the [beta = mu + sigma * draw] structure. */
		if (kfix > 0) BETA = MFIX \ (MRND :+ (SRND*ERR)) 
		else BETA = MRND :+ (SRND*ERR) // since we get all random variables
		 
		/* lognormal distribution */
		if (krln > 0){
            if ((kall-krln) > 0) {
                BETA = BETA[|1,1\(kall-krln),nrep|] \ exp(BETA[|(kall-krln+1),1\kall,nrep|])
                }
                else {
                BETA = exp(BETA)
                }
		}

		/* Object R: will store the product of sequences of choices made by the same individual */
		R = J(1,nrep,1)		
		
		/* Looping over choice situations (t) of individual (n) */
		t = 1 
		nc = mixr_T[m,1] // number of choice sets
        // loop over choice sets (t)
        for(t=1; t<=nc; t++) {
		
		    /* Parse the matrix X and Y at the level of choice situation (t) */
			YMAT = mixr_Y[|m,1\(m+mixr_CSID[m,1]-1),cols(mixr_Y)|]
			XMAT = mixr_X[|m,1\(m+mixr_CSID[m,1]-1),cols(mixr_X)|]
			
			/* shape of block individual n  */	
		    n_rows=rows(XMAT) // Number of faced alternatives
		    n_cols=cols(XMAT) // Number of covariates 
		
			/* Start computing the multinomial probability */
			/* This is the part where he created the probability of equation L_[nit] in mixed random regret */		   
			
			/* ER_rep matrix to store the regret of each draw. */
			ER_rep = J(n_rows, nrep, 0)
			for(rr=1; rr<=nrep; rr++) {
					// Create regret for each alternative 
					regret_draw_r = RRM_log(XMAT,BETA[.,rr]',1,1) 
					// Create regret for each alternative and generate the row sum			
				   // ------  ASC  ------- //		
		           if (cons_demanded=="YES") { 
			           asc_n 	= panelsubmatrix(ASC, n, paninfo) 
			           ASC_prod= asc_n*b_ASC'
			           ER_rep[ .,rr] =rowsum(regret_draw_r) + ASC_prod 
				       }
		           else if (cons_demanded=="NO"){
			           ER_rep[ .,rr] =rowsum(regret_draw_r) 
					 }
					}
	
			ER = exp(-ER_rep) // The negative of the regret
			ER = (ER :/ colsum(ER,1)) 

			/* This generates the product of repeated choices from the same individual */
			R = R :* colsum(YMAT :* ER, 1) 
			 
			/* move the index for the next choice situation */
			 m = m + mixr_CSID[m,1]
			}
		/* Here we average all the probabilities over all the draws */
		P[n, 1] = mean(R',1)
		}
 		
	 /* The vector containing the LL of every individual is simply the log of the probability. */
	 lnfj = ln(P) 
}
end

