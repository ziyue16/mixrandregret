
** old RRM function **
mata:
real matrix RRM_log_old(real matrix 	x_n, 
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
				r_i = ln(gamma :+ exp((b :/ mu) :* (x_n[j, . ] :- x_n[i, . ]))) 	
				r_i
				regret_n[i, . ] = regret_n[i, . ] :+ mu*r_i 		
				} 
			}	  
		}
	return(regret_n)
}
end

** exploration **
mata:
x1 = 23\27\35
x2 = 6\4\3
x_n = x1,x2
x_n

x_n[.,1]
x_n[1,.]
x_n[.,2]

d_12 = x_n[2,.] - x_n[1,.]
d_13 = x_n[3,.] - x_n[1,.]
d_13
d_12

/* expected structure */
d_ij = J(rows(x_n)*cols(x_n), rows(x_n), 0)
d_ij

b = 1, 2

/* test looping */
for(i=1; i <= rows(x_n); ++i) { 
		for(j = 1; j <= rows(x_n); ++j) {  
			d_ijm = x_n[j, . ] :- x_n[i, . ]
			if (j == 1) d_im = d_ijm'
			else d_im = d_im , d_ijm'
			}	
			if (i == 1) d_m = d_im
			else d_m = d_m \ d_im
		}
d_m

/* improved test looping */
d_ijm = J(rows(x_n)*cols(x_n), rows(x_n)-1, 0)
d_ijm

for(i = 1; i <= rows(x_n); ++i) { 
		for(j = 1; j <= rows(x_n); ++j) { 
		    if (i!=j){
			d_ijm = x_n[j, . ] :- x_n[i, . ]

			if ((j == 2 & i == 1) | (j == 1)) d_im = d_ijm
			else d_im = d_im \ d_ijm
			}
		}
		d_im
		regret_i = ln(1 :+ exp(b :* d_im))
		regret_i

		if (i == 1) regret_n = colsum(regret_i) 	
		else regret_n = regret_n \ colsum(regret_i)
}
regret_n

regret_old = RRM_log_old(x_n, b, 1, 1)
regret_old

/* get identical results */

end

** new RRM function **
mata:
real matrix RRM_log(real matrix 	x_n, 
					real rowvector 	b)	
{
	real scalar i, j
	
	for(i = 1; i <= rows(x_n); ++i) { 
		for(j = 1; j <= rows(x_n); ++j) { 
		    if (i!=j){
			d_ijm = x_n[j, . ] :- x_n[i, . ]

			if ((j == 2 & i == 1) | (j == 1)) d_im = d_ijm
			else d_im = d_im \ d_ijm
			}
		}
		d_im
		regret_i = ln(1 :+ exp(b :* d_im))
		regret_i

		if (i == 1) regret_n = colsum(regret_i) 	
		else regret_n = regret_n \ colsum(regret_i)
		}
	
	return(regret_n)
}
end
