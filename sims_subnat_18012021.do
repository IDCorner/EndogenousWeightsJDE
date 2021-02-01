clear all
set maxvar 120000 
set matsize 11000 
global path "/Users/ricardonogales/Dropbox/Endogenous Weights/Datasets and dofiles/"

* Dataset Creation

local cutoff 50 // Select the cutoff for which one wants to create the datset
gen a=""
save "$path/temp1_k`cutoff'",replace
label def method 1 "Exogenous" 2 "Frequency" 3 "PCA"
label def region 0 "National" 1 "Region 1" 2 "Region 2"
foreach x of numlist 1/10 {
	set seed 12`x'
	foreach j of numlist 1/6 {
		local c`j'=runiform(0.2,0.8) 
	}
	local m=runiform(0.4,0.5)
	local m1=0.5*`m'
	local m2=`m'
	local m3=1.5*`m'
	local m4=2*`m'
	matrix C = (1, `c1', `c2', `c3' \ ///
				`c1', 1, `c4', `c5' \ ///
				`c2', `c4', 1, `c6' \ ///
				`c3', `c5', `c6', 1)	

	*corr2data i1 i2 i3 i4, n(500) corr(C) mean(`m1' `m2' `m3' `m4')
	drawnorm i1 i2 i3 i4, n(500) corr(C) mean(`m1' `m2' `m3' `m4')
	foreach i of numlist 1/4 {
		replace i`i'=normal(i`i')
		replace i`i'=(i`i'>=0.8)
	}
	gen region=1
	save "$path/aaa.dta", replace

	* Dataset Creation Region 2
	clear 
	set seed 13`x'
	foreach j of numlist 1/6 {
		local c`j'=runiform(0.2,0.3)
	}
	local m=runiform(0.4,0.5)
	local m1=0.5*`m'
	local m2=`m'
	local m3=1.5*`m'
	local m4=2*`m'
	matrix C = (1, `c1', `c2', `c3' \ ///
				`c1', 1, `c4', `c5' \ ///
				`c2', `c4', 1, `c6' \ ///
				`c3', `c5', `c6', 1)	

	*corr2data i1 i2 i3 i4, n(500) corr(C) mean(`m1' `m2' `m3' `m4')
	drawnorm i1 i2 i3 i4, n(500) corr(C) mean(`m1' `m2' `m3' `m4')
	foreach i of numlist 1/4 {
		replace i`i'=normal(i`i')
		replace i`i'=(i`i'>=0.5)
	}
	gen region=2
	append using "$path/aaa.dta"
	erase "$path/aaa.dta"

	mean i*, over(region)
	bys region: tetrachoric i*
	tetrachoric i*

	egen count=rowtotal(i*)
	bys region: tab count

	** Baseline poverty measures
		* Exogenous weights
		foreach i of numlist 1/4 {
			gen exo_w`i'=1/4
		}

		* Freq weights
		mean i* 
		mat E=e(b)
		local sum=0
		foreach i of numlist 1/4 {
			local sum=`sum'-ln(E[1,`i'])
		}
		foreach i of numlist 1/4 {
			gen freq_w`i'=-ln(E[1,`i'])/`sum'
		}

		* PCA weights
		tetrachoric i*, star(0.05) 
		mat C=r(corr)
		local n=_N
		pcamat C, n(`n') factor(1) forcepsd // we retain only the first principal component
		estat loadings
		mat S1=r(A) // Original loadings = non-normalized weights
		local sum=0
		foreach i of numlist 1/4 {
			local sum=`sum'+S1[`i',1] // Sum of non-normalized weights
		}
		foreach i of numlist 1/4 {
			gen pca_w`i'=S1[`i',1]/`sum' 
		}

		* c scores
		foreach w in exo freq pca {
			foreach i of numlist 1/4 {
				gen g0_`i'_`w'=i`i'*`w'_w`i'
			}
			egen c_`w'=rowtotal(g0*`w')
		}

		* individual poverty status
		foreach k of local cutoff {
			foreach w in exo freq pca {
				gen poor_`w'_`k'=(c_`w'>=`k'/100)
			}

			* intensity among the poor
			foreach w in exo freq pca {
				gen cens_c_`w'_`k'=c_`w'
				replace cens_c_`w'_`k'=0 if poor_`w'_`k'==0
			}
		}
		
		* Initialize matrix to store results
		mat R=J(1,7,.)
		mat colnames R=sim_port sim_ind k w_method mpi subnat rep

		* adjusted poverty headcount ratio = the societal poverty measure of interest
		foreach k of local cutoff {
			local j=0
			foreach w in exo freq pca {
				local j=`j'+1
				
				* National
				sum cens_c_`w'_`k'
				mat B=(0,.,`k',`j',r(mean),0,0) // Last zero is for repetition
				mat R=R\B
				
				* Region 1
				sum cens_c_`w'_`k' if region==1 // Last zero is for repetition
				mat B=(0,.,`k',`j',r(mean),1,0)
				mat R=R\B
				
				* Region 2
				sum cens_c_`w'_`k' if region==2 // Last zero is for repetition
				mat A=e(b)
				mat B=(0,.,`k',`j',r(mean),2,0)
				mat R=R\B
			}
		}
		
		
	** Simulated poverty measures, with synthetic deprivations affecting different proportions of the population
		local sreg 1 // Simulated region
		local rep 10 // number of repetitions for simulation
		local sim 1 2 3 4 // simulated list of indicators
		
		* Clone variables
		foreach i of numlist 1/4 {
			foreach s of local sim { // simulated indicator
				foreach r of numlist 1/`=`rep'' { // number of replications in simulation
					foreach p of numlist 1/19 { // simulated proportion of additional indicators
						clonevar i`i'_s`s'_p`p'_r`r'=i`i' // variables are identical at first
					}
				}
			}
		}

	* Creation of added deprivations
		foreach s of local sim { // simulated indicator
			foreach r of numlist 1/`=`rep'' { // number of replications in simulation
				set seed `r'
				gen u_`r'=runiform() if i`s'==0  & region==`sreg' // random indentifier
				xtile pct_`r' = u_`r', nq(20) // steps of 5%
				foreach p of numlist 1/19 {
					replace i`s'_s`s'_p`p'_r`r'=1 if pct_`r'<=`p' // all deprivations concentrated in region 1
				}
				drop u_`r' pct_`r'
			}
		}
		
	* Recalculation of MPIs
		* Exogenous weights // Do not change
		foreach i of numlist 1/4 {
			foreach s of local sim { // simulated indicator
				foreach r of numlist 1/`=`rep'' { // number of replications in simulation
					foreach p of numlist 1/19 { // simulated proportion of additional indicators
						clonevar exo_w`i'_s`s'_p`p'_r`r'=exo_w`i' // variables are identical at first
					}
				}
			}
		}
		
		
		* Freq weights
		foreach s of local sim { // simulated indicator
			foreach p of numlist 1/19 {
				foreach r of numlist 1/`=`rep'' {
					mean i*s`s'_p`p'_r`r'
					mat E=e(b)
					local sum=0
					foreach i of numlist 1/4 {
						local sum=`sum'-ln(E[1,`i'])
					}
					foreach i of numlist 1/4 {
						gen freq_w`i'_s`s'_p`p'_r`r'=-ln(E[1,`i'])/`sum'
					}
				}
			}
		}
		
		* PCA weights
		foreach s of local sim { // simulated indicator
			foreach p of numlist 1/19 {
				foreach r of numlist 1/`=`rep'' {	
					tetrachoric i*_s`s'_p`p'_r`r'
					mat C=r(corr)
					local n=_N
					pcamat C, n(`n') factor(1) forcepsd // we retain only the first principal component
					estat loadings
					mat S1=r(A) // Original loadings = non-normalized weights
					local sum=0
					foreach i of numlist 1/4 {
						local sum=`sum'+S1[`i',1] // Sum of non-normalized weights
					}
					foreach i of numlist 1/4 {
						gen pca_w`i'_s`s'_p`p'_r`r'=S1[`i',1]/`sum' 
					}
				}
			}
		}
		
		* c scores
		foreach w in exo freq pca {
			foreach s of local sim { // simulated indicator
				foreach p of numlist 1/19 {
					foreach r of numlist 1/`=`rep'' {
						foreach i of numlist 1/4 {
							gen g0_`i'_`w'_s`s'_p`p'_r`r'=i`i'_s`s'_p`p'_r`r'*`w'_w`i'_s`s'_p`p'_r`r'
						}
						egen c_`w'_s`s'_p`p'_r`r'=rowtotal(g0*`w'_s`s'_p`p'_r`r')
					}
				}
			}
		}

		* individual poverty status
		foreach k of local cutoff {
			foreach w in exo freq pca {
				foreach s of local sim { // simulated indicator
					foreach p of numlist 1/19 {
						foreach r of numlist 1/`=`rep'' {
							gen poor_`w'_`k'_s`s'_p`p'_r`r'=(c_`w'_s`s'_p`p'_r`r'>=`k'/100)
						}
					}
				}
			}
		}

		* intensity among the poor
		foreach k of local cutoff {
			foreach w in exo freq pca {
				foreach s of local sim { // simulated indicator
					foreach p of numlist 1/19 {
						foreach r of numlist 1/`=`rep'' {
							gen cens_c_`w'_`k'_s`s'_p`p'_r`r'=c_`w'_s`s'_p`p'_r`r'
							replace cens_c_`w'_`k'_s`s'_p`p'_r`r'=0 if poor_`w'_`k'_s`s'_p`p'_r`r'==0
						}
					}
				}
			}
		}

		* adjusted poverty headcount ratio = the societal poverty measure of interest
		foreach k of local cutoff {
			local j=0
			foreach w in exo freq pca {
				local j=`j'+1
				foreach s of local sim { // simulated indicator
					foreach p of numlist 1/19 {
						foreach r of numlist 1/`=`rep'' {
						* National
							sum cens_c_`w'_`k'_s`s'_p`p'_r`r'
							mat B=(`p',`s',`k',`j',r(mean),0,`r')
							mat R=R\B
							
							* Simulated region
							sum cens_c_`w'_`k'_s`s'_p`p'_r`r' if region==`sreg'
							mat B=(`p',`s',`k',`j',r(mean),`sreg',`r')
							mat R=R\B
						}
					}
				}
			}
		}
		svmat R, names(col)
		keep if k!=.
		keep sim_port sim_ind k w_method mpi subnat rep
		label val w_method method
		gen mc=`x'
		append using "$path/temp1_k`cutoff'"
		save "$path/temp1_k`cutoff'", replace
}


*******************************************
*Graph Creation
*******************************************
	clear all
	set maxvar 120000 
	set matsize 11000 
	
	global path "/Users/ricardonogales/Dropbox/Endogenous Weights/Datasets and dofiles/"
	set scheme s1color
	
	local cutoff 50
	use "$path/temp1_k`cutoff'"
	label def method 1 "Exogenous" 2 "Frequency" 3 "PCA"
	label val w_method method
	drop a // only to initialize
	local sreg 1 // Simulated region
	gen mpi_next=.
	forvalues s=0/19 {
		bys rep mc subnat sim_ind k w_method: gen pair=. if sim_port==`s'
		bys rep mc subnat sim_ind k w_method: replace pair=mpi if sim_port==`s'+1
		bys rep mc subnat sim_ind k w_method: egen aux=max(pair)
		bys rep mc subnat sim_ind k w_method: replace mpi_next=aux if sim_port==`s'
		drop aux pair
	}
	sort rep mc sim_ind k w_method sim_port subnat
	gen dyn_cty=(mpi_next-mpi) if subnat==0
	gen dyn_reg=(mpi_next-mpi) if subnat==`sreg'
	
	bys rep mc sim_ind k w_method sim_port: egen aux=max(dyn_cty) 
	replace dyn_cty=aux if dyn_cty==.
	drop aux
	
	bys rep mc sim_ind k w_method sim_port: egen aux=max(dyn_reg) 
	replace dyn_reg=aux if dyn_reg==.
	drop aux
	
	gen viol=(sign(dyn_cty)!=sign(dyn_reg))
	
	* bar graph
	preserve
		collapse viol, by(sim_port sim_ind w_method)
		foreach i of numlist 1/4 {
			foreach m of numlist 1 2 3 {
				twoway bar viol sim_port if w_method==`m' & sim_ind==`i',  ///
						name(g_`m'_s`i', replace) title("`: label (w_method) `m'', Indicator `i'") nodraw xlabel(1(2)19) ylabel(#10,grid angle(h) format(%9.2g)) ///
						legend(order(1 "Proportion of violations of Subgroup Consistency")) ytitle(%) xtitle(Ventile)			
			}
		}
		
		grc1leg g_1_s1 g_1_s2 g_1_s3 g_1_s4 ///
					  g_2_s1 g_2_s2 g_2_s3 g_2_s4 ///
					  g_3_s1 g_3_s2 g_3_s3 g_3_s4 ///
						, row(3) col(4) ycommon name(g_`k', replace) // title("Poverty cutoff: k=0.`cutoff'")
		graph export "$path/graphs and tables/sim_sgc_v3_k`cutoff'.png", replace
	restore
	
	* line graph
	preserve
		gen mean_mpi=.
		keep if mc==3
		local sreg 1 // Simulated region
		foreach i of numlist 1 2 3 4  {
			foreach m of numlist 1 2 3 {
				foreach r of numlist 0 `sreg' {
					forvalues s=1/19 {
						mean mpi if sim_ind==`i' & w_method==`m' & sim_port==`s' & subnat==`r'
						mat A=e(b)
						replace mean_mpi=A[1,1] if sim_ind==`i' & w_method==`m' & sim_port==`s' & subnat==`r'
					}
				}
			}
		}
		foreach m of numlist 1 2 3 {
			foreach r of numlist 0 `sreg' {
				mean mpi if sim_ind==. & w_method==`m' & subnat==`r'
				mat A=e(b)
				replace mean_mpi=A[1,1] if sim_ind==. & w_method==`m' & subnat==`r'
			}
		}
		drop mpi mpi_next viol dyn_cty dyn_reg
		rename mean_mpi mpi
		
		bys sim_ind w_method sim_port subnat: gen tag=_n
		keep if tag==1
		
		gen mpi_next=.
		forvalues s=0/19 {
			bys mc subnat sim_ind k w_method: gen pair=. if sim_port==`s'
			bys mc subnat sim_ind k w_method: replace pair=mpi if sim_port==`s'+1
			bys mc subnat sim_ind k w_method: egen aux=max(pair)
			bys mc subnat sim_ind k w_method: replace mpi_next=aux if sim_port==`s'
			drop aux pair
		}
		sort rep mc sim_ind k w_method sim_port subnat
		gen dyn_cty=(mpi_next-mpi) if subnat==0
		gen dyn_reg=(mpi_next-mpi) if subnat==`sreg'
		
		bys mc sim_ind k w_method sim_port: egen aux=max(dyn_cty) 
		replace dyn_cty=aux if dyn_cty==.
		drop aux
		
		bys mc sim_ind k w_method sim_port: egen aux=max(dyn_reg) 
		replace dyn_reg=aux if dyn_reg==.
		drop aux
		
		gen viol=(sign(dyn_cty)!=sign(dyn_reg))
		gen bar=mpi_next if viol==1
		gen s_bar=sim_port+1 if viol==1
		
		foreach m of numlist 1 2 3 {
			foreach i of numlist 1 2 3 4 {
					twoway (scatter mpi sim_port if sim_ind==`i' & w_method==`m' & subnat==0, sort ylabel(, angle(h) format(%9.2g) axis(1)) connect(l) mcolor(maroon) lcolor(maroon) yaxis(1)) ///
						   (scatter mpi sim_port if sim_ind==`i' & w_method==`m' & subnat==`sreg', sort ylabel(, angle(h) format(%9.2g) axis(2)) connect(l) mcolor(navy) lcolor(navy) yaxis(2)) ///
						   (bar bar s_bar if sim_ind==`i' & w_method==`m' & subnat==0, color(gray%30) barw(1) yaxis(1)), ///
								name(g_`m'_s`i', replace) ytitle(P, axis(1)) ytitle(P, axis(2)) title("`: label (w_method) `m'', Indicator `i'") nodraw legend(rows(3) order(1 "National (left axis)" 3 "Subnational (right axis)" 2 "Violations of Subgroup Consistency")) 
			}
		}
		grc1leg g_1_s1 g_1_s2 g_1_s3 g_1_s4 ///
					  g_2_s1 g_2_s2 g_2_s3 g_2_s4 ///
					  g_3_s1 g_3_s2 g_3_s3 g_3_s4 ///
						, row(3) col(4) name(g_`k', replace) // title("Poverty cutoff: k=0.`cutoff'")
		graph export "$path/graphs and tables/sim_sgc_example_v3_k`cutoff'.png", replace
	restore

	


			
