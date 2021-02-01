*** Baseline MPI Analysis
	clear all
	set maxvar 120000
	global path "/Users/ricardonogales/Dropbox/Endogenous Weights/Datasets and dofiles"
	cd "$path"
	local rs=2
	foreach country in ECU UGA {
		use mpi_`country'_full, clear
		use mpi_ECU_full, clear
		svyset psu [w=sweight], strata(strata)
		egen nd=rowtotal(d_*)
		svy: prop nd 
		foreach k of numlist 1/9 {
			gen poor_`k'=(nd>=`k')
		}
		svy: mean poor_*
		svy: mean d_* // HIGH PREV: Nutrition
						// MID PREV: Water
						// LOW PREV: Electricity
		tetrachoric d_*, posdef

		global ind_nutr d_cm d_satt d_educ d_wtr d_elct d_sani d_ckfl d_hsg d_asst
		global ind_wtr d_cm d_satt d_educ d_nutr d_elct d_sani d_ckfl d_hsg d_asst
		global ind_elct d_cm d_satt d_educ d_nutr d_wtr d_sani d_ckfl d_hsg d_asst

		local sims nutr wtr elct
		foreach v of local sims {
			* Generate potential deprivation increasing situations (pop. divided in 40 parts)
			set seed 12345
			sort hh_id
			if "`country'"=="ECU" {
				gen u_`v'=runiform() if d_`v'==0 & region==`rs' // Uniform likelihood to be affected: AMAZONIA
				xtile pct_`v' = u_`v', nq(20) // Potentially affected people (steps of 2.5%)
				foreach p of numlist 1/20 {
					clonevar d_`v'_s_`p'=d_`v'
					replace d_`v'_s_`p'=1 if pct_`v'<=`p' & region==`rs' // Deprivation appears
				}
				clonevar d_`v'_s_0=d_`v' // Baseline: real deprivations in data
			}
			if "`country'"=="UGA" {
				gen u_`v'=runiform() if d_`v'==0 & region==`rs' // Uniform likelihood to be affected: EASTERN
				xtile pct_`v' = u_`v', nq(20) // Potentially affected people (steps of 2.5%)
				foreach p of numlist 1/20 {
					clonevar d_`v'_s_`p'=d_`v'
					replace d_`v'_s_`p'=1 if pct_`v'<=`p' & region==`rs' // Deprivation appears
				}
				clonevar d_`v'_s_0=d_`v' // Baseline: real deprivations in data
			}
		
			* Region affected
			levelsof region, local(reg)
			foreach r of local reg {
				if `r'==`rs' {
					* PCA weights: Principal Factor Approach
					foreach p of numlist 0/19 { // We affect up to 50% of the non-deprived population, in steps of 2.5%
						tetrachoric ${ind_`v'} d_`v'_s_`p' if region==`r', posdef
						mat C=r(corr)
						local n=_N
						pcamat C, n(`n') factor(1) // we retain only the first principal component
						estat loadings
						mat S1=r(A) // Original loadings = non-normalized weights
						local sum=0
						foreach i of numlist 1/10 {
							local sum=`sum'+S1[`i',1] // Sum of non-normalized weights
						}
						mat S1_b_s_`p'=J(10,1,.)
						foreach i of numlist 1/10 {
							mat S1_b_s_`p'[`i',1]=S1[`i',1]/`sum' // Matrix of normalized weights
						}
					}
					mat S_`v'=J(10,1,.)
					foreach p of numlist 0/19 {
						mat S_`v'=(S_`v',S1_b_s_`p') // Matrix of non-normalized and normalized weights
					} 
					mat R=J(10,1,`r')
					mat S_`v'=(R,S_`v')
						
					* Frequency based: -ln(.) type
					mat F_`v'=J(10,1,.)
					foreach p of numlist 0/19 { // We affect up to 50% of the non-deprived population, in steps of 2.5%
						svy: mean ${ind_`v'} d_`v'_s_`p' if region==`r'
						mat E=e(b)
						local sum=0
						foreach i of numlist 1/10 {
							local sum=`sum'-ln(E[1,`i']) // Sum of -ln(freq)
						}
						mat F1_s_`p'=J(10,1,.)
						foreach i of numlist 1/10 {
							mat F1_s_`p'[`i',1]=-ln(E[1,`i'])/`sum' // Frequency weight
						}
						mat F_`v'=(F_`v',F1_s_`p') // Matrix of frequency weights
					}
				}
			}
		}	
		* Create weights dataset
		label def ind 0 "cm" ///
					1 "nutr" ///
					2 "satt" ///
					3 "sch" ///
					4 "wtr" ///
					5 "elct"  ///
					6 "sani" ///
					7 "ckfl" ///
					8 "hsg" ///
					9 "asst"
		mat O_nutr=(1\1\1\1\1\1\1\1\1\1)
		mat O_wtr= (4\4\4\4\4\4\4\4\4\4)
		mat O_elct=(5\5\5\5\5\5\5\5\5\5)
		mat I_nutr=(0\2\3\4\5\6\7\8\9\1)
		mat I_wtr= (0\2\3\1\5\6\7\8\9\4)
		mat I_elct=(0\2\3\1\4\6\7\8\9\5)
		foreach v in nutr wtr elct {
			mat W_`v'=(O_`v',I_`v',S_`v',F_`v')
				mat colnames W_`v'=sim ind reg x  sv_0 sv_1 sv_2 sv_3 sv_4 sv_5 sv_6 sv_7 sv_8 sv_9 sv_10 ///
								sv_11 sv_12 sv_13 sv_14 sv_15 sv_16 sv_17 sv_18 sv_19 ///
								xx f_0 f_1 f_2 f_3 f_4 f_5 f_6 f_7 f_8 f_9 f_10 ///
								f_11 f_12 f_13 f_14 f_15 f_16 f_17 f_18 f_19
		}
		mat def W_R=(W_nutr\W_wtr\W_elct)
		
		* Whole country
		local sims nutr wtr elct
		foreach v of local sims {
			* PCA weights: Principal Factor Approach
			foreach p of numlist 0/19 { // We affect up to 50% of the non-deprived population, in steps of 2.5%
				tetrachoric ${ind_`v'} d_`v'_s_`p', posdef 
				mat C=r(corr)
				local n=_N
				pcamat C, n(`n') factor(1) // we retain only the first principal component
				estat loadings
				mat S1=r(A) // Original loadings = non-normalized weights
				local sum=0
				foreach i of numlist 1/10 {
					local sum=`sum'+S1[`i',1] // Sum of non-normalized weights
				}
				mat S1_b_s_`p'=J(10,1,.)
				foreach i of numlist 1/10 {
					mat S1_b_s_`p'[`i',1]=S1[`i',1]/`sum' // Matrix of normalized weights
				}
			}
			mat S_`v'=J(10,1,.)
			foreach p of numlist 0/19 {
				mat S_`v'=(S_`v',S1_b_s_`p') // Matrix of non-normalized and normalized weights
			} 
			mat R=J(10,1,0)
			mat S_`v'=(R,S_`v')
				
			* Frequency based: -ln(.) type
			mat F_`v'=J(10,1,.)
			foreach p of numlist 0/19 { // We affect up to 50% of the non-deprived population, in steps of 2.5%
				svy: mean ${ind_`v'} d_`v'_s_`p' 
				mat E=e(b)
				local sum=0
				foreach i of numlist 1/10 {
					local sum=`sum'-ln(E[1,`i']) // Sum of -ln(freq)
				}
				mat F1_s_`p'=J(10,1,.)
				foreach i of numlist 1/10 {
					mat F1_s_`p'[`i',1]=-ln(E[1,`i'])/`sum' // Frequency weight
				}
				mat F_`v'=(F_`v',F1_s_`p') // Matrix of frequency weights
			}
		}

		* Create weights dataset
		mat O_nutr=(1\1\1\1\1\1\1\1\1\1)
		mat O_wtr= (4\4\4\4\4\4\4\4\4\4)
		mat O_elct=(5\5\5\5\5\5\5\5\5\5)
		mat I_nutr=(0\2\3\4\5\6\7\8\9\1)
		mat I_wtr= (0\2\3\1\5\6\7\8\9\4)
		mat I_elct=(0\2\3\1\4\6\7\8\9\5)
		foreach v in nutr wtr elct {
			mat W_`v'=(O_`v',I_`v',S_`v',F_`v')
				mat colnames W_`v'=sim ind reg x  sv_0 sv_1 sv_2 sv_3 sv_4 sv_5 sv_6 sv_7 sv_8 sv_9 sv_10 ///
								sv_11 sv_12 sv_13 sv_14 sv_15 sv_16 sv_17 sv_18 sv_19 ///
								xx f_0 f_1 f_2 f_3 f_4 f_5 f_6 f_7 f_8 f_9 f_10 ///
								f_11 f_12 f_13 f_14 f_15 f_16 f_17 f_18 f_19
		}
		mat def W_C=(W_nutr\W_wtr\W_elct)
		mat W=(W_R\W_C) // Join regional and country matrices

		preserve
			svmat W, names(col)
			gen country="`country'"
			keep country sim ind reg sv_0 sv_1 sv_2 sv_3 sv_4 sv_5 sv_6 sv_7 sv_8 sv_9 sv_10 ///
								sv_11 sv_12 sv_13 sv_14 sv_15 sv_16 sv_17 sv_18 sv_19 ///
								f_0 f_1 f_2 f_3 f_4 f_5 f_6 f_7 f_8 f_9 f_10 ///
								f_11 f_12 f_13 f_14 f_15 f_16 f_17 f_18 f_19
			keep if sv_0!=.
			label val sim ind
			label val ind ind
			save PCA_w_`country', replace
		restore
	}

	use PCA_w_ECU, clear
	append using PCA_w_UGA
	order country sim
	save PCA_w2, replace
	erase PCA_w_ECU.dta
	erase PCA_w_UGA.dta

*** Create MPIs
clear all
gen ccty=""
save "$path/mpi2.dta", replace

local rs=2
clear all
local cty ECU UGA
foreach c of local cty {
	foreach r of numlist 0 `rs' {
		use PCA_w2.dta, clear
		keep if country=="`c'" & reg==`r'
		egen id=concat(sim ind)
		destring id, gen (id_code)
		drop sim ind id
		reshape wide sv* f*, i(country) j(id_code)
		save temp, replace

		use mpi_`c'_full.dta, clear
		gen country="`c'"
		merge m:1 country using temp.dta
		erase temp.dta
		drop _merge

		local sims nutr wtr elct
		foreach v of local sims {
			* Generate potential additional deprivations 
			set seed 12345
			sort hh_id
			if "`c'"=="ECU" {
				gen u_`v'=runiform() if d_`v'==0 & region==`rs' // Uniform likelihood to be affected: AMAZONIA
				xtile pct_`v' = u_`v', nq(20) //
				foreach p of numlist 1/20 {
					clonevar d_`v'_s_`p'=d_`v'
					replace d_`v'_s_`p'=1 if pct_`v'<=`p' // Deprivation appears
				}
				clonevar d_`v'_s_0=d_`v' // Baseline: real deprivations in data
			}
			if "`c'"=="UGA" {
				gen u_`v'=runiform() if d_`v'==0 & region==`rs' // Uniform likelihood to be affected: EASTERN
				xtile pct_`v' = u_`v', nq(20) // 
				foreach p of numlist 1/20 {
					clonevar d_`v'_s_`p'=d_`v'
					replace d_`v'_s_`p'=1 if pct_`v'<=`p' // Deprivation appears
				}
				clonevar d_`v'_s_0=d_`v' // Baseline: real deprivations in data
			}
		}

		order  d_cm d_nutr d_satt d_educ d_wtr d_elct d_sani d_ckfl d_hsg d_asst
		global ind cm nutr satt educ wtr elct sani ckfl hsg asst
		local i=0
		foreach v of global ind {
			rename d_`v' d_`i'
			local i=`i'+1
		}
		
		* AF Weights
		foreach d of numlist 0/3 {
			gen af_w_`d'=1/6
		}
		foreach d of numlist 4/9 {
			gen af_w_`d'=1/18
		}
		
		foreach s of numlist 0/19 { // Only the deprivations in the simulated scenarios are changed, the rest is unchanged
			foreach d of numlist 0/9 {
				if `d'!=1 {
					gen pca_w_nutr_d`d'_s`s'=d_`d'*sv_`s'1`d' 
					gen af_w_nutr_d`d'_s`s'=d_`d'*af_w_`d' 
					gen f_w_nutr_d`d'_s`s'=d_`d'*f_`s'1`d' 
				}
				if `d'==1 {
					gen pca_w_nutr_d`d'_s`s'=d_nutr_s_`s'*sv_`s'1`d' 
					gen af_w_nutr_d`d'_s`s'=d_nutr_s_`s'*af_w_`d' 
					gen f_w_nutr_d`d'_s`s'=d_nutr_s_`s'*f_`s'1`d' 
				}
				if `d'!=4 {
					gen pca_w_wtr_d`d'_s`s'=d_`d'*sv_`s'4`d'
					gen af_w_wtr_d`d'_s`s'=d_`d'*af_w_`d' 
					gen f_w_wtr_d`d'_s`s'=d_`d'*f_`s'4`d' 
				}
				
				if `d'==4 {
					gen pca_w_wtr_d`d'_s`s'=d_wtr_s_`s'*sv_`s'4`d'
					gen af_w_wtr_d`d'_s`s'=d_wtr_s_`s'*af_w_`d' 
					gen f_w_wtr_d`d'_s`s'=d_wtr_s_`s'*f_`s'4`d' 
				}
				
				if `d'!=5 {
					gen pca_w_elct_d`d'_s`s'=d_`d'*sv_`s'5`d'
					gen af_w_elct_d`d'_s`s'=d_`d'*af_w_`d' 
					gen f_w_elct_d`d'_s`s'=d_`d'*f_`s'5`d'
				}
				
				if `d'==5 {
					gen pca_w_elct_d`d'_s`s'=d_elct_s_`s'*sv_`s'5`d'
					gen af_w_elct_d`d'_s`s'=d_elct_s_`s'*af_w_`d' 
					gen f_w_elct_d`d'_s`s'=d_elct_s_`s'*f_`s'5`d'
				}
			}
			egen pca_c_nutr_`s'=rowtotal(pca_w_nutr_*_s`s')
			egen pca_c_wtr_`s'=rowtotal(pca_w_wtr_*_s`s')
			egen pca_c_elct_`s'=rowtotal(pca_w_elct_*_s`s')
			
			egen af_c_nutr_`s'=rowtotal(af_w_nutr_*_s`s')
			egen af_c_wtr_`s'=rowtotal(af_w_wtr_*_s`s')
			egen af_c_elct_`s'=rowtotal(af_w_elct_*_s`s')
			
			egen f_c_nutr_`s'=rowtotal(f_w_nutr_*_s`s')
			egen f_c_wtr_`s'=rowtotal(f_w_wtr_*_s`s')
			egen f_c_elct_`s'=rowtotal(f_w_elct_*_s`s')
		}
		foreach k of numlist 0 10(10)90 {
			foreach s of numlist 0/19 {
				gen pca_poor_nutr_`s'_`k'=(pca_c_nutr_`s'>`k'/100)
				gen pca_poor_wtr_`s'_`k'=(pca_c_wtr_`s'>`k'/100)
				gen pca_poor_elct_`s'_`k'=(pca_c_elct_`s'>`k'/100)
				
				gen af_poor_nutr_`s'_`k'=(af_c_nutr_`s'>`k'/100)
				gen af_poor_wtr_`s'_`k'=(af_c_wtr_`s'>`k'/100)
				gen af_poor_elct_`s'_`k'=(af_c_elct_`s'>`k'/100)
				
				gen f_poor_nutr_`s'_`k'=(f_c_nutr_`s'>`k'/100)
				gen f_poor_wtr_`s'_`k'=(f_c_wtr_`s'>`k'/100)
				gen f_poor_elct_`s'_`k'=(f_c_elct_`s'>`k'/100)
			}
		}

		foreach k of numlist 0 10(10)90 {
			foreach s of numlist 0/19 {
				gen pca_cens_c_nutr_`s'_`k'=pca_c_nutr_`s'
				replace pca_cens_c_nutr_`s'_`k'=0 if pca_poor_nutr_`s'_`k'==0	
				gen pca_cens_c_wtr_`s'_`k'=pca_c_wtr_`s'
				replace pca_cens_c_wtr_`s'_`k'=0 if pca_poor_wtr_`s'_`k'==0			
				gen pca_cens_c_elct_`s'_`k'=pca_c_elct_`s'
				replace pca_cens_c_elct_`s'_`k'=0 if pca_poor_elct_`s'_`k'==0
				
				gen af_cens_c_nutr_`s'_`k'=af_c_nutr_`s'
				replace af_cens_c_nutr_`s'_`k'=0 if af_poor_nutr_`s'_`k'==0
				gen af_cens_c_wtr_`s'_`k'=af_c_wtr_`s'
				replace af_cens_c_wtr_`s'_`k'=0 if af_poor_wtr_`s'_`k'==0
				gen af_cens_c_elct_`s'_`k'=af_c_elct_`s'
				replace af_cens_c_elct_`s'_`k'=0 if af_poor_elct_`s'_`k'==0
				
				gen f_cens_c_nutr_`s'_`k'=f_c_nutr_`s'
				replace f_cens_c_nutr_`s'_`k'=0 if f_poor_nutr_`s'_`k'==0
				gen f_cens_c_wtr_`s'_`k'=f_c_wtr_`s'
				replace f_cens_c_wtr_`s'_`k'=0 if f_poor_wtr_`s'_`k'==0
				gen f_cens_c_elct_`s'_`k'=f_c_elct_`s'
				replace f_cens_c_elct_`s'_`k'=0 if f_poor_elct_`s'_`k'==0
			}
		}
		* Compute regional MPIs
		local sims nutr wtr elct
		mat C=J(1,18,.)
		local ind=0
		foreach v of local sims {
			local ind=`ind'+1
			foreach k of numlist 0 10(10)50 {
				foreach i of numlist 1/20 {
					local s=`i'-1
					if `r'==`rs' {
						svy: mean pca_cens_c_`v'_0_`k' pca_cens_c_`v'_`s'_`k' if region==`rs'
						_coef_table
						mat A=r(table)'
						mat B1=[1,`rs',`ind',`k',`s',A[2,1],A[2,2],A[2,5],A[2,6],A[1,1]] // Second digit stands for region				
						lincom _b[pca_cens_c_`v'_`s'_`k']-_b[pca_cens_c_`v'_0_`k']
						mat B2=[r(estimate),r(se),r(lb),r(ub)]
						nlcom (_b[pca_cens_c_`v'_`s'_`k']-_b[pca_cens_c_`v'_0_`k'])/_b[pca_cens_c_`v'_0_`k'], post
						_coef_table
						mat A1=r(table)'
						mat B3=[A1[1,1],A1[1,2],A1[1,3],A1[1,4]]
						mat B=[B1,B2,B3]
						mat C=[C\B]
						
						svy: mean af_cens_c_`v'_0_`k' af_cens_c_`v'_`s'_`k' if region==`rs'
						_coef_table
						mat A=r(table)'
						mat B1=[2,`rs',`ind',`k',`s',A[2,1],A[2,2],A[2,5],A[2,6],A[1,1]]	 // second digit stands for region				
						lincom _b[af_cens_c_`v'_`s'_`k']-_b[af_cens_c_`v'_0_`k']
						mat B2=[r(estimate),r(se),r(lb),r(ub)]
						nlcom (_b[af_cens_c_`v'_`s'_`k']-_b[af_cens_c_`v'_0_`k'])/_b[af_cens_c_`v'_0_`k'], post
						_coef_table
						mat A1=r(table)'
						mat B3=[A1[1,1],A1[1,2],A1[1,3],A1[1,4]]
						mat B=[B1,B2,B3]
						mat C=[C\B]
						
						svy: mean f_cens_c_`v'_0_`k' f_cens_c_`v'_`s'_`k' if region==`rs'
						_coef_table
						mat A=r(table)'
						mat B1=[3,`rs',`ind',`k',`s',A[2,1],A[2,2],A[2,5],A[2,6],A[1,1]] // Second digit stands for region
						lincom _b[f_cens_c_`v'_`s'_`k']-_b[f_cens_c_`v'_0_`k']
						mat B2=[r(estimate),r(se),r(lb),r(ub)]
						nlcom (_b[f_cens_c_`v'_`s'_`k']-_b[f_cens_c_`v'_0_`k'])/_b[f_cens_c_`v'_0_`k'], post
						_coef_table
						mat A1=r(table)'
						mat B3=[A1[1,1],A1[1,2],A1[1,3],A1[1,4]]
						mat B=[B1,B2,B3]
						mat C=[C\B]
					}
					if `r'==0 {
						svy: mean pca_cens_c_`v'_0_`k' pca_cens_c_`v'_`s'_`k'
						_coef_table
						mat A=r(table)'
						mat B1=[1,0,`ind',`k',`s',A[2,1],A[2,2],A[2,5],A[2,6],A[1,1]] // Second digit stands for national				
						lincom _b[pca_cens_c_`v'_`s'_`k']-_b[pca_cens_c_`v'_0_`k']
						mat B2=[r(estimate),r(se),r(lb),r(ub)]
						nlcom (_b[pca_cens_c_`v'_`s'_`k']-_b[pca_cens_c_`v'_0_`k'])/_b[pca_cens_c_`v'_0_`k'], post
						_coef_table
						mat A1=r(table)'
						mat B3=[A1[1,1],A1[1,2],A1[1,3],A1[1,4]]
						mat B=[B1,B2,B3]
						mat C=[C\B]
						
						svy: mean af_cens_c_`v'_0_`k' af_cens_c_`v'_`s'_`k' 
						_coef_table
						mat A=r(table)'
						mat B1=[2,0,`ind',`k',`s',A[2,1],A[2,2],A[2,5],A[2,6],A[1,1]]	 // second digit stands for national				
						lincom _b[af_cens_c_`v'_`s'_`k']-_b[af_cens_c_`v'_0_`k']
						mat B2=[r(estimate),r(se),r(lb),r(ub)]
						nlcom (_b[af_cens_c_`v'_`s'_`k']-_b[af_cens_c_`v'_0_`k'])/_b[af_cens_c_`v'_0_`k'], post
						_coef_table
						mat A1=r(table)'
						mat B3=[A1[1,1],A1[1,2],A1[1,3],A1[1,4]]
						mat B=[B1,B2,B3]
						mat C=[C\B]
						
						svy: mean f_cens_c_`v'_0_`k' f_cens_c_`v'_`s'_`k'
						_coef_table
						mat A=r(table)'
						mat B1=[3,0,`ind',`k',`s',A[2,1],A[2,2],A[2,5],A[2,6],A[1,1]] // Second digit stands for region
						lincom _b[f_cens_c_`v'_`s'_`k']-_b[f_cens_c_`v'_0_`k']
						mat B2=[r(estimate),r(se),r(lb),r(ub)]
						nlcom (_b[f_cens_c_`v'_`s'_`k']-_b[f_cens_c_`v'_0_`k'])/_b[f_cens_c_`v'_0_`k'], post
						_coef_table
						mat A1=r(table)'
						mat B3=[A1[1,1],A1[1,2],A1[1,3],A1[1,4]]
						mat B=[B1,B2,B3]
						mat C=[C\B]
					}
				}
			}
		}
		preserve
			label def ind 1 "nutr" 2 "wtr" 3 "elct"
			label def method 1 "PCA" 2 "AF" 3 "FREQ"
			label def reg 0 "Country" `rs' "Region `rs'"
			cap drop reg
			mat colnames C=method reg ind k s b se ul ll b_OR adif adif_se adif_lb adif_ub rdif rdif_se rdif_lb rdif_ub
			svmat C, names(col)
			keep method ind reg k s b se ul ll b_OR adif adif_se adif_lb adif_ub rdif rdif_se rdif_lb rdif_ub
			gen ccty="`c'"
			order ccty, first
			keep if b!=.
			label val ind ind
			label val method method
			label val reg reg
			append using "$path/mpi2.dta", force
			save "$path/mpi2.dta", replace
		restore
	}
}

* Analysis
use "$path/mpi2.dta", clear

set scheme s1color
gen s1=s*5
drop s
rename s1 s
label def method 2 "Exogenous" 1 "PCA" 3 "Frequency", modify
label val method method 

preserve
	keep if k==0
	gen b_next=.
	forvalues s=0(5)90 {
		bys ccty reg method ind: gen pair=. if s==`s'
		bys ccty reg method ind: replace pair=b if s==`s'+5
		bys ccty reg method ind: egen aux=max(pair)
		bys ccty reg method ind: replace b_next=aux if s==`s'
		drop aux pair
	}
	sort ccty reg method ind s
	gen dyn_cty=(b_next-b) if reg==0
	gen dyn_reg=(b_next-b) if reg==2
	
	bys ccty method ind s: egen aux=max(dyn_cty) 
	replace dyn_cty=aux if dyn_cty==.
	drop aux
	
	bys ccty method ind s: egen aux=max(dyn_reg) 
	replace dyn_reg=aux if dyn_reg==.
	drop aux
	

	sort ccty reg method ind s
	
	gen viol=(sign(dyn_cty)!=sign(dyn_reg))

	levelsof method, local(met)
	foreach m of local met {
		levelsof ccty, local(cou)
		foreach c of local cou {
			local kval 0
			foreach k of local kval {
				foreach i of numlist 1 2 3 {
					gen bar=b_next if viol==1
					gen s_bar=s+5 if viol==1
					twoway (scatter b s if k==0 & ind==`i' & ccty=="`c'" & method==`m' & reg==0, connect(l) mcolor(maroon) lcolor(maroon) yaxis(1)) ///
						   (scatter b s if k==0 & ind==`i' & ccty=="`c'" & method==`m' & reg==2, connect(l) mcolor(navy) lcolor(navy) yaxis(2)) ///
						   (bar bar s_bar if k==0 & ind==`i' & ccty=="`c'" & method==`m' & reg==0, color(gray%30) barw(5) yaxis(1)), ///
											ytitle("P (National)", axis(1)) ytitle("P (Regional)", axis(2)) ///
						xtitle(s) xlabel(#5) ylabel(#4,grid axis(2)) ylabel(#4,grid axis(1)) name(`c'_i_`i'_m_`m', replace) nodraw title("`: label (method) `m''") ///
						legend(rows(3) order(1 "National (mean point estimate, left axis)" 3 "Regional (mean point estimate, right axis)" 2 "Violation to Subgroup Consistency"))
					drop bar s_bar
				}
			}
		}
	}
		

	grc1leg  ECU_i_1_m_2 ECU_i_1_m_1 ECU_i_1_m_3, row(1) name(a1, replace) title(Ecuador)
		graph export "$path/graphs and tables\cons_ECU_nutr.png", replace
	grc1leg  ECU_i_2_m_2 ECU_i_2_m_1 ECU_i_2_m_3, row(1) name(b1, replace) title(Ecuador) subtitle(Simulated indicator: Water)
		graph export "$path/graphs and tables\cons_ECU_wtr.png", replace
	grc1leg  ECU_i_3_m_2 ECU_i_3_m_1 ECU_i_3_m_3, row(1) name(c1, replace) title(Ecuador) subtitle(Simulated indicator: Electricity)
		graph export "$path/graphs and tables\cons_ECU_elct.png", replace
	grc1leg  UGA_i_1_m_2 UGA_i_1_m_1 UGA_i_1_m_3, row(1) name(a2, replace) title(Uganda)
		graph export "$path/raphs and tables\cons_UGA_nutr.png", replace
	grc1leg  UGA_i_2_m_2 UGA_i_2_m_1 UGA_i_2_m_3, row(1) name(b2, replace) title(Uganda) subtitle(Simulated indicator: Water)
		graph export "$path/graphs and tables\cons_UGA_wtr.png", replace
	grc1leg  UGA_i_3_m_2 UGA_i_3_m_1 UGA_i_3_m_3, row(1) name(c2, replace) title(Uganda) subtitle(Simulated indicator: Electricity)
		graph export "$path/graphs and tables\cons_UGA_elct.png", replace

	grc1leg  ECU_i_1_m_2 ECU_i_1_m_1 ECU_i_1_m_3 UGA_i_1_m_2 UGA_i_1_m_1 UGA_i_1_m_3, row(2) name(a1, replace) 
	grc1leg  a1 a2, row(2) name(sc_comb_nutr, replace) 
		graph export "$path/graphs and tables\sc_comb_nutr.png", replace
restore


