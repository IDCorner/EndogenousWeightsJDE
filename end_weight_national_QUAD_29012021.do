*** Baseline MPI Analysis
	clear all
	set maxvar 120000
	global path "/Users/ricardonogales/Dropbox/Endogenous Weights/Datasets and dofiles"
	cd "$path"
	
	foreach country in ECU UGA {
		use mpi_`country'_full, clear
		*use mpi_UGA_full, clear
		svyset psu [w=sweight], strata(strata)
		egen nd=rowtotal(d_*)
		svy: prop nd 
		foreach k of numlist 1/9 {
			gen poor_`k'=(nd>=`k')
		}
		svy: mean poor_*
		svy: mean d_*
		tetrachoric d_*, posdef

		global ind_nutr d_cm d_satt d_educ d_wtr d_elct d_sani d_ckfl d_hsg d_asst
		global ind_wtr d_cm d_satt d_educ d_nutr d_elct d_sani d_ckfl d_hsg d_asst
		global ind_elct d_cm d_satt d_educ d_nutr d_wtr d_sani d_ckfl d_hsg d_asst

		local sims nutr wtr elct
		foreach v of local sims {
			* Generate potential deprivation increasing situations (pop. divided in 40 parts)
			set seed 12345
			sort hh_id
			gen u_`v'=runiform() if d_`v'==0 // Uniform likelihood to be affected
			xtile pct_`v' = u_`v', nq(20) // Potentially affected people (steps of 2.5%)
			foreach p of numlist 1/20 {
				clonevar d_`v'_s_`p'=d_`v'
				replace d_`v'_s_`p'=1 if pct_`v'<=`p' // Deprivation appears
			}
			clonevar d_`v'_s_0=d_`v' // Baseline: real deprivations in data
			svy: mean d_`v'_s*
		
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
				mat colnames W_`v'=sim ind x sv_0 sv_1 sv_2 sv_3 sv_4 sv_5 sv_6 sv_7 sv_8 sv_9 sv_10 ///
								sv_11 sv_12 sv_13 sv_14 sv_15 sv_16 sv_17 sv_18 sv_19 ///
								xx f_0 f_1 f_2 f_3 f_4 f_5 f_6 f_7 f_8 f_9 f_10 ///
								f_11 f_12 f_13 f_14 f_15 f_16 f_17 f_18 f_19
		}
		mat def W=(W_nutr\W_wtr\W_elct)

		preserve
			svmat W, names(col)
			gen country="`country'"
			keep country sim ind sv_0 sv_1 sv_2 sv_3 sv_4 sv_5 sv_6 sv_7 sv_8 sv_9 sv_10 ///
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
	save PCA_w, replace
	erase PCA_w_ECU.dta
	erase PCA_w_UGA.dta

*** Create MPIs
clear all
set maxvar 120000 
gen ccty=""
save "$path/mpi_trial.dta", replace

clear all
local cty ECU UGA
foreach c of local cty {
	use PCA_w.dta, clear
	keep if country=="`c'"
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
		* Generate potential additional deprivations (40 groups)
		set seed 12345
		sort hh_id
		gen u_`v'=runiform() if d_`v'==0 // Uniform likelihood
		xtile pct_`v' = u_`v', nq(20) // potential increased deprivations (steps of 2.5%)
		foreach p of numlist 1/20 {
			clonevar d_`v'_s_`p'=d_`v'
			replace d_`v'_s_`p'=1 if pct_`v'<=`p'
		}
		clonevar d_`v'_s_0=d_`v' // Baseline
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
		replace pca_c_nutr_`s'=(pca_c_nutr_`s')^2 
		egen pca_c_wtr_`s'=rowtotal(pca_w_wtr_*_s`s')
		replace pca_c_wtr_`s'=(pca_c_wtr_`s')^2
		egen pca_c_elct_`s'=rowtotal(pca_w_elct_*_s`s')
		replace pca_c_elct_`s'=(pca_c_elct_`s')^2
		
		egen af_c_nutr_`s'=rowtotal(af_w_nutr_*_s`s')
		replace af_c_nutr_`s'=(af_c_nutr_`s')^2
		egen af_c_wtr_`s'=rowtotal(af_w_wtr_*_s`s')
		replace af_c_wtr_`s'=(af_c_wtr_`s')^2
		egen af_c_elct_`s'=rowtotal(af_w_elct_*_s`s')
		replace af_c_elct_`s'=(af_c_elct_`s')^2
		
		egen f_c_nutr_`s'=rowtotal(f_w_nutr_*_s`s')
		replace f_c_nutr_`s'=(f_c_nutr_`s')^2
		egen f_c_wtr_`s'=rowtotal(f_w_wtr_*_s`s')
		replace f_c_wtr_`s'=(f_c_wtr_`s')^2
		egen f_c_elct_`s'=rowtotal(f_w_elct_*_s`s')
		replace f_c_elct_`s'=(f_c_elct_`s')^2
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

	* Compute MPIs
	local sims nutr wtr elct
	mat C=J(1,18,.)
	local ind=0
	foreach v of local sims {
		local ind=`ind'+1
		foreach k of numlist 0 10(10)50 {
			foreach i of numlist 1/20 {
				local s=`i'-1
				* MPI
				svy: mean pca_cens_c_`v'_0_`k' pca_cens_c_`v'_`s'_`k' 
				_coef_table
				mat A=r(table)'
				mat B1=[1,1,`ind',`k',`s',A[2,1],A[2,2],A[2,5],A[2,6],A[1,1]]					
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
				mat B1=[1,2,`ind',`k',`s',A[2,1],A[2,2],A[2,5],A[2,6],A[1,1]]					
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
				mat B1=[1,3,`ind',`k',`s',A[2,1],A[2,2],A[2,5],A[2,6],A[1,1]]					
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
	preserve
		label def ind 1 "nutr" 2 "wtr" 3 "elct"
		label def meas 1 "MPI" 2 "H"
		label def method 1 "PCA" 2 "AF" 3 "FREQ"
		mat colnames C=meas method ind k s b se ul ll b_OG adif adif_se adif_lb adif_ub rdif rdif_se rdif_lb rdif_ub
		svmat C, names(col)
		keep meas method ind k s b se ul ll b_OG adif adif_se adif_lb adif_ub rdif rdif_se rdif_lb rdif_ub
		gen ccty="`c'"
		order ccty, first
		keep if b!=.
		label val ind ind
		label val method method
		label val meas meas
		append using "$path/mpi_trial.dta", force
		save "$path/mpi_trial.dta", replace
	restore
}

*******************
*** Analysis
*******************

 cd "$path"

* 1. Correlations and frequencies
	set scheme s1color
	use mpi_UGA_full, clear
	gen ccty="UGA"
	tostring psu, gen(psu1)
	drop psu
	rename psu1 psu
	
	append using mpi_ECU_full
	replace ccty="ECU" if ccty==""
	label var d_cm "cm"
	label var d_nutr "nutr"
	label var d_satt "satt"
	label var d_educ "educ"
	label var d_elct "elct"
	label var d_wtr "wtr"
	label var d_sani "sani"
	label var d_ckfl "ckfl"
	label var d_hsg "hsg"
	label var d_asst "asst"
	
	svyset psu [w=sweight], strata(strata)
	egen nd=rowtotal(d_*)
	svy: prop nd 
	foreach k of numlist 1/9 {
		gen poor_`k'=(nd>=`k')
		label var poor_`k' "`k'"
	}
	twoway (hist nd if ccty=="ECU", color(navy%60)  discrete ylabel(#10, grid angle(h))) ///
		   (hist nd if ccty=="UGA", color(maroon%60)  discrete ylabel(#10, grid angle(h))), ///
			ytitle(Frequency(%)) legend(order(1 "ECU" 2 "UGA")) xlabel(#10) xtitle("Number of simultaneous deprivations")
	*graph export "$path/graphs and tables/hist.png", replace
	
	order d_nutr d_cm d_educ d_satt d_ckfl d_sani d_wtr d_elct d_hsg d_asst

		graph bar d_* if ccty=="UGA" [w=sweight], blabel(bar, format(%9.2f)) ylabel(#10, grid angle(h)) legend(rows(2)) ///
				bar(1,color(gs8) lcolor(black)) ///
				bar(2,color(brown) lcolor(black)) ///
				bar(3,color(cranberry) lcolor(black)) ///
				bar(4,color(dkgreen) lcolor(black)) ///
				bar(5,color(dknavy) lcolor(black)) ///
				bar(6,color(dkorange) lcolor(black)) ///
				bar(7,color(emerald) lcolor(black)) ///
				bar(8,color(lavender) lcolor(black)) ///
				bar(9,color(olive_teal) lcolor(black)) ///
				bar(10,color(erose) lcolor(black)) ///
				title(UGA) ytitle(%) legend(order(1 "nutr" 2 "cm" 3 "educ" 4 "satt" 5 "ckfl" 6 "sani" 7 "wtr" 8 "elct" 9 "hsg" 10 "asst")) name(un_UGA, replace)
				
		graph bar d_* if ccty=="ECU" [w=sweight], blabel(bar, format(%9.2f)) ylabel(#10, grid angle(h)) legend(rows(2)) ///
				bar(1,color(gs8) lcolor(black)) ///
				bar(2,color(brown) lcolor(black)) ///
				bar(3,color(cranberry) lcolor(black)) ///
				bar(4,color(dkgreen) lcolor(black)) ///
				bar(5,color(dknavy) lcolor(black)) ///
				bar(6,color(dkorange) lcolor(black)) ///
				bar(7,color(emerald) lcolor(black)) ///
				bar(8,color(lavender) lcolor(black)) ///
				bar(9,color(olive_tile) lcolor(black)) ///
				bar(10,color(erose) lcolor(black)) ///
				title(ECU) ytitle(%) legend(order(1 "nutr" 2 "cm" 3 "educ" 4 "satt" 5 "ckfl" 6 "sani" 7 "wtr" 8 "elct" 9 "hsg" 10 "asst")) name(un_ECU, replace)
	
		grc1leg un_ECU un_UGA, ycommon
		*graph export "$path/graphs and tables/uncen1.png", replace

		
	tetrachoric d_* if ccty=="ECU", posdef
	mat T_ECU=r(Rho)
	*esttab matrix(T_ECU, fmt(%9.3f)) using "$path/graphs and tables/T_ECU1.tex", replace 
	tetrachoric d_* if ccty=="UGA", posdef
	mat T_UGA=r(Rho)
	*esttab matrix(T_UGA, fmt(%9.3f)) using "$path/graphs and tables/T_UGA1.tex", replace

	
* 2. Changes in Societal Poverty
clear all
use "$path/mpi_trial.dta"
set scheme s1color

gen s1=s*5
drop s
rename s1 s
label def method 2 "Exogenous" 1 "PCA" 3 "Frequency", modify
label val method method 

// Violation indicator
preserve
	keep if meas==1
	keep if k==0
	gen b_next=.
	forvalues s=0(5)90 {
		bys ccty method ind: gen pair=. if s==`s'
		bys ccty method ind: replace pair=b if s==`s'+5
		bys ccty method ind: egen aux=max(pair)
		bys ccty method ind: replace b_next=aux if s==`s'
		drop aux pair
	}
	sort ccty method ind s
	gen viol=(b_next-b<0)
	*br ccty method ind s b b_next viol // aux pair


	levelsof method, local(met)
	foreach m of local met {
		levelsof ccty, local(cou)
		foreach c of local cou {
			local ind 1 2 3
			foreach i of local ind {
				gen bar=b_next if viol==1
				gen s_bar=s+5 if viol==1
				twoway (scatter b s if k==0 & ind==`i' & ccty=="`c'" & method==`m' & meas==1, legend(off) connect(l) mcolor(maroon) lcolor(maroon)) ///
					   (bar bar s_bar if k==0 & ind==`i' & ccty=="`c'" & method==`m' & meas==1, color(gray%30) barw(5)), ///
					ylabel(#5, grid) ytitle("P") legend(off) ///
						xtitle(s) xlabel(#5) name(`c'_i_`i'_`m', replace) nodraw title("`: label (method) `m''") legend(order(1 "Mean point estimate" 2 "Violation to Monotonicity"))
				drop bar s_bar
			}
		}
	}

	grc1leg ECU_i_1_2 ECU_i_1_1 ECU_i_1_3, row(1) name(a1, replace) title(Ecuador)
		*graph export "$path/graphs and tables/mon_ECU_nutr.png", replace
	grc1leg ECU_i_2_2 ECU_i_2_1 ECU_i_2_3, row(1) name(b1, replace) title(Ecuador) subtitle(Simulated indicator: Water)
		*graph export "$path/graphs and tables/mon_ECU_wrt.png", replace
	grc1leg ECU_i_3_2 ECU_i_3_1 ECU_i_3_3, row(1) name(c1, replace) title(Ecuador) subtitle(Simulated indicator: Electricity)
		*graph export "$path/graphs and tables/mon_ECU_elct.png", replace

	grc1leg UGA_i_1_2 UGA_i_1_1 UGA_i_1_3, row(1) name(a2, replace)  title(Uganda)
		*graph export "$path/graphs and tables/mon_UGA_nutr.png", replace
	grc1leg UGA_i_2_2 UGA_i_2_1 UGA_i_2_3, row(1) name(b2, replace) title(Uganda) subtitle(Simulated indicator: Water)
		*graph export "$path/graphs and tables/mon_UGA_wtr.png", replace
	grc1leg UGA_i_3_2 UGA_i_3_1 UGA_i_3_3, row(1) name(c2, replace) title(Uganda) subtitle(Simulated indicator: Electricity)
		*graph export "$path/graphs and tables/mon_UGA_elct.png", replace

	grc1leg a1 a2, name(comb, replace) rows(2)
		*graph export "$path/graphs and tables/mon_comb_nutr.png", replace
restore

* 3. Actual Weights

	use PCA_w, clear
		
	mat af=[0.16666667\0.16666667\0.16666667\0.16666667\0.0555556\0.0555556\0.0555556\0.0555556\0.0555556\0.0555556 ]
	local country ECU UGA
	foreach c of local country {
		foreach i of numlist 1 4 5 {
			foreach s of numlist 0 10 19 {	
				mean sv_`s' if country=="`c'" & sim==`i', over(ind)
				mat pca_`c'_`i'_`s'=e(b)'
				
				mean f_`s' if country=="`c'" & sim==`i', over(ind)
				mat f_`c'_`i'_`s'=e(b)'
			}
		}
	}
	
	local country ECU UGA
	foreach c of local country {
		foreach i of numlist 1 4 5 {
			mat e_`c'_`i'=(af,pca_`c'_`i'_0,pca_`c'_`i'_10,pca_`c'_`i'_19,f_`c'_`i'_0,f_`c'_`i'_10,f_`c'_`i'_19)
			mat colnames e_`c'_`i'=af pca_0 pca_10 pca_20 f_0 f_10 f_20
			*esttab matrix(e_`c'_`i', fmt(%9.3f)) using "$path/graphs and tables/w_`c'_`i'.tex", replace
		}
	}
	

