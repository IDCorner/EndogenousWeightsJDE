*** This code creates the MPI variables 
*cd "C:\Users\qehs1094\Dropbox\Mismatch Paper"
*cd "/Volumes/GoogleDrive/Mi unidad/Mismatch paper/Datasets/"
cd "C:\Users\qehs1094\Downloads"

clear all
**** MDP INDICATORS 
* The deprivation matrix is created here.
* Let us use the standard notation d_cm d_nutr d_satt d_educ d_elct d_wtr d_sani d_ckfl d_hsg d_asst
*	for the coding. If this is uniform, it will make it easier to do combined analyses or common codes for all countries
label def dep 1 "Deprived" 0 "Non Deprived"

** 1. Child Mortality// Question asked to women aged 12-49. We cannot make the cut 'in the last 5 years'

//Women's Questionnaire
*use UGA_2013_UNPS_v01_M_STATA8/WSEC4.dta, clear
use WSEC4.dta, clear
gen dead_ch=(w4q8==1) if w4q1==1
replace dead_ch=0 if w4q1==2
*drop if w4q1==.
duplicates tag PID, gen(aaa)
drop if aaa==1

sort PID
save temp_wom.dta, replace
use Uganda2013_All, clear
*duplicates tag PID, gen(aaa)
*drop if aaa==1
merge m:1 PID using temp_wom.dta, force
drop _merge

clonevar age=h2q8

gen woman_1549=(h2q3==2 & age>=15 & age<=49) if h2q3!=. & age!=. 
bys HHID: egen hh_woman_1549=sum(woman_1549)
				
gen cm1=1 if (dead_ch>0 & dead_ch<.) & woman_1549==1 // Woman aged 12-49 lost a child at some point
bys HHID: egen hh_cm1=sum(cm1)

replace hh_cm1=0 if hh_woman_1549==0

gen d_cm=(hh_cm1>0 & hh_cm1<.) if hh_cm1!=.
tab d_cm
label val d_cm dep

** 2. Nutrition // We will use the WHO definitions
**                 Source: http://www.who.int/childgrowth/software/en/
**                 Original definition: Any person under 70 years of age for whom there is nutritional information is undernourished. 
** 					Here we have info for children between 6 months and 59 months.
*adopath + "T:\GMPI 2.0\ado\igrowup_stata"
*adopath + "/Users/ricardonogales/igrowup_stata"
adopath + "C:\Users\qehs1094\igrowup_stata"

/* We use 'reflib' to specify the package directory where the .dta files 
containing the WHO Child growth Standards are stored. Note that we use 
strX to specity the length of the path in string. If the path is long, 
you may specify str55 or more, it won't run. */	
*gen str100 reflib = "T:\GMPI 2.0\ado\igrowup_stata"
gen str100 reflib = "C:\Users\qehs1094\igrowup_stata"
*gen str100 reflib="/Users/ricardonogales/igrowup_stata"
lab var reflib "Directory of reference tables"

/* We use datalib to specify the working directory where the input STATA data
set containing the anthropometric measurement is stored. */
gen str100 datalib = "C:\Users\qehs1094\Desktop"
*gen str100 datalib = "/Users/ricardonogales/Desktop/" 
lab var datalib "Directory for datafiles"

/* We use datalab to specify the name that will prefix the output files that 
will be produced from using this ado file (datalab_z_r_rc and datalab_prev_rc)*/
gen str30 datalab = "children_nutri_uga" 
	//Please change the country name extension //
lab var datalab "Working file"

clonevar gender = h2q3
gen age_months=h6q4
gen  str6 ageunit = "months" 

gen	weight=h6q27a if age_months!=.
gen	height=h6q28a if age_months!=. & h6q28a!=. // Length lying down.
replace height=h6q28b if age_months!=. & h6q28b!=. // Height standing up.

gen measure = "l" if age_months!=. & h6q28a!=. // Length lying down.
	//Child measured lying down
replace measure = "h" if age_months!=. & h6q28b!=. // Length lying down.
	//Child measured standing up
replace measure = " " if age_months!=. & h6q28a==. & h6q28b==.
	//Replace with " " if unknown
tab measure

gen  oedema = "n" //We assume no-one has oedema
gen  sw =. // To be checked. I don't believe we have sampling weights

/*We now run the command to calculate the z-scores with the adofile */

igrowup_restricted reflib datalib datalab gender age_months ageunit weight height measure oedema sw

/*We now turn to using the dta file that was created and that contains 
the calculated z-scores */

use "C:\Users\qehs1094\Desktop\children_nutri_uga_z_rc.dta", clear 
*use "/Users/ricardonogales/Desktop/children_nutri_uga_z_rc.dta", clear
gen	z_scorewa = _zwei
replace z_scorewa = . if _fwei==1 
lab var z_scorew "z-score weight-for-age WHO"

** This indicators relate to children aged <5 years
gen	underwa2 = (z_scorewa < -2.0) if age_months<=60
	//Takes value 1 if the child is under 2 stdev below the median & 0 otherwise
replace underwa2 = . if z_scorewa==.
lab var underwa2  "Child is undernourished (weight-for-age) 2sd - WHO"

gen stunting = (_zlen < -2.0) if age_months<=60
replace stunting = . if _zlen == . | _flen==1
lab var stunting "Child is stunted (length/height-for-age) 2sd - WHO"

gen wasting = (_zwfl < - 2.0) if age_months<=60
replace wasting = . if _zwfl == . | _fwfl == 1
lab var wasting  "Child is wasted (weight-for-length/height) 2sd - WHO"

** These indicators relate to children aged +5years
replace _zbmi=. if _fbmi==1 // Incoherent values
gen low_bmi = (_zbmi < -2.0) if age_months>60 & age_months<.
replace low_bmi=. if _zbmi==.
	//Takes value 1 if the child is under 2 stdev below the median & 0 otherwise

** In this case, we have info for children aged +6 months and <9 years
gen child_elig1=(age_months<=60)
bys HHID: egen hh_child_elig1=total(child_elig1)

** This group is FOR THE MOMENT considered deprived if they are either Stunted or Underweight
gen ch_malnour1=1 if underwa2==1 | stunting==1
replace ch_malnour1=0 if underwa2==0 & stunting==0
by HHID: egen hh_ch_malnour1=total(ch_malnour1)

tab hh_ch_malnour1

gen d_nutr=((hh_ch_malnour1>0 & hh_ch_malnour1<.))
tab d_nutr

label val d_nutr dep

** 3. School Attendance // Any school-aged child is not attending school up to the age at which he/she would complete class 8
**                        Note that according to UNESCO (http://data.uis.unesco.org/Index.aspx?queryid=219), compulsory schooling starts at age 7
**                        And that the theoretical duration is 8 years
gen age_714=(age>6 & age<15)
bys HHID: egen hh_age_714=total(age_714)

gen not_atten=1 if h4q5==1 & age_714==1
replace not_atten=0 if h4q5==2 & age_714==1
replace not_atten=0 if h4q5==3 & age_714==1 // Currently attending  
bys HHID: egen hh_not_atten=total(not_atten) // Number of children not attending school
gen d_satt=1 if hh_not_atten>0 & hh_not_atten<.
replace d_satt=0 if hh_not_atten==0
replace d_satt=0 if hh_age_714==0 // no reference population=no risk of deprivation

tab d_satt, m

label val d_satt dep

** 4. Education // No household member aged ten years or older has completed sixyears of schooling
codebook h4q7, tab(100) // This is only for those who are not currently attending
/*
                            45        10  Some schooling but not Completed
                                          P.1
                           126        11  Completed P.1
                           307        12  Completed P.2
                           486        13  Completed P.3
                           532        14  Completed P.4
                           622        15  Completed P.5
                           842        16  Completed P.6
                           996        17  Completed P.7
                            24        21  Completed J.1
                            35        22  Completed J.2
                            12        23  Completed J.3
                           193        31  Completed S.1
                           298        32  Completed S.2
                           227        33  Completed S.3
                           478        34  Completed S.4
                            19        35  Completed S.5
                           158        36  Completed S.6
                           213        41  Completed Post primary
                                          Specialized training or
                                          Certificate
                           358        51  Completed Post secondary
                                          Specialized training or diploma
                           176        61  Completed Degree and above
                            66        99  DK
                        11,329         .  
*/

gen age_10plus=(age>9 & age<.)
bys HHID: egen hh_age_10plus=total(age_10plus), m // Just to check hh composition

gen educ_6plus=(h4q7>16 & h4q7<62) if h4q7!=.  & age_10plus==1 // Many of these missings are currently attending

codebook h4q10, tab(100) // This is only for those who are currently attending

/*
                           848         1  Attending nursery, kindergarten
                                          etc (lower than P.1)
                           882        10  Attending P.1
                           781        11  Attending P.2
                           772        12  Attending P.3
                           730        13  Attending P.4
                           665        14  Attending P.5
                           576        15  Attending P.6
                           416        16  Attending P.7
                           223        30  Attending S.1
                           251        31  Attending S.2
                           222        32  Attending S.3
                           231        33  Attending S.4
                            79        34  Attending S.5
                           109        35  Attending S.6
                            65        40  Attending post primary/junior
                                          specialized training or
                                          certificate or diploma
                            82        50  Attending Post secondary
                                          Specialized training or diploma
                           122        61  Attending Degree and above
                             2        88  Did not attend last year
                             2        99  Don't Know
                        10,484         .  

*/

replace educ_6plus=1 if educ_6plus==. & (h4q10>29 & h4q10<62) & age_10plus==1
replace educ_6plus=0 if educ_6plus==. & h4q10<30 & age_10plus==1 // Are attending but not the right level

replace educ_6plus=0 if h4q5==1 & age_10plus==1 // Never attended school. 
bys HHID: egen hh_educ_6plus=total(educ_6plus)

gen d_educ=(hh_educ_6plus==0)
replace d_educ=0 if hh_age_10plus==0
replace d_educ=. if hh_educ_6plus==.

tab d_educ

label val d_educ dep

** 5. Electricity // Deprived if they don't have access to electricity
tab h10q1, m
gen d_elct=(h10q1==2) if h10q1!=.

tab d_elct, m

label val d_elct dep

** 6. Water //  They do not state explicitly which sources are improved, so we follow standard SDG guidelines
codebook h9q7, tab(100)
/*
                           536      ND  10  Piped water into the dwelling
                         1,168      ND  11  Piped water to the yard
                         1,632      ND  12  Public tap
                           118      ND  13  Borehole in yard/plot
                         6,219      ND  14  Public borehole
                         3,260      ND  15  Protected well/spring
                         3,006      D   16  Unprotected well/spring
                           851      D   17  River/ Stream/ Lake
                           135      D   18  Vendor
                            20      D   19  Tanker truck
                           263      D   20  Gravity flow scheme
                           189      D   21  Rain water
                            19      D   22  Bottled water
                            94      D  96  Other (specify)
                            32         .  

*/
gen imp_wtr=(h9q7<16) if h9q7!=.

gen d_wtr=1 if imp_wtr==0 
replace d_wtr=0 if  imp_wtr==1 

replace d_wtr=1 if d_wtr==0 & (h9q9a>30 & h9q9a<.) // More than 30 minute walk

label val d_wtr dep

tab d_wtr, m 

** 7. Sanitation // The household’s sanitation facility is not improved (SDG) or it is improved but shared with other households (see p. 23 Report)

codebook h9q22, tab(100)
/*
                           288      ND  10  Flush Toilet
                           752      ND  11  VIP Latrine
                         4,572      ND  12  Covered Pit latrine with slab
                         7,718      D   13  Covered Pit latrine without slab
                           423      D   14  Uncovered Pit latrine with slab
                         2,330      D   15  Uncovered Pit latrine without
                                           slab
                            24      ND  16  Ecosan (compost toilet)
                         1,345      D   17  No facility/Bush/Polythene
                                           bags/Bucket
                            58      ND  96  Other (specify)
                            32         .  
*/
		
clonevar d_sani=h9q22
recode d_sani (13=1)(14=1)(15=1)(17=1) ///
			 (10=0)(11=0)(12=0)(16=0)(96=1)
			 
replace d_sani=1 if d_sani==0 & h9q22a==1 // shared facility
label val d_sani dep

tab d_sani

** 8. Cooking fuel // The household cooks with dung, wood or charcoal.

codebook h10q9, tab(100)
/*
                            42         1  Electric
                            15         2  LPG
                            40         3  Kerosene
                            96         4  Wood / Sawdust Burning
                           129         5  Efficient Wood Burning
                         3,109         6  Charcoal
                            20         7  Other Biomass Burning
                        13,915         8  Open fire
                           101         9  Other (specify)
                            75         .  

*/

clonevar d_ckfl=h10q9
recode d_ckfl (1=0)(2=0)(3=0)(5=0) ///
		(4=1)(6=1)(7=1)(8=1)(9=1)
		
label val d_ckfl dep

tab d_ckfl, m // This is increadibly high

** 9. Housing // The household has inadequate housing: the floor is of natural materials or the roof or wall are of rudimentary materials
codebook h9q5, tab(100)
/*
                           144        10  Concrete/stone
                           361        11  Cement blocks
                         7,337        12  Burnt/Stabilised bricks
                           499        13  Unburnt bricks with cement
                         3,384        14  Unburnt bricks with mud
                            69        15  Wood
                         5,488        16  Mud and Poles
                            60        17  Tin/Iron sheets
                           168        96  Other (specify)
                            32         .  


*/
clonevar wall=h9q5
recode wall (13=1)(14=1)(15=1)(16=1) ///
		(10=0)(11=0)(12=0)(17=0)(96=0)

codebook h9q4, tab(100)
/*
                        12,293        10  Iron sheets
                            51        11  Tiles
                            27        12  Asbestos
                            13        13  Concrete
                             7        14  Tin
                         5,079        15  Thatch
                            40        96  Other (specify)
                            32         .  

*/

clonevar roof=h9q4
recode roof (10=0)(11=0)(13=0)(14=0) ///
		    (12=1)(15=1)(96=1)

codebook h9q6, tab(100)
/*
                           519        10  Concrete
                             4        11  Bricks
                            25        12  Stone
                         5,153        13  Cement screed
                        11,576        14  Rammed earth
                            37        15  Wood
                           186        16  Tiles
                            10        96  Other (specify)
                            32         .  
*/

clonevar floor=h9q4		
recode floor (14=1)(15=1)(96=1) ///
			 (10=0)(11=0)(12=0)(13=0)(16=0)
gen d_hsg=(wall==1 | roof==1 | floor==1)
replace d_hsg=. if wall==1 & roof==. & floor==.

label val d_hsg dep
tab d_hsg


** 10. Assets // The household does not own more than one of these assets: radio, TV, telephone, computer, animal cart, 
**               bicycle, motorbike or refrigerator, and does not own a car or truck
gen radio=(AssetDum7==1 & AssetDum7<.)
gen tv=(AssetDum6==1 & AssetDum6<.)
gen phone1=(AssetDum16==1 & AssetDum16<.)
gen animal_cart=. // Anumal cart will not be taken into account
gen bike=(AssetDum10==1 & AssetDum10<.)
gen motorbike=(AssetDum11==1 & AssetDum11<.)
gen refri=(AssetDum5==1 & AssetDum5<.) // Not sure because it is mixed with other appliances

gen car=(AssetDum12==1 & AssetDum12<.)

egen assets=rowtotal(radio tv phone1 bike motorbike refri) // Anumal cart is not taken into account
gen d_asst=(assets<2)
replace d_asst=. if assets==.
replace d_asst=0 if car>0 & car<. // Car is a veto

label val d_asst dep
tab d_asst, m

			
**** INDIVIDUAL DEPRIVATIONS AND POVERTY STATUS
* We will assume equal dimensional weights
* and equal wights within each dimension
	global ind cm nutr satt educ elct wtr sani hsg ckfl asst
	
	global ind_m d_cm d_nutr d_satt d_educ d_elct d_wtr d_sani d_hsg d_ckfl d_asst
	egen extra=rowmiss($ind_m)
	mdesc $ind_m
	keep if extra==0 // Only keep individuals with full information
	
* Ideally, we would have 10 indicators for every country. 
* If, however, some indicators are missing, then the remaining one(s) within the dimenson
* Absorbe its weight. E.g., if child mortality is missing, then nutrition would get 1/3 weight 
* 	instead of 1/6. If, say, cooking fuel was missing, then the other 5 indicators in the living standards dimension
*	would get a 1/15 weight instead of 1/18
* 	For this in the dataprep file the missing indicators can be coded as a column of dots "."
	global ind_16 cm nutr satt educ // Variables with 1/6 weight
	global ind_118 elct wtr sani hsg ckfl asst // Variables with 1/18 weight
	* Just as examples:
	global ind_m_elct wtr sani hsg ckfl asst // In case of missing electricity
	global ind_m_hsg elct wtr sani ckfl asst // In case of missing housing
	global ind_m_ckfl elct wtr sani hsg asst // In case of missing cooking fuel
	foreach v in $ind_16 {
		gen w_`v'=.
		replace w_`v'=1/6 if d_`v'!=.
	}
	foreach v in $ind_118 {
		gen w_`v'=.
		replace w_`v'=1/18 if d_`v'!=.
	}
	replace w_nutr=1/3 if d_cm==. 
	replace w_cm=1/3 if d_nutr==. 
	replace w_educ=1/3 if d_satt==.
	foreach v in $ind_m_elct {
		replace w_`v'=1/15 if d_elct==.
	}
	foreach v in $ind_m_hsg {
		replace w_`v'=1/15 if d_hsg==.
	}
	foreach v in $ind_m_ckfl {
		replace w_`v'=1/15 if d_ckfl==.
	}
	foreach v in $ind {
		count if w_`v'==.
	}
* In the rest of the code, we should only focus on individual level analysis, so I omit the creation of MDP aggregates
* The reason is that this empirical analysis will only feed the `identification' step of poverty measurement, i.e. who is poor?
* and not the `aggregate' step, i.e. How much poverty is there in the country?
* create c-vector
	foreach v in $ind {
		gen w_d_`v'=d_`v'*w_`v'
	}
	egen c_vector=rowtotal(w_d*) // This is the underlying welfare variable in MDP analysis
	* NOTE THERE WAS A poor_13 variable created
	cap drop poor_13
* Identify the poor for different values of k
	foreach k of numlist 1(1)100 {
		gen poor_`k'=(c_vector>=`k'/100)
	}
* Censored c-vector
	foreach k of numlist 1(1)100  {
		gen cens_c_vector_`k'=c_vector
		replace cens_c_vector_`k'=0 if (c_vector<`k'/100) // This could be interesting if we want to focus on `the poor', if we with to identify them
	}
* Accute, extreme and vulnerable
	gen mdp_1=poor_33 // Binary indicator for accutely poor
	gen mdp_2=poor_50 // Binary indicator for Severely poor
	gen mdp_3=(c_vector>=0.2 & c_vector<0.3333) // Binary indicator for  Vulnerable to MDP
	
	gen cens_c_vector_mdp_1=cens_c_vector_33
	gen cens_c_vector_mdp_2=cens_c_vector_50
	gen cens_c_vector_mdp_3=c_vector
		replace cens_c_vector_mdp_3=0 if (c_vector<0.2 | c_vector>=0.3333)
* Censored deprivation indicators, by type of MDP
	foreach i in $ind{
		foreach mdp of numlist 1 2 3 {
			gen c_d_`i'_mdp_`mdp'=d_`i'
			replace c_d_`i'_mdp_`mdp'=0 if mdp_`mdp'==0
			*gen w_c_`i'_mdp_`mdp'=c_d_`i'_mdp_`mdp'*w_`i'
		}
	}
	
** Consumption
*Uganda PPP
gen PPP2012 = 1045.315902287
gen PPP2013 = 1080.759511214
gen PPP2014 = 1096.217413042
gen PPP2015 = 1156.121051707

*replace PPP = PPP2013
gen PPP = PPP2013

*Consumption Variables Per Person Per Day, Nominal and PPP adjusted.
gen Cons_Nom = (nrrexp30 / 30) /  hsize
gen Cons_PPP = Cons_Nom / PPP

* Additional covariates
by HHID: gen hh_size=_N
gen male=(gender==1)
gen int_year=year
gen int_month=month
gen sweight=wgt_X
gen strata=region
gen psu=h1aq4a
gen hh_id=HHID

keep hh_id sweight strata psu d_* region age hh_size male urban 
cd "C:\Users\qehs1094\Dropbox\Endogenous Weights\Datasets and dofiles"
save mpi_UGA_full, replace

