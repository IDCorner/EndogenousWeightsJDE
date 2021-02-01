clear
set more off
cap log close
*****************************
***  ECUADOR 2013 - 2014  ***
*****************************

*** Working Folder Path ***
global workingfolder_in  "C:\Users\qehs1094\Desktop\Ecuador ECV 2013-2014"
global workingfolder_out "C:\Users\qehs1094\Desktop\Ecuador ECV 2013-2014"

** Lower case for all datasets
local files: dir "$workingfolder_in" files "*.dta"
foreach file in `files' {
	cd "$workingfolder_in"
	use `file', clear
	rename *, lower
	save `file', replace
}

*** Merge ***
*------------
clear
*Household's assets database
*---------------------------
use "$workingfolder_in\ecv6r_equipamiento.dta", clear
keep  identif_hog eq00 eqbien  eq01
keep if eq00==4 | eq00==11 | eq00==19 | eq00==25 | eq00==26 | eq00==28 | eq00==29 | eq00==30 | eq00==33 | eq00==34 

replace eqbien="eqpsonido"   if eq00==11
replace eqbien="telfijo"     if eq00==19
replace eqbien="grabadora"   if eq00==25
replace eqbien="TV_bn"       if eq00==28
replace eqbien="TV_plas_lcd" if eq00==29
replace eqbien="TV_color"    if eq00==30
replace eqbien="carro"       if eq00==33
replace eqbien="moto"        if eq00==34
drop eq00
rename eq01 _
reshape wide _, i(identif_hog) j(eqbien) string
sort identif_hog

*dwelling database 
*-----------------
merge 1:1 identif_hog using  "$workingfolder_in\ecv6r_vivienda.dta"
sort identif_hog
drop _merge

*database and individual level 
*----------------------------
merge 1:m identif_hog using "$workingfolder_in\ecv6r_personas.dta"
drop _merge
egen ind_id=concat(identif_hog persona)
**********************************************************************************
*** Step 1: Data preparation 
*** Selecting main variables from IR, MR and KR recode and merging with PR recode 
*********************************************************************************

******************************************************
*** Step 1.1 KR - Children's Recode (Children under 5)
******************************************************
gen child_KR=1 if edad<5

**************************************
*** Step 1.1a Key variable for merging
**************************************

********************************
*** Step 1.1b Children Nutrition
********************************
* Variable "Sex"
gen gender=.
replace gender = sexo if edad<5
label var gender "gender"
label define gender 1 "male" 2"female"
label value gender gender

* Variable "Age"
gen age_months=(edad*12)+pd03
label var age_months "age in months of children <5"
gen str6 ageunit = "months" 
lab var ageunit "Months"

* Variable "body weight" - it must be in kilograms

*correcting some ouliers
gen dif=1 if ps82!=ps82a
replace dif=. if edad>=5
gen dif2= ps82-ps82a if dif==1
gen average=(ps82+ ps82a)/2

replace ps82 =ps82a if ind_id== "08065499900204111"
replace ps82 = average if (dif2<=-0.5) | (dif2>=0.5 & dif2<1.4)
drop dif*

gen	weight = ps82 
replace weight=. if edad>=5 
label var weight "weight in kilogram"

* Variable "height" - it must be in centimeters
*Correting some outliers, when the measures are different we take the average
gen dif=1 if ps83!=ps83a
gen dif2= ps83-ps83a if dif==1
replace ps83 = (ps83 + ps83a)/2 if (dif2<=-0.5) | (dif2>=0.5 & dif2<2.2)
drop dif*

gen dif=1 if ps84!=ps84a
gen dif2= ps84-ps84a if dif==1
replace dif2=. if edad>=5

replace ps84 = ps84a if ind_id== "1409509990170416" | ind_id=="0701520010020614" | ind_id=="0901501230060217" | ind_id=="1601569990130419" | ind_id=="1601569990010919"
replace ps84 = (ps84 + ps84a)/2 if (dif2>=-1.5 & dif2<=-0.5) | (dif2>=0.5 & dif2<2.6)
replace ps84 = ps84a if ind_id=="0901504620051214"
drop dif*

gen	height = ps83 if edad<2
replace  height=ps84 if edad>=2 & edad<5
label var height "height of children under five"

gen	 measure    =1 if ps83!=. & edad<2
replace measure =2 if ps84!=. & edad>=2 & edad<5

* Variable "Oedema" //Only if the variable is available to control for oedema//
gen  oedema = "n"  //It assumes no-one has oedema//

** We now run the command to calculate the z-scores with the adofile
zscore06, s(gender) a(age_months) h(height) w(weight) measure(measure)
*Note: zscore06 assumes that children had no oedema.

gen	z_scorewa = waz06
replace z_scorewa = . if waz06<-6 | waz06>5 
lab var z_scorew "z-score weight-for-age WHO"
/* NOTE: Acording to the following document: http://www.who.int/growthref/tools/readme_stata.pdf
the flags for weight for age are if waz06<-6 | waz06>5
*/

***Now we create the under weight-for-age variable with WHO***

* Standard MPI indicator
gen	underwa2 = (z_scorewa < -2.0) //Takes value 1 if the child is under 2 stdev below the median and 0 otherwise//
replace underwa2 = . if z_scorewa==.
lab var underwa2  "Child is undernourished (weight-for-age) 2sd - WHO"

* Ultrapoverty indicator
gen	underwa3 = (z_scorewa < -3.0) //Takes value 1 if the child is under 3 stdev below the median and 0 otherwise//
replace underwa3 = . if z_scorewa==.
lab var underwa3  "Child is undernourished (weight-for-age) 3sd - WHO"


************************************************************************************************
*** Step 1.2  IR - individual recode (all eligible females between 12-49 years in the household)
************************************************************************************************
* all the woman variables needed are in the databases

******************************************************************************************************
*** Step 1.3  MR - Male Recode (all eligible man - ages varies from survey to survey)
******************************************************************************************************
* all the man variables needed are in the databases 

****************************************
*** Step 1.4  PR Household Member Recode
****************************************
gen country = "Ecuador" 
gen year    = 2014  /*The year of the survey is November,2013- Octubre,2014*/
gen survey  = "ECV"

*** Generate a household unique key variable required to prepare the MPI indicators at the household level ***
gen hh_id =  identif_hog
label var  hh_id "Household ID"

*** Generate individual unique key variable required for data merging
label var  ind_id "Individual ID"

***************************************************************
*** Step 1.5  DATA MERGING (PR, IR, MR, KR) & Control Variables
***************************************************************

*********************
*** 1.5a DATA MERGING
*********************

*NOTE: it is not needed to merge databases.

*** Matrix for quality checks 
mat qcheck = J(1,96,.)

**********************************************************************************************
*** 1.5b Create control variables if household has 'Eligible' members for IR, MR and KR Recode
**********************************************************************************************

*** Eligible Women for IR Recode ***
gen	fem_eligible =  ((edad>=12 & edad<=49) & sexo==2)
replace fem_eligible =. if edad==. | sexo==.
bys	hh_id: egen hh_n_fem_eligible = sum(fem_eligible) //Number of eligible women for interview in the hh//
gen	no_fem_eligible = (hh_n_fem_eligible==0) //Takes value 1 if the household had no eligible females for an interview//
lab var no_fem_eligible "Household has no eligible women"

*** Eligible Women for anthropometric measurements ** 
*As we have information for women beyond 49 years old, we augment the eligible age to 70 
gen	fembmi_eligible = ((edad>=12 & edad<=70) & sexo==2)
replace fembmi_eligible = . if edad==. | sexo==.  
bys hh_id: egen hh_n_fembmi_eligible = sum(fembmi_eligible) //Number of eligible women for hemoglobin measurements//
gen	no_fembmi_eligible = (hh_n_fembmi_eligible==0)  //Takes value 1 if the household had no eligible females for hemoglobin measurements//
lab var no_fembmi_eligible "Household has no eligible women for bmi measurements"

*** Eligible Man for MR Recode ***
gen	male_eligible = ((edad>=12 & edad<=59) & sexo==1)
replace male_eligible =. if edad==. | sexo==.
bys	hh_id: egen hh_n_male_eligible = sum(male_eligible)  //Number of eligible men for interview in the hh//
gen	no_male_eligible = (hh_n_male_eligible==0) //Takes value 1 if the household had no eligible males for an interview//
lab var no_male_eligible "Household has no eligible man"

*** Eligible Man for anthropometric measurements ** 
*As we have information for men beyond 59 years old, we augment the eligible age to 70 
gen	malebmi_eligible = ((edad>=12 & edad<=70) & sexo==1)
replace malebmi_eligible = . if edad==. | sexo==.  
bys hh_id: egen hh_n_malebmi_eligible = sum(malebmi_eligible) //Number of eligible women for hemoglobin measurements//
gen	no_malebmi_eligible = (hh_n_malebmi_eligible==0) //Takes value 1 if the household had no eligible females for hemoglobin measurements//
lab var no_malebmi_eligible "Household has no eligible man for bmi measurements"

*** Eligible Children anthropometrics *** 
gen	child_eligible = ps79!=. 
replace child_eligible =. if edad>=5
bys	hh_id: egen hh_n_children_eligible = sum(child_eligible)  //Number of eligible children for anthropometrics//
gen	no_child_eligible = (hh_n_children_eligible==0) //Takes value 1 if there were no eligible children for anthropometrics//
lab var no_child_eligible "Household has no children eligible"

*** No eligible Children AND no eligible Women (for athropometrics)***
gen	no_child_fem_eligible = (no_child_eligible==1 & no_fembmi_eligible==1 & no_malebmi_eligible==1)
lab var no_child_fem_eligible "Household has no children, women or male eligible for anthropometric measure"

*** No eligible Women and no eligible Man ***
gen	no_eligibles = (no_fem_eligible==1 & no_male_eligible==1) //Takes value 1 if the household had no eligible members for an interview//
lab var no_eligibles "Household has no eligible man or women"

drop   hh_n_fem_eligible  hh_n_male_eligible child_eligible hh_n_children_eligible
sort hh_id
by   hh_id: gen id = _n

* people
count if no_eligible==1 
* households
count if no_eligible==1  & id==1

* people
count if no_fem_eligible==1 
* households
count if no_fem_eligible==1  & id==1

* people
count if no_male_eligible==1 
* households
count if no_male_eligible ==1  & id==1

* people
count if no_child_eligible==1 
* households
count if no_child_eligible==1  & id==1


*****************************************************************************************
***  Data preparation Step 2							     ***
***  Standardization of the 10 indicators and identification of deprived individuals  ***
*****************************************************************************************
drop weight
*** Rename basic variables 
gen weight=fexp //Sample weight//

//Type of place of residency: urban/rural//
generat urban=1 if  area_5000==1
replace urban=0 if  area_5000==2 
label define lab_urban 1 "urban" 0 "rural"
label values urban lab_urban

/*NOTE: The survey offers two variables for area: 
1) area5000=areas with 5000 people or more are urban and areas with less than 5000 are considered rural
2) area2000=areas with 2000 people or more are urban and areas with less than 2000 are considered rural
we choose area5000.
*/

** Quality check
sum urban [w = weight] 
mat qcheck [1,1] = r(mean)

*** Relationship to the head of household
gen relationship=pd04  
label define relationship ///
           1 "head" ///
           2 "wife or husband" ///
           3 "son/daughter" ///
           4 "son/doughter in law" ///
           5 "grandchild" ///
           6 "parents" ///
           7 "parents in law" ///
           8 "brother/sister" ///
           9 "brother/sister in law" ///
          10 "other relative" ///
          11 "maid" ///
          12 "pensioners" ///
          13 "other no relatives/not related" 
label value relationship relationship

*** Sex household member
gen sex= sexo
label define sex 1 "male" 2"female"
label value sex sex

*** Age household member
gen age=edad
replace age = . if age>=98
label define age 98 "98 or more years"
label value age age

*Total number of de jure hh members in the household
gen usual_resid=1 if relationship<12
replace usual_resid=0 if usual_resid==.
bys hh_id: egen hhs=sum(usual_resid) 
/*NOTE: usual residents is understood as de jure population.
all members, except other no relatives, are considered usual residents
Pensioners and other no relatives are'nt considered as usual resident*/

*region
*label drop region
lab define region 1 "Mountains" 2 "Coast" 3 "Amazon" 4"Galapagos Island"
lab values region region
lab var  region "Region for subnational decompositions"

*******************************
*** Step 2.1 Years of Schooling
*******************************
//The entire household is considered deprived if no household member has completed five years of schooling//

    clonevar nivinst = pe47
	clonevar anoinst = pe48
	
 *** total number of years of education//
	gen eduyears = .
	replace eduyears = 0              if nivinst==1 | nivinst==2 | nivinst==3 | nivinst==4
	replace eduyears = anoinst-1      if nivinst==5 /*Substract 1 because they r counting 1 year for kindergarden*/
	replace eduyears = anoinst        if nivinst==6 /*Here is not need substract 1, not take into account kindergarden*/
	replace eduyears = anoinst+9      if nivinst==7
	replace eduyears = anoinst+6      if nivinst==8
	replace eduyears = anoinst+12     if nivinst>=9 & nivinst<=10
	replace eduyears = anoinst+17     if nivinst==11
	replace eduyears =. if (nivinst==. | anoinst==. | anoinst==99) & nivinst!=1
	replace eduyears =0 if eduyears==-1   
	label var eduyears "years of education approved" 
	
	replace eduyears = . if eduyears>=age & age>0
	replace eduyears = 0 if age<10

//Following a control variable is created on whether there is information on years of 
//education for at least 2/3 of the household members - this is the no_missing_edu variable//
gen	temp = (eduyears~=.)
bys	hh_id: egen no_missing_edu = sum(temp)
replace no_missing_edu =  no_missing_edu/hhs   //already considers only usual residents
replace no_missing_edu = (no_missing_edu>=2/3)
drop temp

* Standard MPI indicator
gen	years_edu5 = (eduyears>=6) // note that I have only updated to 2018 specification
replace years_edu5 = . if eduyears==.

bys	hh_id: egen hh_years_edu5_1 = max(years_edu5)
gen	hh_years_edu5 = (hh_years_edu5_1==1)
replace hh_years_edu5 = . if hh_years_edu5_1==.
replace hh_years_edu5 = . if hh_years_edu5==0 & no_missing_edu==0 //Final var missing if hh has info for < 2/3 of members & for the ones it has info, it is less than 5 years//
replace hh_years_edu5 = . if relationship==13 //we exclude non-usual residents//
lab var hh_years_edu5 "Household has at least one member with 5 years of edu"
// The indicator =1 if at least someone in the hh has reported 5 years of edu or more, and =0 if for those hh for which at least 2/3 of the members reported years, 
// no one has 5 years or more. The indicator has a missing value when there is missing info on years of edu for 2/3 or more of the hh members with no one reporting 5 years of 
// education, or when all household members have missing information on years of education. 

*** Ultrapoverty indicator: nobody completed at least 1 year of schooling
gen	years_edu1 = (eduyears>=1)
replace years_edu1 = . if eduyears==.

bys	hh_id: egen hh_years_edu1_1 = max(years_edu1)
gen	hh_years_edu1 = (hh_years_edu1_1==1)
replace hh_years_edu1 = . if hh_years_edu1_1==.
replace hh_years_edu1 = . if hh_years_edu1==0 & no_missing_edu==0
replace hh_years_edu1 = . if relationship==13 //we exclude non-usual residents//
lab var hh_years_edu1 "Household has at least one member with 1 year of edu"

drop hh_years_edu5_1
drop hh_years_edu1_1


** Quality check
gen	years_edu_check = (eduyears>0)
replace years_edu_check = . if eduyears==.

 
sum hh_years_edu5 [w = weight] if usual_resid==1 
mat qcheck [1,4] = 1 - r(mean)
sum years_edu_check [w = weight] if usual_resid==1 & age>4 & age<.
mat qcheck [1,5] = 1 - r(mean)
sum eduyears [w = weight] if usual_resid==1 & age>4 & age<., detail
mat qcheck [1,7] = r(p50)


************************************
*** Step 2.2 Child School Attendance
************************************
//The entire household is considered deprived if any school-aged child is not attending 
//school up to class 8. Data Source for age children start school: United Nations Educational,
//Scientific and Cultural Organization, Institute for Statistics database, Table 1. Education 
//systems [UIS, http://stats.uis.unesco.org/unesco/TableViewer/tableView.aspx?ReportId=163 ].//

//Please check that it corresponds 1=attending, 0=not attending//

gen attendance =0 if pe18!=.
replace attendance =1 if pe18>=1 & pe18<=7 
replace attendance=0 if edad<=4 /*I put no attendance if the person is younger than 5 years*/
label var  attendance "enrrolled in the current scholar year"
label define attendance 1"yes" 0"no"
label value attendance attendance
* Note: It is no posible estimate asistance instead enrrolled in the current scholar year

//IMPORTANT: Please change the age range according to the compulsory schooling age of the 
//particular country. You will need to check this in 
//http://stats.uis.unesco.org/unesco/TableViewer/tableView.aspx?ReportId=163. 
//Look at the starting age and add 8 (regardless of the finalizing compulsory age in the 
//country. So for example, if the compulsory age is 5-11 years, your age range will be 5-13(=5+8) years. 
//If the question on school attendance refers to the PREVIOUS year 
//as it is the case with hv125), please add one more year to the lower & upper bounds in 
//the definition of the child_schoolage var below. If the variable refers to whether the 
//member is CURRENTLY attending school, as it is the case with hv110 and 121 then use the range itself//

* Note: the compulsary age in Ecuador start at 5 years

* Standard MPI indicator
gen	child_schoolage = (age>=5 & age<=13)
replace child_schoolage=. if age==.
bys	hh_id: egen hh_children_schoolage = sum(child_schoolage)
replace hh_children_schoolage = (hh_children_schoolage>0) //Control variable: It takes value 1 if the household has children in school age//
lab var hh_children_schoolage "Household has children in school age"

gen	child_not_atten = (attendance==0) if child_schoolage==1
replace child_not_atten = . if attendance==. & child_schoolage==1
bys	hh_id: egen any_child_not_atten = max(child_not_atten)
gen	hh_all_child_atten = (any_child_not_atten==0) 
replace hh_all_child_atten = . if any_child_not_atten==.
replace hh_all_child_atten = 1 if hh_children_schoolage==0
lab var hh_all_child_atten    "Household has all school age children in school"

//The indicator takes value 1 if ALL children in school age are attending school and 
//0 if there is at least one child not attending. Households with no children receive a 
//value of 1 as non-deprived. The indicator has a missing value only when there are all 
//missing values on children attendance in households that have children in school age. 
//Please create a missing variable only if hv121, hv125 or hv110 are all missing//
*gen hh_all_child_atten=.

* Ultrapoverty indicator
gen	child_schoolage_6 = (age>=5 & age<=11) 
bys	hh_id: egen hh_children_schoolage_6 = sum(child_schoolage_6)
replace hh_children_schoolage_6 = (hh_children_schoolage_6>0) 
lab var hh_children_schoolage_6 "Household has children in school age (6 years of school)"

gen	child_atten_6 = (attendance==1) if child_schoolage_6==1
replace child_atten_6 = . if attendance==. & child_schoolage_6==1
bys	hh_id: egen any_child_atten_6 = max(child_atten_6)
gen	hh_all_child_atten_6 = (any_child_atten_6==1) 
replace hh_all_child_atten_6 = . if any_child_atten_6==.
replace hh_all_child_atten_6 = 1 if hh_children_schoolage_6==0
lab var hh_all_child_atten_6 "Household has at least one school age children (6 years of school) in school"


** Quality check
sum hh_all_child_atten [w = weight] if usual_resid==1
mat qcheck [1,11] = 1 - r(mean)
sum attendance [w = weight] if child_schoolage==1 & usual_resid==1
mat qcheck [1,12] = r(mean)


**********************
*** Step 2.3 Nutrition
**********************

*** Low BMI of mother or daughter 
gen meter2=(ps84/100)^2
gen bmi =ps82/meter2

gen	f_bmi = bmi if fembmi_eligible==1
lab var f_bmi "Women's BMI "
*NOTE: ECV-2006 does not ask about weight and height of adults. ECV 2013 does ask for adults.

* Standard MPI indicator
gen	f_low_bmi = (f_bmi<18.5)
replace f_low_bmi = . if f_bmi==. | f_bmi>=99.97
lab var f_low_bmi "BMI of mother/daughter < 18.5"

//Takes value 1 if no female in the household has low bmi//
bys	hh_id: egen low_bmi = max(f_low_bmi)
gen	hh_no_low_bmi = (low_bmi==0) 
replace hh_no_low_bmi = . if low_bmi==.
replace hh_no_low_bmi = 1 if no_fembmi_eligible==1 
drop	low_bmi


*** Low BMI of male 
gen m_bmi = bmi if malebmi_eligible==1
lab var m_bmi "Male's BMI "
*drop male_eligible fem_eligible

gen m_low_bmi = (m_bmi<18.5)
replace m_low_bmi = . if m_bmi==. | m_bmi>=99.97
lab var m_low_bmi "BMI of male < 18.5"
*NOTE: ECV-2006 does not ask about weight and height of adults. ECV 2013 does ask for adults.

//Replace hh_no_low_bmi with 0 if there's any male with low bmi, and 1 if not any male has low BMI and information was missing for women//
bys hh_id: egen low_bmi = max(m_low_bmi)
replace hh_no_low_bmi = 0 if low_bmi==1

//hh_no_low_bmi = 0 means indvs fr hhs that have at least 1 member who has low bmi
replace hh_no_low_bmi = 1 if hh_no_low_bmi==. & low_bmi==0
lab var hh_no_low_bmi "Household has no adult with low BMI"
drop low_bmi
 
//NOTE that hh_no_low_bmi takes value 1 if: (a) no any eligible adult in the hh has (observed) low BMI or 
//(b) there are no eligible adults in the hh (One has to check and adjust the dofile so all people who are eligible 
//and/or measured are included - it is particularly important to check if male are measured and what age group among 
//males and females). The variable takes values 0 for those households that have at least one adult with observed low BMI. 
//The variable has a missing value only when there is missing info on BMI for ALL eligible adults in the household//


* Ultrapoverty indicator: threshold of 17 instead of 18.5 for adults 
gen	f_low_bmi_17 = (f_bmi<17)
replace f_low_bmi_17 = . if f_bmi==. | f_bmi>=99.97
lab var f_low_bmi_17 "BMI of mother/daughter <17"

//Takes value 1 if no female in the household has low bmi//
bys	hh_id: egen low_bmi = max(f_low_bmi_17)
gen	hh_no_low_bmi_17 = (low_bmi==0) 
replace hh_no_low_bmi_17 = . if low_bmi==.
replace hh_no_low_bmi_17 = 1 if no_fembmi_eligible==1 
drop	low_bmi

//If you have hb40 (male BMI) please activate the following lines so you also include the information about male//

*** Low BMI of male 
gen	m_low_bmi_17 = (m_bmi<17)
replace m_low_bmi_17 = . if m_bmi==. | m_bmi>=99.97
lab var m_low_bmi_17 "BMI of male <17"
 
//Replace hh_no_low_bmi with 0 if there's any male with low bmi, and 1 if not any male has low BMI and information was missing for women//
bys	hh_id: egen low_bmi = max(m_low_bmi_17)
replace	hh_no_low_bmi_17 = 0 if low_bmi==1
replace hh_no_low_bmi_17 = 1 if hh_no_low_bmi_17==. & low_bmi==0
lab var hh_no_low_bmi_17 "Household has no adult with low BMI (<17)"

//NOTE: Not all surveys have adult BMI, if it does not this variable takes missing value//

** Quality check
sum hh_no_low_bmi [w = weight] if usual_resid==1
mat qcheck [1,20] = 1 - r(mean)
sum f_low_bmi [w = weight] if sex==2 & usual_resid==1
mat qcheck [1,24] = r(mean)
sum m_low_bmi [w = weight] if sex==1 & usual_resid==1
mat qcheck [1,26] = r(mean)

**** Household Child Undernutrition Dummy***
//Households with no eligible children to be measured will receive a value of 1//

* Standard MPI indicator
bys	hh_id: egen temp = max(underwa2)
gen	hh_no_underwa2 = (temp==0) //Takes value 1 if no child in the hh is underweighted, 0 if at least one is//
replace hh_no_underwa2 = . if temp==.
replace hh_no_underwa2 = 1 if no_child_eligible==1 
lab var hh_no_underwa2 "Household has no child under weight-for-age - 2 stdev"
drop	temp

/* NOTE that hh_no_underwh2 takes value 1 if: (a) no any eligible children in the
 hh is undernourished or (b) there are no eligible children in the hh. The variable 
 takes values 0 for those households that have at least one measured child undernourished. 
 The variable has missing values only when there is missing info in nutrition for ALL 
 eligible children in the household*/


* Ultrapoverty indicator
bys	hh_id: egen temp = max(underwa3)
gen	hh_no_underwa3 = (temp==0) //Takes value 1 if no child in the hh is underweighted, 0 if at least one is//
replace hh_no_underwa3 = . if temp==.
replace hh_no_underwa3 = 1 if no_child_eligible==1 
lab var hh_no_underwa3 "Household has no child under weight-for-age - 3 stdev"
drop	temp

** Quality check
count if child_KR==1 & usual_resid==1
mat qcheck [1,29] = r(N)
sum hh_no_underwa2 [w = weight] if usual_resid==1
mat qcheck [1,30] = 1 - r(mean)
sum underwa2 [w = weight] if child_KR==1 & usual_resid==1
mat qcheck [1,33] = r(mean)

*** Finally we create the nutrition indicator to be used in the MPI: hh_nutrition ***
//The indicator takes value 1 if there is no undernourished adult or children. 
//It also takes value 1 for the households that had no eligible adult AND no eligible children. 
//The indicator takes value 0 if any adult or child for whom there is nutritional information is 
//undernourished in the household. The indicator takes value missing "." only if all eligible 
//adults and eligible children have missing information in their respective nutrition variable. 
//If the nutritional variable is missing altogether in the dataset, this indicator will not be 
//included in the MPI and the mortality indicator will receive the full health weight//


* Standard MPI indicator
gen	hh_nutrition = 1
replace hh_nutrition = 0 if hh_no_low_bmi==0 | hh_no_underwa2==0
replace hh_nutrition = . if hh_no_low_bmi==. & hh_no_underwa2==.
replace hh_nutrition = 1 if no_child_fem_eligible==1 
lab var hh_nutrition "Household has no women or child undernourished (weight-for-age)"

* Ultrapoverty indicator
gen	hh_nutrition_3_17 = 1
replace hh_nutrition_3_17 = 0 if hh_no_low_bmi_17==0 | hh_no_underwa3==0
replace hh_nutrition_3_17 = . if hh_no_low_bmi_17==. & hh_no_underwa3==.
replace hh_nutrition_3_17 = 1 if no_child_fem_eligible==1 
lab var hh_nutrition_3_17 "Household has no women or child undernourished (weight-for-age)"

//If there is no nutritional information, then we create an empty variable//
*gen hh_nutrition = .

** Quality check
sum hh_nutrition [w = weight] if usual_resid==1
mat qcheck [1,35] = r(mean)



****************************
*** Step 2.4 Child Mortality
****************************
//The entire household is considered deprived if any child has died in the family. 
//Mortality at any age was considered to make the estimations comparable to those which use 
//MICS datasets. Alternative estimations were performed considering whether there had been 
//death occurrence in children under five years of age only and this does not alter the results 
//significantly. For further details, see: Alkire and Santos (2010).//

//If there is not a male recode create the following empty variables//
gen mv201 = . 
gen mv206 = . 
gen mv207 = . 

//If there is not info on child mortality in the women recode create the following 
//empty variables//
gen v201 = . 
gen v206 = . 
gen v207 = . 

* Standard MPI indicator 
gen child_mortality =pf15-pf17 /*NOTE: the var child_mortality is the number of sons and daughters from women of 12-49 years that have died*/
replace child_mortality=0 if pf05==0 | pf15==0 
lab var child_mortality "Occurrence of child mortality in the household"
bys	hh_id: egen temp = max(child_mortality)
gen	hh_no_child_mortality = (temp==0) 
replace hh_no_child_mortality = . if temp==.
replace hh_no_child_mortality = 1 if no_fem_eligible==1 
lab var hh_no_child_mortality "Household had no child mortality"
drop	temp
//The final indicator takes value 1 if the household was free of child mortality and 0 if at least one children died - according to the individual recode - females. It is missing if there was missing info on both sons and daughters (v206 and v207). Households with no children/no women being interviewed receive a value of 1//


* Ultrapoverty indicator: at least 2 children or more died in the hh 
bys	hh_id: egen child_mortality_f = sum(child_mortality), missing

gen	hh_no_child_mortality_2 = (child_mortality_f<2)
replace hh_no_child_mortality_2 = . if child_mortality_f==.
replace hh_no_child_mortality_2 = 1 if no_fem_eligible==1 
lab var hh_no_child_mortality_2 "Household has less than two children died"
drop child_mortality_f 


** Quality check
sum  hhs [weight=weight] if id==1 
gen  mean_hhs = r(mean)
egen tot_hh = sum(id==1)
gen  i=1
egen tot_pop = sum(i)
drop i

sum hh_no_child_mortality [w = weight] if usual_resid==1
mat qcheck [1,36] = 1 - r(mean)
sum tot_hh [w = weight] if usual_resid==1
mat qcheck [1,37] = r(mean)
sum tot_pop [w = weight] if usual_resid==1
mat qcheck [1,38] = r(mean)
sum mean_hhs [w = weight] if usual_resid==1
mat qcheck [1,39] = r(mean)

************************
*** Step 2.5 Electricity
************************
//Members of the household are considered deprived if the household has no electricity//

* Standard MPI indicator = Ultrapoverty indicator
gen electricity= 0 if vi26<6
replace electricity=1 if (vi26==1 | vi26==2| vi26==3) & electricity==0
label define electricity 0 "no" 1 "yes" //Deprived if no electricity//
label var  electricity "Electricity"
label values electricity lab_yes_no
gen electricity_ultra = electricity

** Quality check
sum electricity [w = weight] if usual_resid==1
mat qcheck [1,44] = 1 - r(mean)

********************************
*** Step 2.6 Improved Sanitation
********************************
//Members of the household are considered deprived if the household's sanitation facility is not 
//improved, according to MDG guidelines, or it is improved but shared with other household. 
//Following the definition of the MDG indicators: "A household is considered to have access to 
//improved sanitation if it uses: Flush or pour flush to piped sewer system, septic tank or pit 
//latrine; Pit latrine with slab; Composting toilet; Ventilated improved pit latrine. The excreta 
//disposal system is considered improved if it is private or shared by a reasonable number of 
//households." Source: "The Challenge of Slums: Global Report on Human Settlements 2003 
//Revised version, April 2010), Chapter 1, p. 16. 
//http://www.unhabitat.org/pmss/listItemDetails.aspx?publicationID=1156"//

gen toilet=vi14
label define toilet ///
1 "toilet connected to sewer" ///
2 "toilet connected to septic well" ///
3 "toilet connected to septic tank" ///
4 "latrine" ///
5 "no toilet facility" 
label value toilet toilet

gen shared_toilet=. //0=no, 1=yes//
replace shared_toilet=1 if vi15a==vi15c
replace shared_toilet=0 if vi15a==vi15b
replace shared_toilet=0 if (vi15a>0 &vi15a!=.) & (vi15b>0 & vi15b!=.)
replace shared_toilet=. if vi15a==.
/*NOTE: If the household has more than one bathroom and has at least one for exclusive use then
the household is considered as not deprived*/

* Standard MPI indicator
gen	toilet_mdg = (toilet<5 & shared_toilet!=1)  
replace toilet_mdg = 0 if toilet<5  & shared_toilet==1 
replace toilet_mdg = . if toilet==.  | toilet==99
lab var toilet_mdg "Household has improved sanitation with MDG Standards"

* Ultrapoverty indicator
gen	toilet_ultra = .
replace toilet_ultra = 0 if toilet==5  //31 open defecation, 96 Other//
replace toilet_ultra = 1 if toilet!=5 & toilet!=. 

** Quality check
sum toilet_mdg [w = weight] if usual_resid==1
mat qcheck [1,48] = 1 - r(mean)
sum shared_toilet [w = weight] if usual_resid==1
mat qcheck [1,49] = r(mean)


********************************
*** Step 2.7 Safe Drinking Water
********************************
//Members of the household are considered deprived if the household does not have access to safe 
//drinking water according to MDG guidelines, or safe drinking water is more than a 30-minute walk from home roundtrip. 
//"A household has improved drinking water supply if it uses water from sources that include: 
//piped water into dwelling, plot or yard; public tap/ stand pipe; tube well/borehole; protected 
//dug well; protected spring; rain water collection. (...) Households using bottled water are 
//only considered to be using improved water when they use water from an improved source for 
//cooking and personal hygiene." Source: "The Challenge of Slums: Global Report on Human 
//Settlements 2003 (Revised version, April 2010), Chapter 1, p. 16, 21. 
//http://www.unhabitat.org/pmss/listItemDetails.aspx?publicationID=1156"
//Access to safe water refers to the percentage of the population with reasonable access to 
//an adequate supply of safe water in their dwelling or within a convenient distance of their 
//dwelling. The Global Water Supply and Sanitation Assessment 2000 Report defines reasonable 
//access as "the availability of 20 litres per capita per day at a distance no longer than 1,000 
//metres". "Indicators for Monitoring the Millennium Development Goals", p. 64-65. 
//[As distance is not available, we convert the 1000 metres distance into 30 minutes. 
//This is the DHS standard too.]//

gen water= vi17
label define water ///
1 "piped water from utility company" ///
2 "other piped water" ///
3 "tanker truck" ///
4 "well" ///
5 "river/stream/spring" ///
6 "other"
label value water water

gen timeminutes= vi21a *60  /*Two household said 10 hours!! */
egen timetowater=rowtotal(timeminutes vi21b), m //Save the original variable//

* Standard MPI indicator
gen	water_mdg = 1 if water==1 | water==2 | water==4 //Non deprived if water is "piped into dwelling", "piped to yard/plot", "public tap/standpipe", "tube well or borehole", "protected well", "protected spring", "rainwater", "bottled water"//
replace water_mdg = 0 if water==3 | water==5 | water==6 //Deprived if it is "unprotected well", "unprotected spring", "surface water (river/lake, etc)", "tanker truck", "cart with small tank", "other"//
replace water_mdg = 0 if water_mdg==1 & timetowater>=30 & timetowater~=.  //Deprived if water is at more than 30 minutes' walk (roundtrip). Please check the value assigned to 'in premises' and if this is different from 996 or 995, add the condition: & timetowater~=XXX accordingly//
lab var water_mdg "Household has drinking water with MDG standards (considering distance)"

* Ultrapoverty indicator
gen	water_mdg_45 = .
replace	water_mdg_45 = 1 if water==1 | water==2 | water==4 //Non deprived if water is "piped into dwelling", "piped to yard/plot", "public tap/standpipe", "tube well or borehole", "protected well", "protected spring", "rainwater", "bottled water"//
replace water_mdg_45 = 0 if water==3 | water==5 | water==6 //Deprived if it is "unprotected well", "unprotected spring", "surface water (river/lake, etc)", "tanker truck", "cart with small tank", "other"//
replace water_mdg_45 = 0 if water_mdg_45==1 & timetowater>45 & timetowater~=. 
lab var water_mdg_45 "Household has drinking water with MDG standards (45 minutes distance)"

** Quality check
gen	water_2 = 1 if water==1 | water==2 | water==4 //Non deprived if water is "piped into dwelling", "piped to yard/plot", "public tap/standpipe", "tube well or borehole", "protected well", "protected spring", "rainwater", "bottled water"//
replace water_2 = 0 if water==3 | water==5 | water==6 //Deprived if it is "unprotected well", "unprotected spring", "surface water (river/lake, etc)", "tanker truck", "cart with small tank", "other"//
replace water_2 = . if water==. | water==99

gen	bottled_water =. /*It is not possible identify bottled of water*/

gen	distance = (timetowater>=30 & timetowater<=.)
replace distance = . if timetowater==.

sum water_2 [w = weight] if usual_resid==1
mat qcheck [1,52] = 1 - r(mean)
sum bottled_water [w = weight] if usual_resid==1
mat qcheck [1,53] = r(mean)
sum water_mdg [w = weight] if usual_resid==1
mat qcheck [1,54] = 1 - r(mean)
sum distance [w = weight] if usual_resid==1
mat qcheck [1,59] = r(mean)



*********************
*** Step 2.8 Flooring
*********************
//Members of the household are considered deprived if the household has a dirt, sand or dung floor//

* Standard MPI indicator = Ultrapoverty indicator
gen floor=vi05 
label define floor 1 "stave / parquet / plank"  ///
2 "ceramic / tile / vinyl"  ///
3 "Marble / marmeton"  ///
4 "cement / brick"  ///
5 "wooden board /  plank untreated"  ///
6 "cane"  ///
7 "sand"  ///
8 "other"
label value floor floor

/*
gen	floor_imp = 1
replace floor_imp = 0 if floor<=6 //Deprived if "earth/sand",  "other"//
replace floor_imp = . if floor==. //Please check that missing values remain missing//
lab var floor_imp "Household has floor that it is not earth/sand/dung"
*/
gen	floor_imp = 1
replace floor_imp = 0 if floor>=7 
replace floor_imp = . if floor==. //Please check that missing values remain missing//
lab var floor_imp "Household has floor that it is not earth/sand/dung"
gen floor_ultra = floor_imp

** Quality check
sum floor_imp [w = weight] if usual_resid==1
mat qcheck [1,61] = 1 - r(mean)


*************************
*** Step 2.9 Cooking Fuel
*************************
//Members of the household are considered deprived if the household cooks with solid fuels: wood, charcoal, crop residues or dung. "Indicators for Monitoring the Millennium Development Goals", p. 63//

* Standard MPI indicator
gen cookingfuel=vi13
label define cookingfuel 1 "gas" 2 "electricity" 3 "wood/charcoal"  4 "other"
label value cookingfuel cookingfuel

gen	cooking_mdg = 1
replace cooking_mdg = 0 if cookingfuel==3 | cookingfuel==4
replace cooking_mdg = . if cookingfuel==. 
lab var cooking_mdg "Househod has cooking fuel according to MDG standards"
//Non deprived if: 1 "electricity", 2 "lpg", 3 "natural gas", 4 "biogas", 5 "kerosene", 96 "other"//
//Deprived     if: 6 "coal/lignite", 7 "charcoal", 8 "wood", 9 "straw/shrubs/grass", 10 "agricultural crop", 11 "animal dung", 95 "no food cooked in household"//

* Ultrapoverty indicator: similar to Cooking Fuel MDG but coal/lignite and charcoal are now not deprived 
gen	cooking_ultra = cooking_mdg
replace cooking_ultra = 0 if cookingfuel==3 | cookingfuel==4 //cannot disentangle wood from charcoal
lab var cooking_ultra "Househod has cooking fuel according to MDG standards (charcoal and coal are not deprived"


** Quality check
sum cooking_mdg [w = weight] if usual_resid==1
mat qcheck [1,65] = 1 - r(mean)

******************************
*** Step 2.10 Assets ownership
******************************
//Members of the household are considered deprived if the household does not own more than one of: radio, TV, telephone, bike, motorbike or refrigerator and does not own a car or truck//
//Check that for all assets in living standards: "no"==0 and yes=="1"//

*Television
*-----------
gen television=.
replace television=1 if  _TV_bn==1 | _TV_color==1 | _TV_plas_lcd==1
replace television=0 if  _TV_bn==2 & _TV_color==2 & _TV_plas_lcd==2
replace television=. if  _TV_bn==. & _TV_color==. & _TV_plas_lcd==.
*Note: if the household has a black/white or color TV it is considered not deprived.

*Radio
*-----------
gen radio=.
replace radio=1 if _grabadora==1 | _eqpsonido==1
replace radio=0 if _grabadora==2 & _eqpsonido==2
replace radio=. if _grabadora==. & _eqpsonido==.
*Note: if the household has a radio (radio equipment or sound equipment) it is considered not deprived.

*Fix Telephone at home
*-----------
gen telephone= 1 if  _telfijo==1
replace telephone= 0 if  _telfijo==2

*Mobile phone
*------------
gen  mobiletelephone =   1 if ph09a==1
replace mobiletelephone= 0 if ph09a==2

*Refrigeradora
*-------------
gen  refrigerator     = 1 if  _Refrigeradora==1
replace  refrigerator = 0 if  _Refrigeradora==2

*Car
*----
gen  car = 1 if     _carro==1
replace  car = 0 if _carro==2

*Bicycle
*-----------
gen bicycle     = 1 if  _Bicicleta==1
replace bicycle = 0 if  _Bicicleta==2

*Motorcycle
*-----------
gen  motorbike   = 1 if _moto==1
replace motorbike= 0 if _moto==2 

*Check not missing values or othe codes than 1 to 2 (excep cellular)
*tab1 TV_blan_negr TV_color Radio_grab Equipo_sonido Lín_telefón celula_1 Refrigeradora Carro Bicicleta Motocicleta
*codebook TV_blan_negr TV_color Radio_grab Equipo_sonido Lín_telefón celula_1 Refrigeradora Carro Bicicleta Motocicleta

gen bw_television=.
replace bw_television=1 if _TV_bn==1 
replace bw_television=0 if _TV_bn==2 
replace bw_television=. if _TV_bn==. 


** Quality check
sum radio [w = weight] if usual_resid==1 & id==1
mat qcheck [1,69] = r(mean)
sum telephone [w = weight] if usual_resid==1 & id==1
mat qcheck [1,72] = r(mean)
sum mobiletelephone [w = weight] if usual_resid==1 & id==1
mat qcheck [1,75] = r(mean)
sum television [w = weight] if usual_resid==1 & id==1
mat qcheck [1,78] = r(mean)
sum bw_television [w = weight] if usual_resid==1 & id==1 // it is already included in television
mat qcheck [1,81] = r(mean)
sum refrigerator [w = weight] if usual_resid==1 & id==1
mat qcheck [1,84] = r(mean)
sum bicycle [w = weight] if usual_resid==1 & id==1
mat qcheck [1,87] = r(mean)
sum motorbike [w = weight] if usual_resid==1 & id==1
mat qcheck [1,90] = r(mean)
sum car [w = weight] if usual_resid==1 & id==1 
mat qcheck [1,93] = r(mean)

//Skip the following lines if black and white tv or mobile phone were missing//
replace telephone=1 if telephone==0 & mobiletelephone==1
replace telephone=1 if telephone==. & mobiletelephone==1

*** Combined Assets Indicator
egen	n_small_assets = rowtotal(television radio telephone refrigerator bicycle motorbike), missing
lab var n_small_assets "Household Number of Small Assets Owned"

* Standard MPI indicator
gen	hh_assets = (car==1 | n_small_assets>1) 
replace hh_assets = . if car==. & n_small_assets==.
lab var hh_assets "Household Asset Ownership: HH has car or more than 1 of small assets"

* Ultrapoverty indicator: only "No Assets" is deprived 
gen	hh_assets_ultra = (car==1 | n_small_assets>0)
replace hh_assets_ultra = . if car==. & n_small_assets==.
lab var hh_assets_ultra "Household Asset Ownership: HH has car or 1 small assets"


** Quality check
sum hh_assets [w = weight] if usual_resid==1 & id==1
mat qcheck [1,96] = r(mean)

***************************************************
*** List with the 10 indicators included in the MPI
***************************************************
local varlist_pov hh_years_edu5 hh_all_child_atten hh_no_child_mortality hh_nutrition electricity toilet_mdg water_mdg floor_imp cooking_mdg hh_assets

************************************************************************
*** Define deprivation vector `d_' 
*** which takes values 1 if individual is deprived in the particular 
*** indicator according to deprivation cutoff z as defined during step 2 
************************************************************************
gen d_cm=(hh_no_child_mortality==0) if hh_no_child_mortality!=.
gen d_nutr=(hh_nutrition==0) if hh_nutrition!=.
gen d_satt=(hh_all_child_atten==0) if hh_all_child_atten!=.
gen d_educ=(hh_years_edu5==0) if hh_years_edu5!=.
gen d_wtr=(water_mdg==0) if water_mdg!=.
gen d_elct=(electricity==0) if electricity!=.
gen d_sani=(toilet_mdg==0) if toilet_mdg!=.
gen d_ckfl=(cooking_mdg==0) if cooking_mdg!=.
gen d_hsg=(floor_imp==0) if floor_imp!=.
gen d_asst=(hh_assets==0) if hh_assets!=.

**** INDIVIDUAL DEPRIVATIONS AND POVERTY STATUS
* We will assume equal dimensional weights
* and equal wights within each dimension
	global ind cm nutr satt educ elct wtr sani hsg ckfl asst
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

gen psu=identif_sect
gen strata=zdp
gen sweight=fexp
	
keep psu strata sweight hh_id identif_hog provincia region ciudad zona sector vivienda hogar d_*
cd "C:\Users\qehs1094\Dropbox\Endogenous Weights\Datasets and dofiles"
save mpi_ECU_full, replace
