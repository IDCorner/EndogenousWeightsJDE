clear all
mat def M1=(0.3,  0  , 0 \  ///
            0,    0.3, 0.4 \ ///
		    0,    0.4, 0.2)

mat def M2=(0.3,  0  , 0.02 \  ///
            0,    0.3, 0.4 \ ///
		    0.02, 0.4, 0.2)

mat def M3=(0.3,  0.1, 0.02 \  ///
            0.1,  0.3, 0.4 \ ///
		    0.02, 0.4, 0.2)

mat def M4=(0.3,  0.1, 0.2 \  ///
            0.1,  0.3, 0.4 \ ///
		    0.2,  0.4, 0.2)
			
foreach i of numlist 1 2 3 4 {
	mat A`i'=M`i'*M`i'
	mat a`i'=vecdiag(A`i')
	mat B`i'=diag(a`i')
	foreach j of numlist 1 2 3 {
		mat B`i'[`j',`j']=sqrt(B`i'[`j',`j'])
	}
	mat C`i'=inv(B`i')*A`i'*inv(B`i')
	foreach j of numlist 1 2 3 {
		foreach k of numlist 1 2 3 {
			if `k'>`j' {
				mat C`i'[`j',`k']=C`i'[`k',`j']
			}
		}
	}
	pcamat C`i', n(1000) comp(1) // names(r1 r2 r3) 
	mat W_`i'=e(L)
	mat NW_`i'=J(3,1,.)
	foreach h of numlist 1 2 3 {
		mat NW_`i'[`h',1]=(W_`i'[`h',1]^2)*100
	}
}

mat NW=(NW_1,NW_2,NW_3,NW_4)
mat list NW
foreach i of numlist 1 2 3 4 {
	esttab matrix(C`i', fmt(%9.3f)) using "C:\Users\qehs1094\Dropbox\Endogenous Weights\Datasets and dofiles\graphs and tables\C`i'.tex", replace
}
esttab matrix(NW, fmt(%9.2f)) using "C:\Users\qehs1094\Dropbox\Endogenous Weights\Datasets and dofiles\graphs and tables\NW.tex", replace

	