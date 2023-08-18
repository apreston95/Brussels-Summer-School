clear all

import excel "/Users/andrewpreston/Library/Mobile Documents/com~apple~CloudDocs/Brussels Summer School/Practice session I/PSI1/recession_prob.xls", sheet("Sheet1") firstrow

g k = _n
tsset k 

g b_recession_us = .

g t_recession_us = .

g recession_us_k_3 = (l3.recession_us==1) | (l2.recession_us==1) | (l.recession_us==1) | (recession_us==1) | (f.recession_us==1) | (f2.recession_us==1) | (f3.recession_us==1) 


forvalues i = 12/90 {
reg f`i'.recession_us_k_3 recession_us
replace b_recession_us = _b[_cons] + _b[recession_us] if _n == `i'
replace t_recession_us =  _b[recession_us]/_se[recession_us] if _n == `i'
}

label variable b_recession_us "Recession Probability"
label variable t_recession_us "t-statistic for 	β_(recession)"

tsline b_recession_us if inrange(k,12,90), lw(medthick) lc(red) title("US")

tsline t_recession_us if inrange(k,12,90), lw(medthick) lc(navy) title("US") yline(1.96) yline(-1.96) saving(t_US.gph, replace)

g b_recession_uk = .

g t_recession_uk = .

g recession_uk_k_3 = (l3.recession_uk==1) | (l2.recession_uk==1) | (l.recession_uk==1) | (recession_uk==1) | (f.recession_uk==1) | (f2.recession_uk==1) | (f3.recession_uk==1) 


forvalues i = 12/90 {
reg f`i'.recession_uk_k_3 recession_uk
replace b_recession_uk = _b[_cons] + _b[recession_uk] if _n == `i'
replace t_recession_uk =  _b[recession_uk]/_se[recession_uk] if _n == `i'
}

label variable b_recession_uk "Recession Probability"
label variable t_recession_uk "t-statistic for β_(recession)"

tsline b_recession_uk if inrange(k,12,90), lw(medthick) lc(red) title("UK")

tsline t_recession_uk if inrange(k,12,90), lw(medthick) lc(navy) title("UK") yline(1.96) yline(-1.96) saving(t_UK.gph, replace)



g b_recession_ca = .

g t_recession_ca = .

g recession_ca_k_3 = (l3.recession_ca==1) | (l2.recession_ca==1) | (l.recession_ca==1) | (recession_ca==1) | (f.recession_ca==1) | (f2.recession_ca==1) | (f3.recession_ca==1) 


forvalues i = 12/90 {
reg f`i'.recession_ca_k_3 recession_ca
replace b_recession_ca = _b[_cons] + _b[recession_ca] if _n == `i'
replace t_recession_ca =  _b[recession_ca]/_se[recession_ca] if _n == `i'
}

label variable b_recession_ca "Recession Probability"
label variable t_recession_ca "t-statistic for β_(recession)"

tsline b_recession_ca if inrange(k,12,90), lw(medthick) lc(red) title("Canada")

tsline t_recession_ca if inrange(k,12,90), lw(medthick) lc(navy) title("Canada") yline(1.96) yline(-1.96) saving(t_CA.gph, replace)

graph combine t_US.gph t_UK.gph t_CA.gph, rows(1)

graph export t_stats.pdf


/// Could also use a logistic regression to estimate the regression

forvalues i = 12/90 {
logistic f`i'.recession_us_k_3 recession_us

}
