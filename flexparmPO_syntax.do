
use data_sample, clear

count
global nobs= r(N)

global nknots = 3 //Number of spline knots including boundary knots
global maxtime = 5 //Administrative censoring time
global integration_precision = 1000 //Number of points in which to evaluate integration (in likelihood function)



capture program drop find_spline_knots
program define find_spline_knots, rclass
	version 16.1
	args t1 d1 t2 d2

	tempvar x1 x2
	generate double `x1' = ln(`t1')
	generate double `x2' = ln(`t2')

	// Spline for transition 01
	{
		quietly: summarize `x1' if `d1'==1, meanonly
		local lowerknot `r(min)'
		local upperknot `r(max)'
		local splineknots "`lowerknot'"
		local percentilelist ""
		forvalues k = 1/`=$nknots-2' {
			local percentilelist = "`percentilelist'" +  " `=100*(`k'/($nknots-1))'"
		}
		quietly: _pctile `x1' if `d1'==1, p(`percentilelist')
		forvalues k = 1/`=$nknots-2' {
			local splineknots = "`splineknots'" + " `" + "r(r" + "`k'" + ")" + "'"
		}
	}
	local splineknots = "`splineknots'" + " `upperknot'"
	local splineknots01 = "`splineknots'"
	local drop splineknots

	// Spline for transition 02
	{
		quietly: summarize `x2' if `d1'==0 & `d2'==1, meanonly
		local lowerknot `r(min)'
		local upperknot `r(max)'
		local splineknots "`lowerknot'"
		local percentilelist ""
		forvalues k = 1/`=$nknots-2' {
			local percentilelist = "`percentilelist'" +  " `=100*(`k'/($nknots-1))'"
		}
		quietly: _pctile `x2' if  `d1'==0 & `d2'==1, p(`percentilelist')
		forvalues k = 1/`=$nknots-2' {
			local splineknots = "`splineknots'" + " `" + "r(r" + "`k'" + ")" + "'"
		}
	}
	local splineknots = "`splineknots'" + " `upperknot'"
	local splineknots02 = "`splineknots'"
	local drop splineknots

	// Spline for transition 12
	{
		quietly: summarize `x2' if `d1'==1 & `d2'==1, meanonly
		local lowerknot `r(min)'
		local upperknot `r(max)'
		local splineknots "`lowerknot'"
		local percentilelist ""
		forvalues k = 1/`=$nknots-2' {
			local percentilelist = "`percentilelist'" +  " `=100*(`k'/($nknots-1))'"
		}
		quietly: _pctile `x2' if  `d1'==1 & `d2'==1, p(`percentilelist')
		forvalues k = 1/`=$nknots-2' {
			local splineknots = "`splineknots'" + " `" + "r(r" + "`k'" + ")" + "'"
		}
	}
	local splineknots = "`splineknots'" + " `upperknot'"
	local splineknots12 = "`splineknots'"
	local drop splineknots

	return local splineknots01 = "`splineknots01'"
	return local splineknots02 = "`splineknots02'"
	return local splineknots12 = "`splineknots12'"
end


capture program drop spline_prep
program define spline_prep, nclass
	version 16.1
	args t1 d1 t2 d2
	capture: drop rcs01_t1_* d_rcs01_t1_*
	capture: drop rcs12_t1_* d_rcs12_t1_*
	capture: drop rcs02_t1_* d_rcs02_t1_*
	capture: drop rcs01_t2_* d_rcs01_t2_*
	capture: drop rcs12_t2_* d_rcs12_t2_*
	capture: drop rcs02_t2_* d_rcs02_t2_*

	tempvar x1 x2
	generate double `x1' = ln(`t1')
	generate double `x2' = ln(`t2')

	// Spline for transition 01
	quietly: rcsgen `x1', knots($splineknots01) gen(rcs01_t1_) dgen(d_rcs01_t1_)
	quietly: rcsgen `x2', knots($splineknots01) gen(rcs01_t2_) dgen(d_rcs01_t2_)

	// Spline for transition 02
	quietly: rcsgen `x1', knots($splineknots02) gen(rcs02_t1_) dgen(d_rcs02_t1_)
	quietly: rcsgen `x2', knots($splineknots02) gen(rcs02_t2_) dgen(d_rcs02_t2_)

	// Spline for transition 12
	quietly: rcsgen `x1', knots($splineknots12) gen(rcs12_t1_) dgen(d_rcs12_t1_)
	quietly: rcsgen `x2', knots($splineknots12) gen(rcs12_t2_) dgen(d_rcs12_t2_)
end



capture program drop integration_setup
program define integration_setup, nclass
	version 16.1
	args int_nobs eofu

	quietly: capture: frame create integration
	frame integration {
		quietly: clear
		quietly: set obs `int_nobs'
		quietly: generate u = .
		quietly: replace u = ((_n)/`int_nobs')*`eofu'
		quietly: generate ln_u = ln(u)
		quietly: rcsgen ln_u, knots($splineknots01) gen(rcs01_u_) dgen(d_rcs01_u_)
		quietly: rcsgen ln_u, knots($splineknots12) gen(rcs12_u_) dgen(d_rcs12_u_)
		quietly: rcsgen ln_u, knots($splineknots02) gen(rcs02_u_) dgen(d_rcs02_u_)
		quietly: set obs `=`int_nobs'+1'
		quietly: replace u = 0 in `=`int_nobs'+1'
		sort u
	}
end	




capture program drop ic_initvals_fpmodel
program define ic_initvals_fpmodel, rclass
	version 16.1
	args t1 d1 l1 t2 d2 scenario

	tempvar trans01_t trans01_d trans12_t0 trans12_t trans12_d trans02_t trans02_d

	quietly: generate `trans01_t' = .
	quietly: generate `trans01_d' = .
	quietly: generate `trans02_t' = .
	quietly: generate `trans02_d' = .
	quietly: generate `trans12_t0' = .
	quietly: generate `trans12_t' = .
	quietly: generate `trans12_d' = .

	// Transition 01
	{
		quietly: replace `trans01_t' = `t1' 					if inlist(`scenario',1,4)
		quietly: replace `trans01_d' = 1 						if inlist(`scenario',1,4)
		quietly: replace `trans01_t' = `l1'						if inlist(`scenario',2,5)
		quietly: replace `trans01_d' = 0 						if inlist(`scenario',2,5)
		quietly: replace `trans01_t' = `l1' + ((`t1'-`l1')/2) 	if inlist(`scenario',3,6)
		quietly: replace `trans01_d' = 1 						if inlist(`scenario',3,6)

		quietly: stset `trans01_t', failure(`trans01_d'==1)
		quietly: merlin (_t, family(rp, failure(_d) ltruncated(_t0) df(`=$nknots-1') noorthog))
		matrix initmat01 = e(b)
		forvalues k = 0/`=$nknots-1' {
			local init01_`k' = initmat01[1,`=`k'+1']
			return local init01_`k' `init01_`k''
		}
		quietly: predict survest01, survival
		quietly: capture: frame drop survest
		quietly: frame put _t survest01, into(survest)
		quietly: drop survest01
	}

	// Transition 02
	{
		// We consider some of those who had a negative examination as having had an event of type 1.
		// The proportion is determined by the estimated difference in survival for the 01 transition
		quietly: replace `trans02_t' = `t1' 					if inlist(`scenario',1,4)
		quietly: replace `trans02_d' = 0 						if inlist(`scenario',1,4)
		quietly: replace `trans02_t' = `l1' + ((`t1'-`l1')/2) 	if inlist(`scenario',3,6)
		quietly: replace `trans02_d' = 0 						if inlist(`scenario',3,6)
		quietly: count
		local nobs = r(N)
		forvalues i = 1/`nobs' {
			if `scenario'[`i']==2 | `scenario'[`i']==5 {
				scalar l1i = `l1'[`i']
				scalar t2i = `t2'[`i']
				frame survest {
					quietly: summarize _t, meanonly
					scalar min_t = r(min)
					quietly: summarize survest01 if _t <= scalar(l1i), meanonly
					scalar survest_l1i = r(min)
					quietly: summarize survest01 if _t <= max(scalar(t2i),scalar(min_t)), meanonly
					scalar survest_t2i = r(min)
				}
				if scalar(l1i) == 0 | scalar(l1i) == . {
					scalar survdiff = 1 - scalar(survest_t2i)
				}
				else {
					scalar survdiff = scalar(survest_l1i) - scalar(survest_t2i)
				}
				if scalar(survdiff) == 0 {
					quietly: replace `trans02_d' = `d2' in `i'
					quietly: replace `trans02_t' = `t2' in `i'
				}
				else {
					local temp = rbinomial(1,scalar(survdiff))
					if `temp' == 0 {
						quietly: replace `trans02_t' = `t2' in `i'
						quietly: replace `trans02_d' = `d2' in `i'
					}
					else if `temp' == 1 {
						quietly: replace `trans02_t' = `l1' + ((`t2'-`l1')/2) in `i'
						quietly: replace `trans02_d' = 0 in `i'
						quietly: replace `trans01_d' = 1 in `i'
						quietly: replace `trans01_t' = `l1' + ((`t2'-`l1')/2) in `i'
						quietly: replace `trans12_d' = `d2' in `i'
						quietly: replace `trans12_t0' = `l1' + ((`t2'-`l1')/2) in `i'
						quietly: replace `trans12_t' = `t2' in `i'
					}
				}
			}
		}

		quietly: stset `trans02_t', failure(`trans02_d'==1)
		quietly: merlin (_t, family(rp, failure(_d) ltruncated(_t0) df(`=$nknots-1') noorthog))
		matrix initmat02 = e(b)
		forvalues k = 0/`=$nknots-1' {
			local init02_`k' = initmat02[1,`=`k'+1']
			return local init02_`k' `init02_`k''
		}
	}

	// Transition 12
	{
		quietly: replace `trans12_t0' = `t1' 					if inlist(`scenario',1,4)
		quietly: replace `trans12_t' = `t2' 					if inlist(`scenario',1,4)
		quietly: replace `trans12_d' = `d2' 					if inlist(`scenario',1,4)
		quietly: replace `trans12_t0' = `l1' + ((`t1'-`l1')/2) 	if inlist(`scenario',3,6)
		quietly: replace `trans12_t' = `t2' 					if inlist(`scenario',3,6)
		quietly: replace `trans12_d' = `d2' 					if inlist(`scenario',3,6)

		quietly: stset `trans12_t' if `trans01_d'==1, failure(`trans12_d'==1) entry(time `trans12_t0')
		quietly: merlin (_t, family(rp, failure(_d) ltruncated(_t0) df(`=$nknots-1') noorthog))
		matrix initmat12 = e(b)
		forvalues k = 0/`=$nknots-1' {
			local init12_`k' = initmat12[1,`=`k'+1']
			return local init12_`k' `init12_`k''
		}
	}

	quietly: stset, clear
	quietly: drop _rcs1_*
	quietly: frame drop survest
end




capture program drop estimate_F01
program define estimate_F01, rclass
	args timepoint

	matrix define B = e(b)

	quietly: capture: frame create F01_estimation
	frame F01_estimation {

		clear
		set obs 100000
		generate time = .
		replace time = ((_n)/100000)*`timepoint'
		generate double ln_time = ln(time)
		rcsgen ln_time, knots($splineknots01) gen(rcs01_time_) dgen(d_rcs01_time_)
		rcsgen ln_time, knots($splineknots02) gen(rcs02_time_)

		local eq01 = `"B[1,"gamma01_0:_cons"]"'
		local eq02 = `"B[1,"gamma02_0:_cons"]"'
		forvalues k = 1/`=$nknots-1' {
			local eq01 = `"`eq01'"' + `" + B[1,"gamma01_`k':_cons"]*rcs01_time_`k'"'
			local eq02 = `"`eq02'"' + `" + B[1,"gamma02_`k':_cons"]*rcs02_time_`k'"'
		}

		generate double spline01 = `eq01'
		generate double spline02 = `eq02'

		generate double cumhaz01 = exp(spline01)
		generate double cumhaz02 = exp(spline02)

		local d_eq01 = `"B[1,"gamma01_1:_cons"]*d_rcs01_time_1"'
		forvalues k = 2/`=$nknots-1' {
			local d_eq01 = `"`d_eq01'"' + `" + B[1,"gamma01_`k':_cons"]*d_rcs01_time_`k'"'
		}
		generate double hazard01 = (1/time)*(`d_eq01')*cumhaz01

		generate double integrand = hazard01 * exp(-cumhaz01-cumhaz02)
		integ integrand time

	}
	return scalar F01_estimate = r(integral)
end




capture program drop intcens_lf0
program define intcens_lf0
	version 16.1
	args todo b lnfj

	local gamma01str = "gamma01_0"
	local gamma12str = "gamma12_0"
	local gamma02str = "gamma02_0"
	forvalues k = 1/`=$nknots-1' {
		local gamma01str = "`gamma01str'" + " gamma01_`k'"
		local gamma12str = "`gamma12str'" + " gamma12_`k'"
		local gamma02str = "`gamma02str'" + " gamma02_`k'"
	}
	tempname `gamma01str' `gamma12str' `gamma02str'
	forvalues k = 1/`=$nknots' {
		mleval `gamma01_`=`k'-1'' = `b', eq(`=(`k'-1)*3+1') scalar
		mleval `gamma12_`=`k'-1'' = `b', eq(`=(`k'-1)*3+2') scalar
		mleval `gamma02_`=`k'-1'' = `b', eq(`=(`k'-1)*3+3') scalar
	}
	local t1 "$ML_y1"
	local d1 "$ML_y2"
	local l1 "$ML_y3"
	local t2 "$ML_y4"
	local d2 "$ML_y5"

	/*Functions evaluated at t1*/
	tempvar s01_t1 ddx_s01_t1 s12_t1 s02_t1 H01_t1 H12_t1 H02_t1 h01_t1 S_t1

	local s01_t1_eq = "`" + "gamma01_0" + "'"
	local s12_t1_eq = "`" + "gamma12_0" + "'"
	local s02_t1_eq = "`" + "gamma02_0" + "'"
	forvalues k = 1/`=$nknots-1' {
		local s01_t1_eq = "`s01_t1_eq'" + " + `" + "gamma01_`k'" + "'" + "*rcs01_t1_`k'"
		local s12_t1_eq = "`s12_t1_eq'" + " + `" + "gamma12_`k'" + "'" + "*rcs12_t1_`k'"
		local s02_t1_eq = "`s02_t1_eq'" + " + `" + "gamma02_`k'" + "'" + "*rcs02_t1_`k'"
	}
	local ddx_s01_t1_eq = "`" + "gamma01_1" + "'" + "*d_rcs01_t1_1"
	forvalues k = 2/`=$nknots-1' {
		local ddx_s01_t1_eq = "`ddx_s01_t1_eq'" + " + `" + "gamma01_`k'" + "'" + "*d_rcs01_t1_`k'"
	}

	quietly: generate double `s01_t1' = `s01_t1_eq'
	quietly: generate double `ddx_s01_t1' = `ddx_s01_t1_eq'
	quietly: generate double `s12_t1' = `s12_t1_eq'
	quietly: generate double `s02_t1' = `s02_t1_eq'

	quietly: generate double `H01_t1' = exp(`s01_t1')
	quietly: generate double `H12_t1' = exp(`s12_t1')
	quietly: generate double `H02_t1' = exp(`s02_t1')
	quietly: generate double `h01_t1' = (1/`t1')*`ddx_s01_t1'*`H01_t1'

	quietly: generate double `S_t1' = exp(-`H01_t1'-`H02_t1')

	/*Functions evaluated at t2*/
	tempvar s01_t2 s12_t2 ddx_s12_t2 s02_t2 ddx_s02_t2 H01_t2 H02_t2 H12_t2 h12_t2 h02_t2 S_t2

	local s01_t2_eq = "`" + "gamma01_0" + "'"
	local s12_t2_eq = "`" + "gamma12_0" + "'"
	local s02_t2_eq = "`" + "gamma02_0" + "'"
	forvalues k = 1/`=$nknots-1' {
		local s01_t2_eq = "`s01_t2_eq'" + " + `" + "gamma01_`k'" + "'" + "*rcs01_t2_`k'"
		local s12_t2_eq = "`s12_t2_eq'" + " + `" + "gamma12_`k'" + "'" + "*rcs12_t2_`k'"
		local s02_t2_eq = "`s02_t2_eq'" + " + `" + "gamma02_`k'" + "'" + "*rcs02_t2_`k'"
	}
	local ddx_s12_t2_eq = "`" + "gamma12_1" + "'" + "*d_rcs12_t2_1"
	local ddx_s02_t2_eq = "`" + "gamma02_1" + "'" + "*d_rcs02_t2_1"
	forvalues k = 2/`=$nknots-1' {
		local ddx_s12_t2_eq = "`ddx_s12_t2_eq'" + " + `" + "gamma12_`k'" + "'" + "*d_rcs12_t2_`k'"
		local ddx_s02_t2_eq = "`ddx_s02_t2_eq'" + " + `" + "gamma02_`k'" + "'" + "*d_rcs02_t2_`k'"
	}

	quietly: generate double `s01_t2' = `s01_t2_eq'
	quietly: generate double `s12_t2' = `s12_t2_eq'
	quietly: generate double `ddx_s12_t2' = `ddx_s12_t2_eq'
	quietly: generate double `s02_t2' = `s02_t2_eq'
	quietly: generate double `ddx_s02_t2' = `ddx_s02_t2_eq'

	quietly: generate double `H01_t2' = exp(`s01_t2')
	quietly: generate double `H12_t2' = exp(`s12_t2')
	quietly: generate double `H02_t2' = exp(`s02_t2')
	quietly: generate double `h12_t2' = (1/`t2')*`ddx_s12_t2'*`H12_t2'
	quietly: generate double `h02_t2' = (1/`t2')*`ddx_s02_t2'*`H02_t2'

	quietly: generate double `S_t2' = exp(-`H01_t2'-`H02_t2')

	forvalues k = 0/`=$nknots-1' {
		scalar g01_`k' = `gamma01_`k''
		scalar g12_`k' = `gamma12_`k''
		scalar g02_`k' = `gamma02_`k''
	}

	frame integration {
		tempvar s01_u ddx_s01_u s12_u s02_u H01_u H12_u H02_u h01_u
		local s01_u_eq = "scalar(g01_0)"
		local s12_u_eq = "scalar(g12_0)"
		local s02_u_eq = "scalar(g02_0)"
		forvalues k = 1/`=$nknots-1' {
			local s01_u_eq = "`s01_u_eq'" + " + scalar(g01_`k')*rcs01_u_`k'"
			local s12_u_eq = "`s12_u_eq'" + " + scalar(g12_`k')*rcs12_u_`k'"
			local s02_u_eq = "`s02_u_eq'" + " + scalar(g02_`k')*rcs02_u_`k'"
		}
		local ddx_s01_u_eq = "scalar(g01_1)*d_rcs01_u_1"
		forvalues k = 2/`=$nknots-1' {
			local ddx_s01_u_eq = "`ddx_s01_u_eq'" + " + scalar(g01_`k')*d_rcs01_u_`k'"
		}
		quietly: generate double `s01_u' = `s01_u_eq'
		quietly: generate double `ddx_s01_u' = `ddx_s01_u_eq'
		quietly: generate double `s12_u' = `s12_u_eq'
		quietly: generate double `s02_u' = `s02_u_eq'
		quietly: generate double `H01_u' = exp(`s01_u')
		quietly: generate double `H12_u' = exp(`s12_u')
		quietly: generate double `H02_u' = exp(`s02_u')
		quietly: generate double `h01_u' = (1/u)*`ddx_s01_u'*`H01_u'
		capture: drop integrand integral
		quietly: generate double integrand = exp(-`H01_u' - `H02_u') * `h01_u' * (1 / exp(-`H12_u'))
		quietly: replace integrand = 0 if u == 0
		quietly: integ integrand u, generate(integral) double
		quietly: replace integral = 0 if u == 0
	}

	quietly: count
	local nobs = r(N)
	forvalues i = 1/`nobs' {
		capture: scalar drop integral_value
		if scenario[`i']==2 | scenario[`i']==5 {
			scalar l1i = `l1'[`i']
			scalar t2i = `t2'[`i']
			frame integration {
				quietly: summarize integral if u <= scalar(l1i), meanonly
				scalar integral_l1i = r(max)
				quietly: summarize integral if u <= scalar(t2i), meanonly
				scalar integral_t2i = r(max)
			}
			scalar integral_value = scalar(integral_t2i) - scalar(integral_l1i)
		}
		else if scenario[`i']==3 | scenario[`i']==6 {
			scalar l1i = `l1'[`i']
			scalar t1i = `t1'[`i']
			frame integration {
				quietly: summarize integral if u <= scalar(l1i), meanonly
				scalar integral_l1i = r(max)
				quietly: summarize integral if u <= scalar(t1i), meanonly
				scalar integral_t1i = r(max)
			}
			scalar integral_value = scalar(integral_t1i) - scalar(integral_l1i)
		}
		if scenario[`i']==1 {
			quietly: replace `lnfj' = ln(`S_t1' * `h01_t1' * (exp(-`H12_t2')/exp(-`H12_t1'))) in `i'
		}
		else if scenario[`i']==2 {
			quietly: replace `lnfj' = ln(`S_t2' + (exp(-`H12_t2') * scalar(integral_value))) in `i'
		}
		else if scenario[`i']==3 {
			quietly: replace `lnfj' = ln(exp(-`H12_t2') * scalar(integral_value)) in `i'
		}
		else if scenario[`i']==4 {
			quietly: replace `lnfj' = ln(`S_t1' * `h01_t1' * (exp(-`H12_t2')/exp(-`H12_t1')) * `h12_t2') in `i'
		}
		else if scenario[`i']==5 {
			quietly: replace `lnfj' = ln((`S_t2' * `h02_t2') + (exp(-`H12_t2') * scalar(integral_value) * `h12_t2')) in `i'
		}
		else if scenario[`i']==6 {
			quietly: replace `lnfj' = ln(exp(-`H12_t2') * scalar(integral_value) * `h12_t2') in `i'
		}
	}
end



capture program drop estimate_splines
program define estimate_splines, nclass
	args t1 d1 l1 t2 d2 scenario iter
	quietly {
		integration_setup $integration_precision $maxtime
		ic_initvals_fpmodel `t1' `d1' `l1' `t2' `d2' `scenario'
		local init01_string = ""
		local init12_string = ""
		local init02_string = ""
		local eq_string = ""
		forvalues k = 0/`=$nknots-1' {
			local init01_`k' = r(init01_`k')
			local init01_string = "`init01_string'" + " gamma01_`k':_cons=" + "`" + "init01_`k'" + "'"
			local init12_`k' = r(init12_`k')
			local init12_string = "`init12_string'" + " gamma12_`k':_cons=" + "`" + "init12_`k'" + "'"
			local init02_`k' = r(init02_`k')
			local init02_string = "`init02_string'" + " gamma02_`k':_cons=" + "`" + "init02_`k'" + "'"
			local eq_string = "`eq_string'" + " (gamma01_`k': `t1' `d1' `l1' `t2' `d2' = )"
			local eq_string = "`eq_string'" + " (gamma12_`k': `t1' `d1' `l1' `t2' `d2' = )"
			local eq_string = "`eq_string'" + " (gamma02_`k': `t1' `d1' `l1' `t2' `d2' = )"
		}
	}
	quietly {
		ml model lf0 intcens_lf0 ///
			`eq_string' ///
			, init(`init01_string' `init12_string' `init02_string') missing nopreserve
	}
	quietly: ml maximize, difficult iter(`iter') search(norescale)
	quietly: frame drop integration
end


capture program drop estimation_general_setup
program define estimation_general_setup, nclass
	args t1 d1 l1 t2 d2 scenario iter time1 time2

	find_spline_knots `t1' `d1' `t2' `d2'
	global splineknots01 = r(splineknots01)
	global splineknots12 = r(splineknots12)
	global splineknots02 = r(splineknots02)
	spline_prep `t1' `d1' `t2' `d2'
	capture noisily {
		estimate_splines `t1' `d1' `l1' `t2' `d2' `scenario' `iter'
		matrix full_sample_estimates = r(table)
		forvalues time = `time1'/`time2' {
			estimate_F01 `time'
			quietly: generate _F01_`time'_full = r(F01_estimate)
			quietly: generate _F01_`time'_jackknife = .
		}
		forvalues j = 1/$nobs {
			capture noisily {
				capture: frame drop jackknife_frame
				frame put if _n != `j', into(jackknife_frame)
				frame jackknife_frame {
					integration_setup $integration_precision $maxtime
					local init01_string = ""
					local init12_string = ""
					local init02_string = ""
					local eq_string = ""
					forvalues k = 0/`=$nknots-1' {
						local init01_`k' = full_sample_estimates[rownumb(full_sample_estimates,"b"),colnumb(full_sample_estimates,"gamma01_`k':_cons")]
						local init01_string = "`init01_string'" + " gamma01_`k':_cons=" + "`" + "init01_`k'" + "'"
						local init12_`k' = full_sample_estimates[rownumb(full_sample_estimates,"b"),colnumb(full_sample_estimates,"gamma12_`k':_cons")]
						local init12_string = "`init12_string'" + " gamma12_`k':_cons=" + "`" + "init12_`k'" + "'"
						local init02_`k' = full_sample_estimates[rownumb(full_sample_estimates,"b"),colnumb(full_sample_estimates,"gamma02_`k':_cons")]
						local init02_string = "`init02_string'" + " gamma02_`k':_cons=" + "`" + "init02_`k'" + "'"
						local eq_string = "`eq_string'" + " (gamma01_`k': `t1' `d1' `l1' `t2' `d2' = )"
						local eq_string = "`eq_string'" + " (gamma12_`k': `t1' `d1' `l1' `t2' `d2' = )"
						local eq_string = "`eq_string'" + " (gamma02_`k': `t1' `d1' `l1' `t2' `d2' = )"
					}
					quietly {
						ml model lf0 intcens_lf0 ///
							`eq_string' ///
							, init(`init01_string' `init12_string' `init02_string') missing
					}

					quietly: capture: ml maximize, difficult iter(`iter') search(norescale)
					scalar ml_max_converged = e(converged)
					quietly: capture: frame drop integration
				}
				frame drop jackknife_frame
			}
			if scalar(ml_max_converged) == 1 {
				capture {
					forvalues time = `time1'/`time2' {
						estimate_F01 `time'
						quietly: replace _F01_`time'_jackknife = r(F01_estimate) in `j'
					}
				}
			}
		}
		forvalues time = `time1'/`time2' {
			quietly: generate pseudo_`time' = $nobs*_F01_`time'_full - ($nobs-1)*_F01_`time'_jackknife
			drop _F01_`time'_full _F01_`time'_jackknife
		}
	}

end


timer on 1
estimation_general_setup t1 d1 l1 t2 d2 scenario 25 3 3
timer off 1
timer list 1
