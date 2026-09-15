*###############################################################################
* DO-FILE (Station-Specific Pass-Through)
* Spatial Competition and Fuel Tax Pass-Through
* Hourly Prices & Two-Stage RDiT Approach (incl. Bootstrapping) - Revision
* Frederik von Waldow
* 10.04.2025
*###############################################################################
clear all
ssc install estout, replace
ssc install outreg2
ssc install coefplot, replace
net from http://www.stata-journal.com/software/sj3-2/
net describe st0039
net install st0039
graph set window fontface "Times New Roman"
*-------------------------------------------------------------------------------
log using "sta_spec_passthru", replace
*-------------------------------------------------------------------------------
* FIRST INTERVENTION
*-------------------------------------------------------------------------------
timer on 1

foreach i in mean_e5 mean_diesel {
	use "C:\Users\frede\Desktop\Bachelor-Thesis\Daten\tankstellenpreise\roadside_price_2022_0309_openandassume_revised.dta", clear
	cd "C:\Users\frede\Desktop\Code_TwoStage_RDiT_vWaldow\station_sec_passthru_rev"
	*drop stations in 10km distance to german boarder:
		drop if d_gborder_10_BEL ==1 | ///
					d_gborder_10_CHE  ==1 |  ///
					d_gborder_10_CZE  ==1 | ///
					d_gborder_10_DNK  ==1 | ///
					d_gborder_10_FRA  ==1 | ///
					d_gborder_10_LUX  ==1 | ///
					d_gborder_10_NLD  ==1 | ///
					d_gborder_10_POL  ==1 | ///
					d_gborder_10_AUT  ==1
	gen id_new = station_id
	xtset station_id time
	keep if inrange(time,9961,13321) // 20 Weeks (10 Weeks pre treatment, 10 Week with treatment)
	bysort station_id (time): gen L24_BRENT_FOB = L24.BRENT_FOB
	bysort station_id (time): gen L24_usd_to_euro = L24.usd_to_euro
	bysort station_id (time): gen L1_kfz_menge_B = L1.kfz_menge_B

	*** RDiT Indicators ***
	gen hour_to_treat_1 = time
	replace hour_to_treat_1 = (time- 11641) // Treatment June 1st 2022
	gen tankrabatt_1 = 0
	replace tankrabatt_1 = 1 if inrange(month,6,8)
	scalar da = 9 // Donut 9 hours (baseline)
	scalar db = -da
	scalar h = 336 // Bandwidth 336 hours (baseline)
	gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h)) // triangular kernel (non-parametric estimation)
		replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	global var_seas i.hour i.week_days // seasonality and price cycles covariates
	global var_sup BRENT_FOB L24_BRENT_FOB usd_to_euro L24_usd_to_euro // supply covariates
	global var_dem kfz_menge_B L1_kfz_menge_B temperature precipitation holiday_state day_before_holiday_state // demand covariates

	* First Stage
			reghdfe `i' $var_sup $var_dem $var_seas , absorb(station_id) resid(res_`i')
	keep if inrange(hour_to_treat_1, scalar(-h), scalar(h))
	bysort station_id: egen count_bef=count(res_`i') if inrange(hour_to_treat_1, scalar(-h),scalar(db))
	bysort station_id: egen count_bef_max = max(count_bef)
	bysort station_id: egen count_af=count(res_`i') if inrange(hour_to_treat_1,scalar(da),scalar(h))
	bysort station_id: egen count_af_max = max(count_af)
	drop count_bef count_af
	replace count_bef_max = 0 if count_bef_max ==.
	replace count_af_max = 0 if count_af_max ==.
	drop if count_af_max <25 | count_bef_max <25
	* Second Stage
	encode station_uuid, gen(station_id2)
			gen TR_`i'_first =.
			xtset station_id time
			distinct station_id2 if !mi(res_`i')
			forval x = 1/ `r(ndistinct)' {
				qui: reg res_`i' 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
									if station_id2 ==`x' & (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
					  replace TR_`i'_first = _b[1.tankrabatt_1] if station_id2 ==`x'
		}	
	save "TR_ind_hour_first_`i'.dta", replace
}

*-------------------------------------------------------------------------------
* SECOND INTERVENTION
*-------------------------------------------------------------------------------
foreach i in mean_diesel {
	use "C:\Users\frede\Desktop\Bachelor-Thesis\Daten\tankstellenpreise\roadside_price_2022_0309_openandassume_revised.dta", clear
	cd "C:\Users\frede\Desktop\Code_TwoStage_RDiT_vWaldow\station_sec_passthru_rev"
	*drop stations in 10km distance to german boarder:
		drop if d_gborder_10_BEL ==1 | ///
					d_gborder_10_CHE  ==1 |  ///
					d_gborder_10_CZE  ==1 | ///
					d_gborder_10_DNK  ==1 | ///
					d_gborder_10_FRA  ==1 | ///
					d_gborder_10_LUX  ==1 | ///
					d_gborder_10_NLD  ==1 | ///
					d_gborder_10_POL  ==1 | ///
					d_gborder_10_AUT  ==1
	gen id_new = station_id
	xtset station_id time
	keep if inrange(time,12169,15529) // 20 Weeks (10 Weeks pre treatment, 10 Week with treatment)
	bysort station_id (time): gen L24_BRENT_FOB = L24.BRENT_FOB
	bysort station_id (time): gen L24_usd_to_euro = L24.usd_to_euro
	bysort station_id (time): gen L1_kfz_menge_B = L1.kfz_menge_B

	*** RDiT Indicators ***
	gen hour_to_treat_2 = time
	replace hour_to_treat_2 = (time- 13849) // Treatment September 1st 2022
	gen tankrabatt_2 = 1
	replace tankrabatt_2 = 0 if inrange(month,6,8)
	scalar da = 9 // Donut 9 hours (baseline)
	scalar db = -da
	scalar h = 336 // Bandwidth 336 hours (baseline)
	gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h)) // triangular kernel (non-parametric estimation)
		replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	global var_seas i.hour i.week_days // seasonality and price cycles covariates
	global var_sup BRENT_FOB L24_BRENT_FOB usd_to_euro L24_usd_to_euro // supply covariates
	global var_dem kfz_menge_B L1_kfz_menge_B temperature precipitation holiday_state day_before_holiday_state // demand covariates

	* First Stage
			reghdfe `i' $var_sup $var_dem $var_seas , absorb(station_id) resid(res_`i')
	keep if inrange(hour_to_treat_2, scalar(-h), scalar(h))
	bysort station_id: egen count_bef=count(res_`i') if inrange(hour_to_treat_2, scalar(-h),scalar(db))
	bysort station_id: egen count_bef_max = max(count_bef)
	bysort station_id: egen count_af=count(res_`i') if inrange(hour_to_treat_2,scalar(da),scalar(h))
	bysort station_id: egen count_af_max = max(count_af)
	drop count_bef count_af
	replace count_bef_max = 0 if count_bef_max ==.
	replace count_af_max = 0 if count_af_max ==.
	drop if count_af_max <25 | count_bef_max <25
	* Second Stage
	encode station_uuid, gen(station_id2)
			gen TR_`i'_second =.
			xtset station_id time
			distinct station_id2 if !mi(res_`i')
			forval x = 1/ `r(ndistinct)' {
				qui: reg res_`i' 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
									if station_id2 ==`x' & (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
					  replace TR_`i'_second = _b[1.tankrabatt_2] if station_id2 ==`x'
		}	
	save "TR_ind_hour_second_`i'.dta", replace
}
timer off 1
timer list
log close