*###############################################################################
* DO-FILE: GRAPH COMBINE 
* Spatial Competition and Fuel Tax Pass-Through
* Hourly Prices & Two-Stage RDiT Approach (incl. Bootstrapping)
* Frederik von Waldow
* 26.08.2024
*###############################################################################
clear all
ssc install estout, replace
ssc install outreg2
ssc install reghdfe
ssc install ftools
ssc install coefplot, replace
ssc install rdrobust
net from http://www.stata-journal.com/software/sj3-2/
net describe st0039
net install st0039
graph set window fontface "Times New Roman"

cd "Z:\!_Analysen\vonWaldow_TaxIncidence_Market\Two_Stage_RDiT\graphs"
use "Z:\!_Analysen\vonWaldow_TaxIncidence_Market\Two_Stage_RDiT\roadside_price_2022_0309_openandassume.dta", clear
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
keep if inrange(time,9961,11809) // 11 Weeks (10 Weeks pre treatment, 1 Week with treatment)
gen lnBRENT_FOB = ln(BRENT_FOB)
bysort station_id (time): gen L24_lnBRENT_FOB = L24.lnBRENT_FOB
bysort station_id (time): gen L24_usd_to_euro = L24.usd_to_euro
bysort station_id (time): gen L1_kfz_menge_B = L1.kfz_menge_B
*encode NUTS, gen(NUTS3)

*** RDiT Indicators ***
gen hour_to_treat_1 = time
replace hour_to_treat_1 = (time- 11641) // Treatment June 1st 2022
gen tankrabatt_1 = 0
replace tankrabatt_1 = 1 if inrange(month,6,8)
scalar da = 9 // Donut 9 hours
scalar ba = 120 // Bandwidth 120 hours
scalar db = -da
scalar bb = -ba
global var_seas i.hour i.week_days
global var_sup lnBRENT_FOB L24_lnBRENT_FOB usd_to_euro L24_usd_to_euro	 
global var_dem kfz_menge_B L1_kfz_menge_B temperature precipitation holiday_state day_before_holiday_state


*** RDiT - PLOTS ***
/*
foreach i in mean_e10 mean_diesel {
	reghdfe `i' $var_seas, absorb(station_id) resid(res_`i'_seas)
	rdplot res_`i'_seas hour_to_treat_1 if  (inrange(hour_to_treat_1, -120,-9) | inrange(hour_to_treat_1,9,120)), c(0) h(120) p(1) nbins(28 28) masspoints(adjust) graph_options(ytitle(Residuals `i' (ct/l)) xtitle(Hours to Treatment) legend(off) scale(1.3) name(`i'_first, replace))
	graph save rd_plot_first_`i', replace
	drop res_`i'_seas
}	
*/

*** HHI and Market Boundaries ***
* First Intervention
scalar scal_tax_mean_diesel = 1/16.71 * 100 // (1/+-16.71)*100 = +-5.984440455
scalar scal_tax_mean_e5 = 1/35.16 * 100 // (1/+-35.16)*100 = +-2.844141069
scalar scal_tax_mean_e10 = 1/35.16 * 100 // (1/+-35.16)*100 = +-2.844141069
foreach i in e5 e10 diesel {
	* tax reduction
		cd "Z:\!_Analysen\vonWaldow_TaxIncidence_Market\Two_Stage_RDiT\FIRST_`i'"
		estimates use "E3_first_`i'_hhi_rad_3"
			scalar scal_tankrabatt_3 = _b[1.tankrabatt_1]
			scalar scal_tankrabatt_3_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_`i'
		estimates use "E4_first_`i'_hhi_rad_4" 
			scalar scal_tankrabatt_4 = _b[1.tankrabatt_1]
			scalar scal_tankrabatt_4_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_`i'
		estimates use "E5_first_`i'_hhi_rad_5"
			scalar scal_tankrabatt_5 = _b[1.tankrabatt_1]
			scalar scal_tankrabatt_5_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_`i'
		graph use "HHI_SENS_first_`i'_per", scheme(stgcolor) play(Thick_Dash)
		cd "Z:\!_Analysen\vonWaldow_TaxIncidence_Market\Two_Stage_RDiT\graphs"
		graph export HHI_SENS_first_`i'_per.png, replace
	* tax increase
		cd "Z:\!_Analysen\vonWaldow_TaxIncidence_Market\Two_Stage_RDiT\SECOND_`i'"
		estimates use "E3_second_`i'_hhi_rad_3"
			scalar scal_tankrabatt_3 = _b[1.tankrabatt_1]
			scalar scal_tankrabatt_3_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_`i'
		estimates use "E4_second_`i'_hhi_rad_4"
			scalar scal_tankrabatt_4 = _b[1.tankrabatt_1]
			scalar scal_tankrabatt_4_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_`i'
		estimates use "E5_second_`i'_hhi_rad_5"
			scalar scal_tankrabatt_5 = _b[1.tankrabatt_1]
			scalar scal_tankrabatt_5_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_`i'
		graph use "HHI_SENS_second_`i'_per", scheme(stgcolor) play(Thick_Dash)
		cd "Z:\!_Analysen\vonWaldow_TaxIncidence_Market\Two_Stage_RDiT\graphs"
		graph export HHI_SENS_second_`i'_per.png, replace
	}


*** HHI and Vertical Integration ***
foreach i in e5 e10 diesel {
	* tax reduction
		cd "Z:\!_Analysen\vonWaldow_TaxIncidence_Market\Two_Stage_RDiT\FIRST_`i'"
		estimates use "first_`i'_hhi_vert_1"
			scalar scal_tankrabatt_1 = _b[1.tankrabatt_1]
			scalar scal_tankrabatt_1_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_`i'
		estimates use "first_`i'_hhi_vert_2"
			scalar scal_tankrabatt_2 = _b[1.tankrabatt_1]
			scalar scal_tankrabatt_2_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_`i'
		estimates use "first_`i'_hhi_vert_3"
			scalar scal_tankrabatt_3 = _b[1.tankrabatt_1]
			scalar scal_tankrabatt_3_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_`i'
		graph use "HHI_VERT_first_`i'_per", scheme(stgcolor) play(Thick_Dash)
		cd "Z:\!_Analysen\vonWaldow_TaxIncidence_Market\Two_Stage_RDiT\graphs"
		graph export HHI_VERT_first_`i'_per.png, replace
	* tax increase
		cd "Z:\!_Analysen\vonWaldow_TaxIncidence_Market\Two_Stage_RDiT\SECOND_`i'"
		estimates use "second_`i'_hhi_vert_1"
			scalar scal_tankrabatt_1 = _b[1.tankrabatt_1]
			scalar scal_tankrabatt_1_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_`i'
		estimates use "second_`i'_hhi_vert_2"
			scalar scal_tankrabatt_2 = _b[1.tankrabatt_1]
			scalar scal_tankrabatt_2_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_`i'
		estimates use "second_`i'_hhi_vert_3"
			scalar scal_tankrabatt_3 = _b[1.tankrabatt_1]
			scalar scal_tankrabatt_3_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_`i'
		graph use "HHI_VERT_second_`i'_per", scheme(stgcolor) play(Thick_Dash)
		cd "Z:\!_Analysen\vonWaldow_TaxIncidence_Market\Two_Stage_RDiT\graphs"
		graph export HHI_VERT_second_`i'_per.png, replace
}

*** Station-Specific Pass-Through - Kernel Density smoothed Distribution ***
cd "C:\Users\frede\Desktop\Code_TwoStage_RDiT_vWaldow\station_sec_passthru_rev"
use "TR_ind_hour_first_mean_diesel.dta", clear
rename count_bef_max count_bef_max_first_diesel
rename count_af_max count_af_max_first_diesel
merge m:m station_uuid using "TR_ind_hour_first_mean_e5.dta", keepusing(TR_mean_e5_first count_bef_max count_af_max)
drop _merge
rename count_bef_max count_bef_max_first_e5
rename count_af_max count_af_max_first_e5
merge m:m station_uuid using "TR_ind_hour_second_mean_diesel.dta", keepusing(TR_mean_diesel_second count_bef_max count_af_max)
drop _merge
rename count_bef_max count_bef_max_second_diesel
rename count_af_max count_af_max_second_diesel
merge m:m station_uuid using "TR_ind_hour_second_mean_e5.dta", keepusing(TR_mean_e5_second count_bef_max count_af_max)
drop _merge
rename count_bef_max count_bef_max_second_e5
rename count_af_max count_af_max_second_e5
keep station_uuid station_id TR_mean_* count_*
duplicates drop

* Pass_Trhough in %
	gen passthrough_e5_first = (-TR_mean_e5_first/35.16)*100
	*gen passthrough_e10_first = (-TR_sta_mean_e10_first/35.16)*100
	gen passthrough_diesel_first = (-TR_mean_diesel_first/16.71)*100
	gen passthrough_e5_second = (TR_mean_e5_second/35.16)*100
	*gen passthrough_e10_second = (TR_sta_mean_e10_sec/35.16)*100
	gen passthrough_diesel_second = (TR_mean_diesel_second/16.71)*100

twoway (kdensity passthrough_e5_first if inrange(count_bef_max_first_e5,70,336) & inrange(count_af_max_first_e5,70,336) & inrange(passthrough_e5_first,-20,145), lwidth(medthick) kernel(epanechnikov) xtitle(Tax Pass-Through (%)) ytitle(Density) lcolor(stc1) legend(label(1 "Super E5 (Tax Reduction)"))) || ///
		   (kdensity passthrough_diesel_first if inrange(count_bef_max_first_diesel,70,336) & inrange(count_af_max_first_diesel,70,336) & inrange(passthrough_diesel_first,-20,145),lwidth(medthick) kernel(epanechnikov) lcolor(stc2) legend(label(2 "Diesel (Tax Reduction)"))) || ///
		   (kdensity passthrough_e5_second if inrange(count_bef_max_second_e5,70,336) & inrange(count_af_max_second_e5,70,336) & inrange(passthrough_e5_second,-20,145), lwidth(medthick) kernel(epanechnikov) lcolor(stc1) lpattern(dash) legend(label(3 "Super E5 (Tax Increase)"))) || ///
		   (kdensity passthrough_diesel_second if inrange(count_bef_max_second_diesel,70,336) & inrange(count_af_max_second_diesel,70,336) & inrange(passthrough_diesel_second,-20,145),lwidth(medthick) kernel(epanechnikov) lcolor(stc2) xlab(-25(25) 150) lpattern(dash) legend(col(2) pos(6) order(1 2 3 4) label(4 "Diesel (Tax Increase)")) xline(0 100) scale(1.18) xsize(20cm) ysize(10cm) name(graph_1, replace))
		   		   


