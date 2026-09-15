*###############################################################################
* DO-FILE (e5)
* Spatial Competition and Fuel Tax Pass-Through
* Hourly Prices & Two-Stage RDiT (Bootstrap) - Revision
* Frederik von Waldow
* 08.04.2024
*###############################################################################
clear all
set more off
ssc install estout, replace
ssc install outreg2
ssc install reghdfe
ssc install ftools
ssc install coefplot, replace
net from http://www.stata-journal.com/software/sj3-2/
net describe st0039
net install st0039
graph set window fontface "Times New Roman"

log using "BOOT_SECOND_E5", replace

*-------------------------------------------------------------------------------
* Bootstrap Iterations:
global boot_rep reps(1000)
*-------------------------------------------------------------------------------

*###############################################################################
* ROADSIDE STATIONS
*###############################################################################

*-------------------------------------------------------------------------------
* SECOND INTERVENTION
*-------------------------------------------------------------------------------
cd "X:\prj-avr\!_Analysen\vonWaldow_TaxIncidence_Market\Two_Stage_RDiT_Revision\SECOND_e5"
use "X:\prj-avr\!_Analysen\vonWaldow_TaxIncidence_Market\Two_Stage_RDiT_Revision\roadside_price_2022_0309_openandassume_revised.dta", clear

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

* keep relevant variables
drop mean_diesel mean_e10 /// 
     d_gborder_10_BEL d_gborder_10_CHE d_gborder_10_CZE d_gborder_10_DNK d_gborder_10_FRA d_gborder_10_LUX d_gborder_10_NLD d_gborder_10_POL d_gborder_10_AUT /// 
	 own_share_3km totalSize_4km own_share_4km own_share_5km /// 
	 month land station_uuid NUTS
	 
*-------------------------------------------------------------------------------
*FUEL TYPE : e5*
*------------------------------------------------------------------------------- 
timer on 1
*######################### AVERAGE EFFECTS ####################################*
	* Seasonality
est clear
	capture program drop _all
		program second_e5_season
			reghdfe mean_e5 $var_seas , absorb(station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5
		end 	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_season_file, replace): second_e5_season
	estimates save second_e5_season , replace
	* Supply & Seasonality
		program second_e5_sup
			reghdfe mean_e5 $var_sup $var_seas , absorb(station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new) saving(second_e5_sup_file, replace): second_e5_sup
	estimates save second_e5_sup , replace
	* Demand & Seasonality
		program second_e5_dem
			reghdfe mean_e5 $var_dem $var_seas , absorb(station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5
		end
		xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new) saving(second_e5_sem_file, replace): second_e5_dem
	estimates save second_e5_dem , replace
	* Supply & Demand & Seasonality
	capture program drop _all
		program second_e5_allcov
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_allcov_file, replace): second_e5_allcov
	estimates save second_e5_allcov , replace
	esttab using "Main_e5_second_boot.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Main specification e5 second Intervention) nomtitles addnotes(Notes follow)

*######################### HETEROGENEOUS EFFECTS ##############################*
est clear
gen d_comp =0
replace d_comp = 1 if inrange(d_next_competitor,1,3)
replace d_comp = 2 if d_next_competitor >3
xtile hhi_km_quantile = hhi_4km, nquantiles(4)
	*Vertical Integration:
	capture program drop _all
		program second_e5_vert_int
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.brand_vert_int [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_vert_int_file, replace): second_e5_vert_int
	estimates save second_e5_vert_int , replace
	* Region Type:
	capture program drop _all
		program second_e5_reg_typ
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.urban_rural_cat [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_reg_typ_file, replace): second_e5_reg_typ
	estimates save second_e5_reg_typ , replace
	* Distance next competitor:
	capture program drop _all
		program second_e5_dist_comp
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.d_comp [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_dist_comp_file, replace): second_e5_dist_comp
	estimates save second_e5_dist_comp , replace
	* HHI Quartiles (4km):
		capture program drop _all
		program second_e5_hhi_4km
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_hhi_4km_file, replace): second_e5_hhi_4km
	estimates save E4_second_e5_hhi_rad_4 , replace
	esttab using "Het_e5_second_boot.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Heterogeneous Effects second Intervention e5) keep(1.tankrabatt_2 1.tankrabatt_2#2.brand_vert_int 1.tankrabatt_2#3.brand_vert_int 1.tankrabatt_2#2.urban_rural_cat 1.tankrabatt_2#3.urban_rural_cat 1.tankrabatt_2#2.hhi_km_quantile 1.tankrabatt_2#3.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile 1.tankrabatt_2#1.d_comp 1.tankrabatt_2#2.d_comp _cons) nomtitles addnotes(Notes follow)

*######################### HHI and VERT-INT ###################################*
est clear
scalar scal_tax_mean_e5 = 1/35.16 * 100 // (1/+-35.16)*100 = +-2.844141069

	* OLIGOPOLISTIC PLAYER
	capture program drop _all
		program second_e5_hhi_vert_1 
			reghdfe mean_e5 $var_sup $var_dem $var_seas if brand_vert_int == 1, absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if brand_vert_int == 1 & (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_hhi_vert_1_file, replace): second_e5_hhi_vert_1
	ereturn list
	estimates store E_1_second_e5
	estimates save second_e5_hhi_vert_1, replace
	scalar scal_tankrabatt_1 = _b[1.tankrabatt_2]
	scalar scal_tankrabatt_1_per = abs(_b[1.tankrabatt_2])*scal_tax_mean_e5
	
	* OTHER INTEGRATED PLAYER
	capture program drop _all
		program second_e5_hhi_vert_2 
			reghdfe mean_e5 $var_sup $var_dem $var_seas if brand_vert_int == 2, absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if brand_vert_int == 2 & (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_hhi_vert_2_file, replace): second_e5_hhi_vert_2
	ereturn list
	estimates store E_2_second_e5
	estimates save second_e5_hhi_vert_2, replace
	scalar scal_tankrabatt_2 = _b[1.tankrabatt_2]
	scalar scal_tankrabatt_2_per = abs(_b[1.tankrabatt_2])*scal_tax_mean_e5
	
	* INDEPENDENT PLAYER
	capture program drop _all
		program second_e5_hhi_vert_3 
			reghdfe mean_e5 $var_sup $var_dem $var_seas if brand_vert_int == 3, absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if brand_vert_int == 3 & (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_hhi_vert_3_file, replace): second_e5_hhi_vert_3
	ereturn list
	estimates store E_3_second_e5
	estimates save second_e5_hhi_vert_3 , replace
	scalar scal_tankrabatt_3 = _b[1.tankrabatt_2]
	scalar scal_tankrabatt_3_per = abs(_b[1.tankrabatt_2])*scal_tax_mean_e5
		esttab using "HHI_VERT_e5_second_boot.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Heterogeneous Effects second Intervention) nomtitles addnotes(Notes follow)
	
coefplot (E_1_second_e5, transform(1.tankrabatt_2#2.hhi_km_quantile 1.tankrabatt_2#3.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile:= @+scal_tankrabatt_1) recast(bar) base(scal_tankrabatt_1) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(oligopolistic)) /// 
		 (E_2_second_e5, transform(1.tankrabatt_2#2.hhi_km_quantile 1.tankrabatt_2#3.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile:= @+scal_tankrabatt_2) recast(bar) base(scal_tankrabatt_2) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(other integrated)) /// 
		 (E_3_second_e5, transform(1.tankrabatt_2#2.hhi_km_quantile 1.tankrabatt_2#3.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile:= @+scal_tankrabatt_3) recast(bar) base(scal_tankrabatt_3) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(independent)), /// 
												  vertical /// 
												  keep(1.tankrabatt_2 1.tankrabatt_2#2.hhi_km_quantile 1.tankrabatt_2#3.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile) /// 
												  legend(pos(6) col(3)) xlabel(1 "({it:Ref.}) Q1" 2 "Q2" 3 "Q3" 4 "Q4") /// 
												  ytitle("Tax Pass-Through E5 (ct/l)") xtitle("HHI Quartiles") /// 
												  addplot(function second=scal_tankrabatt_1, ra(0.5 4.5) lwidth(medthick) lcolor(stc1) lpattern(shortdash) || /// 
												  function second=scal_tankrabatt_2, ra(0.5 4.5) lwidth(medthick) lcolor(stc2) lpattern(shortdash) || ///
												  function thrid=scal_tankrabatt_3, ra(0.5 4.5) lwidth(medthick) lcolor(stc3) lpattern(shortdash)) /// 
												  scale(1.3) ysize(5cm) xsize(8cm) name(e5_brand_hhi_second, replace)
	graph save "HHI_VERT_second_e5", replace
coefplot (E_1_second_e5,  rescale(2.844141069) transform(1.tankrabatt_2#2.hhi_km_quantile 1.tankrabatt_2#3.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile:= @+scal_tankrabatt_1_per) recast(bar) base(scal_tankrabatt_1_per) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(oligopolistic)) /// 
		 (E_2_second_e5, rescale(2.844141069) transform(1.tankrabatt_2#2.hhi_km_quantile 1.tankrabatt_2#3.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile:= @+scal_tankrabatt_2_per) recast(bar) base(scal_tankrabatt_2_per) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(other integrated)) /// 
		 (E_3_second_e5, rescale(2.844141069) transform(1.tankrabatt_2#2.hhi_km_quantile 1.tankrabatt_2#3.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile:= @+scal_tankrabatt_3_per) recast(bar) base(scal_tankrabatt_3_per) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(independent)), /// 
												  vertical /// 
												  keep(1.tankrabatt_2 1.tankrabatt_2#2.hhi_km_quantile 1.tankrabatt_2#3.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile) /// 
												  legend(pos(6) col(3)) xlabel(1 "({it:Ref.}) Q1" 2 "Q2" 3 "Q3" 4 "Q4") /// 
												  ytitle("Tax Pass-Through E5 (%)") xtitle("HHI Quartiles") /// 
												  addplot(function second=scal_tankrabatt_1_per, ra(0.5 4.5) lwidth(medthick) lcolor(stc1) lpattern(shortdash) || /// 
												  function second=scal_tankrabatt_2_per, ra(0.5 4.5) lwidth(medthick) lcolor(stc2) lpattern(shortdash) || ///
												  function thrid=scal_tankrabatt_3_per, ra(0.5 4.5) lwidth(medthick) lcolor(stc3) lpattern(shortdash)) /// 
												  scale(1.3) ysize(5cm) xsize(8cm) name(e5_brand_hhi_second_per, replace)
	graph save "HHI_VERT_second_e5_per", replace
	graph export "HHI_VERT_second_e5_per.emf", replace
	
*##################### SENSITIVITY (HHI) - MARKET BOUNDARIES #################*
est clear
	* 3km RADIUS
	drop hhi_km_quantile
	xtile hhi_km_quantile = hhi_3km, nquantiles(4)
	capture program drop _all
		program second_e5_hhi_rad_3 
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_hhi_rad_3_file, replace): second_e5_hhi_rad_3 
	ereturn list
	estimates store E3_second_e5_hhi_rad_3
	estimates save E3_second_e5_hhi_rad_3, replace
	scalar scal_tankrabatt_3 = _b[1.tankrabatt_2]
	scalar scal_tankrabatt_3_per = abs(_b[1.tankrabatt_2])*scal_tax_mean_e5
	
	* 4km RADIUS
	estimates use "E4_second_e5_hhi_rad_4"
	eststo: estimates replay "E4_second_e5_hhi_rad_4"
	estimates store E4_second_e5_hhi_rad_4
	scalar scal_tankrabatt_4 = _b[1.tankrabatt_2]
	scalar scal_tankrabatt_4_per = abs(_b[1.tankrabatt_2])*scal_tax_mean_e5
	
	* 5km RADIUS
	drop hhi_km_quantile
	xtile hhi_km_quantile = hhi_5km, nquantiles(4)
	capture program drop _all
		program second_e5_hhi_rad_5 
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_hhi_rad_5_file, replace): second_e5_hhi_rad_5 
	ereturn list
	estimates store E5_second_e5_hhi_rad_5
	estimates save E5_second_e5_hhi_rad_5, replace
	scalar scal_tankrabatt_5 = _b[1.tankrabatt_2]
	scalar scal_tankrabatt_5_per = abs(_b[1.tankrabatt_2])*scal_tax_mean_e5
		esttab using "HHI_SENS_second_boot.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Heterogeneous Effects second Intervention) nomtitles addnotes(Notes follow)
	
coefplot (E3_second_e5_hhi_rad_3, transform(1.tankrabatt_2#2.hhi_km_quantile 1.tankrabatt_2#3.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile:= @+scal_tankrabatt_3) recast(bar) base(scal_tankrabatt_3) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(3km Radius)) /// 
		 (E4_second_e5_hhi_rad_4, transform(1.tankrabatt_2#2.hhi_km_quantile 1.tankrabatt_2#3.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile:= @+scal_tankrabatt_4) recast(bar) base(scal_tankrabatt_4) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(4km Radius)) /// 
		 (E5_second_e5_hhi_rad_5, transform(1.tankrabatt_2#2.hhi_km_quantile 1.tankrabatt_2#3.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile:= @+scal_tankrabatt_5) recast(bar) base(scal_tankrabatt_5) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(5km Radius)), /// 
												  vertical /// 
												  keep(1.tankrabatt_2 1.tankrabatt_2#2.hhi_km_quantile 1.tankrabatt_2#3.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile) /// 
												  legend(pos(6) col(3)) xlabel(1 "({it:Ref.}) Q1" 2 "Q2" 3 "Q3" 4 "Q4") /// 
												  ytitle("Tax Pass-Through E5 (ct/l)") xtitle("HHI Quartiles") /// 
												  addplot(function second=scal_tankrabatt_3, ra(0.5 4.5) lwidth(medthick) lcolor(stc1) lpattern(shortdash) || /// 
												  function second=scal_tankrabatt_4, ra(0.5 4.5) lwidth(medthick) lcolor(stc2) lpattern(shortdash) || ///
												  function thrid=scal_tankrabatt_5, ra(0.5 4.5) lwidth(medthick) lcolor(stc3) lpattern(shortdash)) /// 
												  scale(1.3) ysize(5cm) xsize(8cm) name(e5_sens_hhi_second, replace)
	graph save "HHI_SENS_second_e5", replace
coefplot (E3_second_e5_hhi_rad_3,  rescale(2.844141069) transform(1.tankrabatt_2#2.hhi_km_quantile 1.tankrabatt_2#3.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile:= @+scal_tankrabatt_3_per) recast(bar) base(scal_tankrabatt_3_per) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(3km Radius)) /// 
		 (E4_second_e5_hhi_rad_4, rescale(2.844141069) transform(1.tankrabatt_2#2.hhi_km_quantile 1.tankrabatt_2#3.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile:= @+scal_tankrabatt_4_per) recast(bar) base(scal_tankrabatt_4_per) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(4km Radius)) /// 
		 (E5_second_e5_hhi_rad_5, rescale(2.844141069) transform(1.tankrabatt_2#2.hhi_km_quantile 1.tankrabatt_2#3.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile:= @+scal_tankrabatt_5_per) recast(bar) base(scal_tankrabatt_5_per) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(5km Radius)), /// 
												  vertical /// 
												  keep(1.tankrabatt_2 1.tankrabatt_2#2.hhi_km_quantile 1.tankrabatt_2#3.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile) /// 
												  legend(pos(6) col(3)) xlabel(1 "({it:Ref.}) Q1" 2 "Q2" 3 "Q3" 4 "Q4") /// 
												  ytitle("Tax Pass-Through E5 (%)") xtitle("HHI Quartiles") /// 
												  addplot(function second=scal_tankrabatt_3_per, ra(0.5 4.5) lwidth(medthick) lcolor(stc1) lpattern(shortdash) || /// 
												  function second=scal_tankrabatt_4_per, ra(0.5 4.5) lwidth(medthick) lcolor(stc2) lpattern(shortdash) || ///
												  function thrid=scal_tankrabatt_5_per, ra(0.5 4.5) lwidth(medthick) lcolor(stc3) lpattern(shortdash)) /// 
												  scale(1.3) ysize(5cm) xsize(8cm) name(e5_sens_hhi_second_per, replace)
	graph save "HHI_SENS_second_e5_per", replace	
	graph export "HHI_SENS_second_e5_per.emf", replace

* #################### ROBUSTNESS - BANDWIDTH & DONUT #########################*
est clear
	* 240 hour (10 day) bandwidth
	scalar h = 240 // Bandwidth 240 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program second_10days
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_10days, replace): second_10days
	estimates save second_10days, replace
	* 288 hour (12 day) bandwidth
	scalar h = 288 // Bandwidth 288 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program second_12days
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_12days, replace): second_12days
	estimates save second_12days, replace
	* 336 hour (14 day) bandwidth (no donut)
	scalar da = 0 // Donut 0 hours
	scalar db = -da
	scalar h = 336 // Bandwidth 336 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program second_14days_0d
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_14days_0d, replace): second_14days_0d
	estimates save second_14days_0d, replace
	* 336 hour bandwidth (9 hours donut) siehe Average Effects
	estimates use "second_e5_allcov"
	eststo: estimates replay "second_e5_allcov"
	* 336 hour (14 day) bandwidth (18 hours donut)
	scalar da = 18 // Donut 18 hours
	scalar db = -da
	scalar h = 336 // Bandwidth 336 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program second_14days_18d
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_14days_18d, replace): second_14days_18d
	estimates save second_14days_18d, replace
	* 384 hour (16 day) bandwidth
	scalar da = 9 // Donut 9 hours
	scalar db = -da
	scalar h = 384 // Bandwidth 384 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program second_16days
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_16days, replace): second_16days
	estimates save second_16days, replace
		esttab using "ROBUST_e5_second_boot.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Heterogeneous Effects second Intervention) nomtitles addnotes(Notes follow)

*######################### BRAND AFFILIATION ##################################*
	scalar h = 336 // Bandwidth 336 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
est clear	
	* ARAL
	capture program drop _all
		program second_e5_ARAL
			reghdfe mean_e5 $var_sup $var_dem $var_seas if brand == "ARAL", absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if brand == "ARAL" & (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_ARAL, replace): second_e5_ARAL
	estimates save second_e5_ARAL , replace
	* SHELL
	capture program drop _all
		program second_e5_SHELL
			reghdfe mean_e5 $var_sup $var_dem $var_seas if brand == "SHELL", absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if brand == "SHELL" & (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_SHELL, replace): second_e5_SHELL
	estimates save second_e5_SHELL , replace
	* TOTALENERGIES
	capture program drop _all
		program second_e5_TOTALEnergies
			reghdfe mean_e5 $var_sup $var_dem $var_seas if brand == "TOTALEnergies", absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if brand == "TOTALEnergies" & (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_TOTALEnergies, replace): second_e5_TOTALEnergies
	estimates save second_e5_TOTALEnergies , replace
	* ESSO
	capture program drop _all
		program second_e5_ESSO
			reghdfe mean_e5 $var_sup $var_dem $var_seas if brand == "ESSO", absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if brand == "ESSO" & (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_ESSO, replace): second_e5_ESSO
	estimates save second_e5_ESSO , replace
	* JET
	capture program drop _all
		program second_e5_JET
			reghdfe mean_e5 $var_sup $var_dem $var_seas if brand == "JET", absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if brand == "JET" & (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_JET, replace): second_e5_JET
	estimates save second_e5_JET , replace
	* BFT
	capture program drop _all
		program second_e5_BFT
			reghdfe mean_e5 $var_sup $var_dem $var_seas if brand == "BFT", absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if brand == "BFT" & (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_BFT, replace): second_e5_BFT
	estimates save second_e5_BFT , replace
	* AVIA
	capture program drop _all
		program second_e5_AVIA
			reghdfe mean_e5 $var_sup $var_dem $var_seas if brand == "AVIA", absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if brand == "AVIA" & (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_AVIA, replace): second_e5_AVIA
	estimates save second_e5_AVIA , replace
	esttab using "BRAND_e5_second_boot.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Main specification e5 second Intervention) nomtitles addnotes(Notes follow)
		
*#################### EDGEWORTH CYCLES (no hour FE) ###########################*
drop hhi_km_quantile
xtile hhi_km_quantile = hhi_4km, nquantiles(4)
est clear	
	* Average Effect
	capture program drop _all
		program second_e5_allcov_nohour
			reghdfe mean_e5 $var_sup $var_dem i.week_days, absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
									if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_allcov_file, replace): second_e5_allcov_nohour
	estimates save second_e5_nohour , replace
	* HHI Quartiles (4km):
	capture program drop _all
		program second_e5_hhi_4km_nohour
			reghdfe mean_e5 $var_sup $var_dem i.week_days, absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_hhi_4km_file, replace): second_e5_hhi_4km_nohour
	estimates save second_e5_nohour_hhi , replace
	* Average Effect: Morning Peak (05:00 to 10:59)
	capture program drop _all
		program second_e5_nohour_5_10
			reghdfe mean_e5 $var_sup $var_dem i.week_days, absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h))) & inrange(hour,5,10)
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_nohour_5_10_file, replace): second_e5_nohour_5_10
	estimates save second_e5_nohour_5_10, replace
	* Average Effect: Mid Day (11:00 to 16:59)
	capture program drop _all
		program second_e5_nohour_11_16
			reghdfe mean_e5 $var_sup $var_dem i.week_days, absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h))) & inrange(hour,11,16)
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_nohour_11_16_file, replace): second_e5_nohour_11_16
	estimates save second_e5_nohour_11_16, replace
	* Average Effect: Evening Off-Peak (17:00 to 22:59)
	capture program drop _all
		program second_e5_nohour_17_22
			reghdfe mean_e5 $var_sup $var_dem i.week_days, absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h))) & inrange(hour,17,22)
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_nohour_17_22_file, replace): second_e5_nohour_17_22
	estimates save second_e5_nohour_17_22, replace
	* Average Effect: Night (23:00 to 04:59)
	capture program drop _all
		program second_e5_nohour_23_4
			reghdfe mean_e5 $var_sup $var_dem i.week_days, absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h))) & (inrange(hour,0,4) | hour == 23)
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_nohour_23_4_file, replace): second_e5_nohour_23_4
	estimates save second_e5_nohour_23_4, replace
	* HHI Quartiles (4km): Morning Peak (05:00 to 10:59)
	capture program drop _all
		program second_e5_nohour_5_10_hhi
			reghdfe mean_e5 $var_sup $var_dem i.week_days, absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h))) & inrange(hour,5,10)
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_nohour_5_10_file_hhi, replace): second_e5_nohour_5_10_hhi
	estimates save second_e5_nohour_5_10_hhi, replace
	* HHI Quartiles (4km): Mid Day (11:00 to 16:59)
	capture program drop _all
		program second_e5_nohour_11_16_hhi
			reghdfe mean_e5 $var_sup $var_dem i.week_days, absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h))) & inrange(hour,11,16)
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_nohour_11_16_file_hhi, replace): second_e5_nohour_11_16_hhi
	estimates save second_e5_nohour_11_16_hhi, replace
	* HHI Quartiles (4km): Evening Off-Peak (17:00 to 22:59)
	capture program drop _all
		program second_e5_nohour_17_22_hhi
			reghdfe mean_e5 $var_sup $var_dem i.week_days, absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h))) & inrange(hour,17,22)
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_nohour_17_22_file_hhi, replace): second_e5_nohour_17_22_hhi
	estimates save second_e5_nohour_17_22_hhi, replace
	* HHI Quartiles (4km): Night (23:00 to 04:59)
	capture program drop _all
		program second_e5_nohour_23_4_hhi
			reghdfe mean_e5 $var_sup $var_dem i.week_days, absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h))) & (inrange(hour,0,4) | hour == 23)
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_nohour_23_4_file_hhi, replace): second_e5_nohour_23_4_hhi
	estimates save second_e5_nohour_23_4_hhi, replace
	esttab using "nohour_e5_second_boot.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) keep(1.tankrabatt_2 1.tankrabatt_2#2.hhi_km_quantile 1.tankrabatt_2#3.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile _cons) title(NO hourly FE e5 second Intervention) nomtitles addnotes(Notes follow)

*################### HHI and VERT-INT (one regression) ########################*
est clear	
	* Average Effect
	capture program drop _all
		program second_e5_het_included
			reghdfe mean_e5 $var_sup $var_dem i.week_days, absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.hhi_km_quantile##i.brand_vert_int [aw=kernel_wgt] ///
									if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_het_included, replace): second_e5_het_included
	estimates save second_e5_het_included , replace
	esttab using "hhi_vert_tog_e5_second_boot.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) /// 
	keep(1.tankrabatt_2 /// 
	     2.brand_vert_int 3.brand_vert_int /// 
		 2.hhi_km_quantile 3.hhi_km_quantile 4.hhi_km_quantile ///
	     1.tankrabatt_2#2.hhi_km_quantile 1.tankrabatt_2#3.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile /// 
		 1.tankrabatt_2#2.brand_vert_int 1.tankrabatt_2#3.brand_vert_int /// 
		 1.tankrabatt_2#2.hhi_km_quantile#2.brand_vert_int 1.tankrabatt_2#2.hhi_km_quantile#3.brand_vert_int /// 
		 1.tankrabatt_2#3.hhi_km_quantile#2.brand_vert_int 1.tankrabatt_2#3.hhi_km_quantile#3.brand_vert_int ///
		 1.tankrabatt_2#4.hhi_km_quantile#2.brand_vert_int 1.tankrabatt_2#4.hhi_km_quantile#3.brand_vert_int ///
		 _cons) title(HHI and VERT e5 second Intervention) nomtitles addnotes(Notes follow)

			
timer off 1
timer list
log close