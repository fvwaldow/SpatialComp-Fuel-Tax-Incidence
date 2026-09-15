*###############################################################################
* DO-FILE (e10)
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

log using "BOOT_FIRST_E10", replace

*-------------------------------------------------------------------------------
* Bootstrap Iterations:
global boot_rep reps(1000)
*-------------------------------------------------------------------------------

*###############################################################################
* ROADSIDE STATIONS
*###############################################################################

*-------------------------------------------------------------------------------
* FIRST INTERVENTION
*-------------------------------------------------------------------------------
cd "X:\prj-avr\!_Analysen\vonWaldow_TaxIncidence_Market\Two_Stage_RDiT_Revision\FIRST_e10"
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

* keep relevant variables
drop mean_diesel mean_e5 /// 
     d_gborder_10_BEL d_gborder_10_CHE d_gborder_10_CZE d_gborder_10_DNK d_gborder_10_FRA d_gborder_10_LUX d_gborder_10_NLD d_gborder_10_POL d_gborder_10_AUT /// 
	 own_share_3km totalSize_4km own_share_4km own_share_5km /// 
	 month land station_uuid NUTS
	 
*-------------------------------------------------------------------------------
*FUEL TYPE : e10*
*------------------------------------------------------------------------------- 
timer on 1
*######################### AVERAGE EFFECTS ####################################*
	* Seasonality
est clear
	capture program drop _all
		program first_e10_season
			reghdfe mean_e10 $var_seas , absorb(station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10
		end 	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_season_file, replace): first_e10_season
	estimates save first_e10_season , replace
	* Supply & Seasonality
		program first_e10_sup
			reghdfe mean_e10 $var_sup $var_seas , absorb(station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new) saving(first_e10_sup_file, replace): first_e10_sup
	estimates save first_e10_sup , replace
	* Demand & Seasonality
		program first_e10_dem
			reghdfe mean_e10 $var_dem $var_seas , absorb(station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10
		end
		xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new) saving(first_e10_sem_file, replace): first_e10_dem
	estimates save first_e10_dem , replace
	* Supply & Demand & Seasonality
	capture program drop _all
		program first_e10_allcov
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_allcov_file, replace): first_e10_allcov
	estimates save first_e10_allcov , replace
	esttab using "Main_e10_first_boot.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Main specification e10 First Intervention) nomtitles addnotes(Notes follow)

*######################### HETEROGENEOUS EFFECTS ##############################*
est clear
gen d_comp =0
replace d_comp = 1 if inrange(d_next_competitor,1,3)
replace d_comp = 2 if d_next_competitor >3
xtile hhi_km_quantile = hhi_4km, nquantiles(4)
	*Vertical Integration:
	capture program drop _all
		program first_e10_vert_int
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.brand_vert_int [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_vert_int_file, replace): first_e10_vert_int
	estimates save first_e10_vert_int , replace
	* Region Type:
	capture program drop _all
		program first_e10_reg_typ
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.urban_rural_cat [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_reg_typ_file, replace): first_e10_reg_typ
	estimates save first_e10_reg_typ , replace
	* Distance next competitor:
	capture program drop _all
		program first_e10_dist_comp
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.d_comp [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_dist_comp_file, replace): first_e10_dist_comp
	estimates save first_e10_dist_comp , replace
	* HHI Quartiles (4km):
		capture program drop _all
		program first_e10_hhi_4km
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_hhi_4km_file, replace): first_e10_hhi_4km
	estimates save E4_first_e10_hhi_rad_4 , replace
	esttab using "Het_e10_first_boot.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Heterogeneous Effects First Intervention e10) keep(1.tankrabatt_1 1.tankrabatt_1#2.brand_vert_int 1.tankrabatt_1#3.brand_vert_int 1.tankrabatt_1#2.urban_rural_cat 1.tankrabatt_1#3.urban_rural_cat 1.tankrabatt_1#2.hhi_km_quantile 1.tankrabatt_1#3.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile 1.tankrabatt_1#1.d_comp 1.tankrabatt_1#2.d_comp _cons) nomtitles addnotes(Notes follow)

*######################### HHI and VERT-INT ###################################*
est clear
scalar scal_tax_mean_e10 = 1/35.16 * 100 // (1/+-35.16)*100 = +-2.844141069

	* OLIGOPOLISTIC PLAYER
	capture program drop _all
		program first_e10_hhi_vert_1 
			reghdfe mean_e10 $var_sup $var_dem $var_seas if brand_vert_int == 1, absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if brand_vert_int == 1 & (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_hhi_vert_1_file, replace): first_e10_hhi_vert_1
	ereturn list
	estimates store E_1_first_e10
	estimates save first_e10_hhi_vert_1, replace
	scalar scal_tankrabatt_1 = _b[1.tankrabatt_1]
	scalar scal_tankrabatt_1_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_e10
	
	* OTHER INTEGRATED PLAYER
	capture program drop _all
		program first_e10_hhi_vert_2 
			reghdfe mean_e10 $var_sup $var_dem $var_seas if brand_vert_int == 2, absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if brand_vert_int == 2 & (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_hhi_vert_2_file, replace): first_e10_hhi_vert_2
	ereturn list
	estimates store E_2_first_e10
	estimates save first_e10_hhi_vert_2, replace
	scalar scal_tankrabatt_2 = _b[1.tankrabatt_1]
	scalar scal_tankrabatt_2_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_e10
	
	* INDEPENDENT PLAYER
	capture program drop _all
		program first_e10_hhi_vert_3 
			reghdfe mean_e10 $var_sup $var_dem $var_seas if brand_vert_int == 3, absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if brand_vert_int == 3 & (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_hhi_vert_3_file, replace): first_e10_hhi_vert_3
	ereturn list
	estimates store E_3_first_e10
	estimates save first_e10_hhi_vert_3 , replace
	scalar scal_tankrabatt_3 = _b[1.tankrabatt_1]
	scalar scal_tankrabatt_3_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_e10
		esttab using "HHI_VERT_e10_first_boot.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Heterogeneous Effects First Intervention) nomtitles addnotes(Notes follow)
	
coefplot (E_1_first_e10, transform(1.tankrabatt_1#2.hhi_km_quantile 1.tankrabatt_1#3.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile:= @+scal_tankrabatt_1) recast(bar) base(scal_tankrabatt_1) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(oligopolistic)) /// 
		 (E_2_first_e10, transform(1.tankrabatt_1#2.hhi_km_quantile 1.tankrabatt_1#3.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile:= @+scal_tankrabatt_2) recast(bar) base(scal_tankrabatt_2) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(other integrated)) /// 
		 (E_3_first_e10, transform(1.tankrabatt_1#2.hhi_km_quantile 1.tankrabatt_1#3.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile:= @+scal_tankrabatt_3) recast(bar) base(scal_tankrabatt_3) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(independent)), /// 
												  vertical /// 
												  keep(1.tankrabatt_1 1.tankrabatt_1#2.hhi_km_quantile 1.tankrabatt_1#3.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile) /// 
												  legend(pos(6) col(3)) xlabel(1 "({it:Ref.}) Q1" 2 "Q2" 3 "Q3" 4 "Q4") /// 
												  ytitle("Tax Pass-Through E10 (ct/l)") xtitle("HHI Quartiles") /// 
												  addplot(function first=scal_tankrabatt_1, ra(0.5 4.5) lwidth(medthick) lcolor(stc1) lpattern(shortdash) || /// 
												  function second=scal_tankrabatt_2, ra(0.5 4.5) lwidth(medthick) lcolor(stc2) lpattern(shortdash) || ///
												  function thrid=scal_tankrabatt_3, ra(0.5 4.5) lwidth(medthick) lcolor(stc3) lpattern(shortdash)) /// 
												  scale(1.3) ysize(5cm) xsize(8cm) name(e10_brand_hhi_first, replace)
	graph save "HHI_VERT_first_e10", replace
coefplot (E_1_first_e10,  rescale(-2.844141069) transform(1.tankrabatt_1#2.hhi_km_quantile 1.tankrabatt_1#3.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile:= @+scal_tankrabatt_1_per) recast(bar) base(scal_tankrabatt_1_per) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(oligopolistic)) /// 
		 (E_2_first_e10, rescale(-2.844141069) transform(1.tankrabatt_1#2.hhi_km_quantile 1.tankrabatt_1#3.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile:= @+scal_tankrabatt_2_per) recast(bar) base(scal_tankrabatt_2_per) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(other integrated)) /// 
		 (E_3_first_e10, rescale(-2.844141069) transform(1.tankrabatt_1#2.hhi_km_quantile 1.tankrabatt_1#3.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile:= @+scal_tankrabatt_3_per) recast(bar) base(scal_tankrabatt_3_per) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(independent)), /// 
												  vertical /// 
												  keep(1.tankrabatt_1 1.tankrabatt_1#2.hhi_km_quantile 1.tankrabatt_1#3.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile) /// 
												  legend(pos(6) col(3)) xlabel(1 "({it:Ref.}) Q1" 2 "Q2" 3 "Q3" 4 "Q4") /// 
												  ytitle("Tax Pass-Through E10 (%)") xtitle("HHI Quartiles") /// 
												  addplot(function first=scal_tankrabatt_1_per, ra(0.5 4.5) lwidth(medthick) lcolor(stc1) lpattern(shortdash) || /// 
												  function second=scal_tankrabatt_2_per, ra(0.5 4.5) lwidth(medthick) lcolor(stc2) lpattern(shortdash) || ///
												  function thrid=scal_tankrabatt_3_per, ra(0.5 4.5) lwidth(medthick) lcolor(stc3) lpattern(shortdash)) /// 
												  scale(1.3) ysize(5cm) xsize(8cm) name(e10_brand_hhi_first_per, replace)
	graph save "HHI_VERT_first_e10_per", replace
	graph export "HHI_VERT_first_e10_per.emf", replace
	
*##################### SENSITIVITY (HHI) - MARKET BOUNDARIES #################*
est clear
	* 3km RADIUS
	drop hhi_km_quantile
	xtile hhi_km_quantile = hhi_3km, nquantiles(4)
	capture program drop _all
		program first_e10_hhi_rad_3 
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_hhi_rad_3_file, replace): first_e10_hhi_rad_3 
	ereturn list
	estimates store E3_first_e10_hhi_rad_3
	estimates save E3_first_e10_hhi_rad_3, replace
	scalar scal_tankrabatt_3 = _b[1.tankrabatt_1]
	scalar scal_tankrabatt_3_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_e10
	
	* 4km RADIUS
	estimates use "E4_first_e10_hhi_rad_4"
	eststo: estimates replay "E4_first_e10_hhi_rad_4"
	estimates store E4_first_e10_hhi_rad_4
	scalar scal_tankrabatt_4 = _b[1.tankrabatt_1]
	scalar scal_tankrabatt_4_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_e10
	
	* 5km RADIUS
	drop hhi_km_quantile
	xtile hhi_km_quantile = hhi_5km, nquantiles(4)
	capture program drop _all
		program first_e10_hhi_rad_5 
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_hhi_rad_5_file, replace): first_e10_hhi_rad_5 
	ereturn list
	estimates store e10_first_e10_hhi_rad_5
	estimates save e10_first_e10_hhi_rad_5, replace
	scalar scal_tankrabatt_5 = _b[1.tankrabatt_1]
	scalar scal_tankrabatt_5_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_e10
		esttab using "HHI_SENS_first_boot.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Heterogeneous Effects First Intervention) nomtitles addnotes(Notes follow)
	
coefplot (E3_first_e10_hhi_rad_3, transform(1.tankrabatt_1#2.hhi_km_quantile 1.tankrabatt_1#3.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile:= @+scal_tankrabatt_3) recast(bar) base(scal_tankrabatt_3) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(3km Radius)) /// 
		 (E4_first_e10_hhi_rad_4, transform(1.tankrabatt_1#2.hhi_km_quantile 1.tankrabatt_1#3.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile:= @+scal_tankrabatt_4) recast(bar) base(scal_tankrabatt_4) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(4km Radius)) /// 
		 (e10_first_e10_hhi_rad_5, transform(1.tankrabatt_1#2.hhi_km_quantile 1.tankrabatt_1#3.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile:= @+scal_tankrabatt_5) recast(bar) base(scal_tankrabatt_5) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(5km Radius)), /// 
												  vertical /// 
												  keep(1.tankrabatt_1 1.tankrabatt_1#2.hhi_km_quantile 1.tankrabatt_1#3.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile) /// 
												  legend(pos(6) col(3)) xlabel(1 "({it:Ref.}) Q1" 2 "Q2" 3 "Q3" 4 "Q4") /// 
												  ytitle("Tax Pass-Through E10 (ct/l)") xtitle("HHI Quartiles") /// 
												  addplot(function first=scal_tankrabatt_3, ra(0.5 4.5) lwidth(medthick) lcolor(stc1) lpattern(shortdash) || /// 
												  function second=scal_tankrabatt_4, ra(0.5 4.5) lwidth(medthick) lcolor(stc2) lpattern(shortdash) || ///
												  function thrid=scal_tankrabatt_5, ra(0.5 4.5) lwidth(medthick) lcolor(stc3) lpattern(shortdash)) /// 
												  scale(1.3) ysize(5cm) xsize(8cm) name(e10_sens_hhi_first, replace)
	graph save "HHI_SENS_first_e10", replace
coefplot (E3_first_e10_hhi_rad_3,  rescale(-2.844141069) transform(1.tankrabatt_1#2.hhi_km_quantile 1.tankrabatt_1#3.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile:= @+scal_tankrabatt_3_per) recast(bar) base(scal_tankrabatt_3_per) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(3km Radius)) /// 
		 (E4_first_e10_hhi_rad_4, rescale(-2.844141069) transform(1.tankrabatt_1#2.hhi_km_quantile 1.tankrabatt_1#3.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile:= @+scal_tankrabatt_4_per) recast(bar) base(scal_tankrabatt_4_per) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(4km Radius)) /// 
		 (e10_first_e10_hhi_rad_5, rescale(-2.844141069) transform(1.tankrabatt_1#2.hhi_km_quantile 1.tankrabatt_1#3.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile:= @+scal_tankrabatt_5_per) recast(bar) base(scal_tankrabatt_5_per) barwidth(0.2) lwidth(medthick) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop label(5km Radius)), /// 
												  vertical /// 
												  keep(1.tankrabatt_1 1.tankrabatt_1#2.hhi_km_quantile 1.tankrabatt_1#3.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile) /// 
												  legend(pos(6) col(3)) xlabel(1 "({it:Ref.}) Q1" 2 "Q2" 3 "Q3" 4 "Q4") /// 
												  ytitle("Tax Pass-Through E10 (%)") xtitle("HHI Quartiles") /// 
												  addplot(function first=scal_tankrabatt_3_per, ra(0.5 4.5) lwidth(medthick) lcolor(stc1) lpattern(shortdash) || /// 
												  function second=scal_tankrabatt_4_per, ra(0.5 4.5) lwidth(medthick) lcolor(stc2) lpattern(shortdash) || ///
												  function thrid=scal_tankrabatt_5_per, ra(0.5 4.5) lwidth(medthick) lcolor(stc3) lpattern(shortdash)) /// 
												  scale(1.3) ysize(5cm) xsize(8cm) name(e10_sens_hhi_first_per, replace)
	graph save "HHI_SENS_first_e10_per", replace	
	graph export "HHI_SENS_first_e10_per.emf", replace

* #################### ROBUSTNESS - BANDWIDTH & DONUT #########################*
est clear
	* 240 hour (10 day) bandwidth
	scalar h = 240 // Bandwidth 240 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program first_10days
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_10days, replace): first_10days
	estimates save first_10days, replace
	* 288 hour (12 day) bandwidth
	scalar h = 288 // Bandwidth 288 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program first_12days
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_12days, replace): first_12days
	estimates save first_12days, replace
	* 336 hour (14 day) bandwidth (no donut)
	scalar da = 0 // Donut 0 hours
	scalar db = -da
	scalar h = 336 // Bandwidth 336 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program first_14days_0d
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_14days_0d, replace): first_14days_0d
	estimates save first_14days_0d, replace
	* 336 hour bandwidth (9 hours donut) siehe Average Effects
	estimates use "first_e10_allcov"
	eststo: estimates replay "first_e10_allcov"
	* 336 hour (14 day) bandwidth (18 hours donut)
	scalar da = 18 // Donut 18 hours
	scalar db = -da
	scalar h = 336 // Bandwidth 336 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program first_14days_18d
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_14days_18d, replace): first_14days_18d
	estimates save first_14days_18d, replace
	* 384 hour (16 day) bandwidth
	scalar da = 9 // Donut 9 hours
	scalar db = -da
	scalar h = 384 // Bandwidth 384 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program first_16days
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_16days, replace): first_16days
	estimates save first_16days, replace
		esttab using "ROBUST_e10_first_boot.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Heterogeneous Effects First Intervention) nomtitles addnotes(Notes follow)

*######################### BRAND AFFILIATION ##################################*
	scalar h = 336 // Bandwidth 336 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
est clear	
	* ARAL
	capture program drop _all
		program first_e10_ARAL
			reghdfe mean_e10 $var_sup $var_dem $var_seas if brand == "ARAL", absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if brand == "ARAL" & (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_ARAL, replace): first_e10_ARAL
	estimates save first_e10_ARAL , replace
	* SHELL
	capture program drop _all
		program first_e10_SHELL
			reghdfe mean_e10 $var_sup $var_dem $var_seas if brand == "SHELL", absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if brand == "SHELL" & (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_SHELL, replace): first_e10_SHELL
	estimates save first_e10_SHELL , replace
	* TOTALENERGIES
	capture program drop _all
		program first_e10_TOTALEnergies
			reghdfe mean_e10 $var_sup $var_dem $var_seas if brand == "TOTALEnergies", absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if brand == "TOTALEnergies" & (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_TOTALEnergies, replace): first_e10_TOTALEnergies
	estimates save first_e10_TOTALEnergies , replace
	* ESSO
	capture program drop _all
		program first_e10_ESSO
			reghdfe mean_e10 $var_sup $var_dem $var_seas if brand == "ESSO", absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if brand == "ESSO" & (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_ESSO, replace): first_e10_ESSO
	estimates save first_e10_ESSO , replace
	* JET
	capture program drop _all
		program first_e10_JET
			reghdfe mean_e10 $var_sup $var_dem $var_seas if brand == "JET", absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if brand == "JET" & (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_JET, replace): first_e10_JET
	estimates save first_e10_JET , replace
	* BFT
	capture program drop _all
		program first_e10_BFT
			reghdfe mean_e10 $var_sup $var_dem $var_seas if brand == "BFT", absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if brand == "BFT" & (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_BFT, replace): first_e10_BFT
	estimates save first_e10_BFT , replace
	* AVIA
	capture program drop _all
		program first_e10_AVIA
			reghdfe mean_e10 $var_sup $var_dem $var_seas if brand == "AVIA", absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if brand == "AVIA" & (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_AVIA, replace): first_e10_AVIA
	estimates save first_e10_AVIA , replace
	esttab using "BRAND_e10_first_boot.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Main specification e10 First Intervention) nomtitles addnotes(Notes follow)
		
*#################### EDGEWORTH CYCLES (no hour FE) ###########################*
drop hhi_km_quantile
xtile hhi_km_quantile = hhi_4km, nquantiles(4)
est clear	
	* Average Effect
	capture program drop _all
		program first_e10_allcov_nohour
			reghdfe mean_e10 $var_sup $var_dem i.week_days, absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
									if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_allcov_file, replace): first_e10_allcov_nohour
	estimates save first_e10_nohour , replace
	* HHI Quartiles (4km):
	capture program drop _all
		program first_e10_hhi_4km_nohour
			reghdfe mean_e10 $var_sup $var_dem i.week_days, absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_hhi_4km_file, replace): first_e10_hhi_4km_nohour
	estimates save first_e10_nohour_hhi , replace
	* Average Effect: Morning Peak (05:00 to 10:59)
	capture program drop _all
		program first_e10_nohour_5_10
			reghdfe mean_e10 $var_sup $var_dem i.week_days, absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h))) & inrange(hour,5,10)
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_nohour_5_10_file, replace): first_e10_nohour_5_10
	estimates save first_e10_nohour_5_10, replace
	* Average Effect: Mid Day (11:00 to 16:59)
	capture program drop _all
		program first_e10_nohour_11_16
			reghdfe mean_e10 $var_sup $var_dem i.week_days, absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h))) & inrange(hour,11,16)
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_nohour_11_16_file, replace): first_e10_nohour_11_16
	estimates save first_e10_nohour_11_16, replace
	* Average Effect: Evening Off-Peak (17:00 to 22:59)
	capture program drop _all
		program first_e10_nohour_17_22
			reghdfe mean_e10 $var_sup $var_dem i.week_days, absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h))) & inrange(hour,17,22)
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_nohour_17_22_file, replace): first_e10_nohour_17_22
	estimates save first_e10_nohour_17_22, replace
	* Average Effect: Night (23:00 to 04:59)
	capture program drop _all
		program first_e10_nohour_23_4
			reghdfe mean_e10 $var_sup $var_dem i.week_days, absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h))) & (inrange(hour,0,4) | hour == 23)
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_nohour_23_4_file, replace): first_e10_nohour_23_4
	estimates save first_e10_nohour_23_4, replace
	* HHI Quartiles (4km): Morning Peak (05:00 to 10:59)
	capture program drop _all
		program first_e10_nohour_5_10_hhi
			reghdfe mean_e10 $var_sup $var_dem i.week_days, absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h))) & inrange(hour,5,10)
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_nohour_5_10_file_hhi, replace): first_e10_nohour_5_10_hhi
	estimates save first_e10_nohour_5_10_hhi, replace
	* HHI Quartiles (4km): Mid Day (11:00 to 16:59)
	capture program drop _all
		program first_e10_nohour_11_16_hhi
			reghdfe mean_e10 $var_sup $var_dem i.week_days, absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h))) & inrange(hour,11,16)
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_nohour_11_16_file_hhi, replace): first_e10_nohour_11_16_hhi
	estimates save first_e10_nohour_11_16_hhi, replace
	* HHI Quartiles (4km): Evening Off-Peak (17:00 to 22:59)
	capture program drop _all
		program first_e10_nohour_17_22_hhi
			reghdfe mean_e10 $var_sup $var_dem i.week_days, absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h))) & inrange(hour,17,22)
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_nohour_17_22_file_hhi, replace): first_e10_nohour_17_22_hhi
	estimates save first_e10_nohour_17_22_hhi, replace
	* HHI Quartiles (4km): Night (23:00 to 04:59)
	capture program drop _all
		program first_e10_nohour_23_4_hhi
			reghdfe mean_e10 $var_sup $var_dem i.week_days, absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h))) & (inrange(hour,0,4) | hour == 23)
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_nohour_23_4_file_hhi, replace): first_e10_nohour_23_4_hhi
	estimates save first_e10_nohour_23_4_hhi, replace
	esttab using "nohour_e10_first_boot.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) keep(1.tankrabatt_1 1.tankrabatt_1#2.hhi_km_quantile 1.tankrabatt_1#3.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile _cons) title(NO hourly FE e10 First Intervention) nomtitles addnotes(Notes follow)

*################### HHI and VERT-INT (one regression) ########################*
est clear	
	* Average Effect
	capture program drop _all
		program first_e10_het_included
			reghdfe mean_e10 $var_sup $var_dem i.week_days, absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.hhi_km_quantile##i.brand_vert_int [aw=kernel_wgt] ///
									if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_het_included, replace): first_e10_het_included
	estimates save first_e10_het_included , replace
	esttab using "hhi_vert_tog_e10_first_boot.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) /// 
	keep(1.tankrabatt_1 /// 
	     2.brand_vert_int 3.brand_vert_int /// 
		 2.hhi_km_quantile 3.hhi_km_quantile 4.hhi_km_quantile ///
	     1.tankrabatt_1#2.hhi_km_quantile 1.tankrabatt_1#3.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile /// 
		 1.tankrabatt_1#2.brand_vert_int 1.tankrabatt_1#3.brand_vert_int /// 
		 1.tankrabatt_1#2.hhi_km_quantile#2.brand_vert_int 1.tankrabatt_1#2.hhi_km_quantile#3.brand_vert_int /// 
		 1.tankrabatt_1#3.hhi_km_quantile#2.brand_vert_int 1.tankrabatt_1#3.hhi_km_quantile#3.brand_vert_int ///
		 1.tankrabatt_1#4.hhi_km_quantile#2.brand_vert_int 1.tankrabatt_1#4.hhi_km_quantile#3.brand_vert_int ///
		 _cons) title(HHI and VERT e10 First Intervention) nomtitles addnotes(Notes follow)

			
timer off 1
timer list
log close