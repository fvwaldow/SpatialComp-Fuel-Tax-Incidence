* Bootstrap Iterations:
global boot_rep reps(1000)
*-------------------------------------------------------------------------------

*###############################################################################
* ROADSIDE STATIONS
*###############################################################################

*-------------------------------------------------------------------------------
* FIRST INTERVENTION
*-------------------------------------------------------------------------------
* DIESEL
cd "C:\Users\frede\Desktop\Bachelor-Thesis\Daten\tankstellenpreise"
use "roadside_price_2022_0309_openandassume_revised.dta", clear

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
drop mean_e5 mean_e10 /// 
     d_gborder_10_BEL d_gborder_10_CHE d_gborder_10_CZE d_gborder_10_DNK d_gborder_10_FRA d_gborder_10_LUX d_gborder_10_NLD d_gborder_10_POL d_gborder_10_AUT /// 
	 own_share_3km totalSize_4km own_share_4km own_share_5km /// 
	 month land station_uuid NUTS

xtile hhi_km_quantile = hhi_4km, nquantiles(4)

*################### HHI and VERT-INT (one regression) ########################*
cd "C:\Users\frede\Desktop\RDiT JTEP Results\full interaction HHI VERT"
est clear	
est clear
gen d_comp =0
replace d_comp = 1 if inrange(d_next_competitor,1,3)
replace d_comp = 2 if d_next_competitor >3
	* Region Type & Vert Int:
	capture program drop _all
		program first_diesel_reg_typ
			reghdfe mean_diesel $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.brand_vert_int /// 
			                1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.urban_rural_cat [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_diesel 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_diesel_reg_typ_file, replace): first_diesel_reg_typ
	estimates save first_diesel_reg_typ , replace
	* Distance next competitor & Vert Int:
	capture program drop _all
		program first_diesel_dist_comp
			reghdfe mean_diesel $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.brand_vert_int /// 
			                1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.d_comp [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_diesel 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_diesel_dist_comp_file, replace): first_diesel_dist_comp
	estimates save first_diesel_dist_comp , replace
	* HHI Quartiles (4km) & Vert Int:
		capture program drop _all
		program first_diesel_hhi_4km
			reghdfe mean_diesel $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.brand_vert_int /// 
			                1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_diesel 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_diesel_hhi_4km_file, replace): first_diesel_hhi_4km
	estimates save E4_first_diesel_hhi_rad_4 , replace
	esttab using "Het_inter_diesel_first_boot.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Heterogeneous Effects First Intervention diesel) keep(1.tankrabatt_1 1.tankrabatt_1#2.brand_vert_int 1.tankrabatt_1#3.brand_vert_int 1.tankrabatt_1#2.urban_rural_cat 1.tankrabatt_1#3.urban_rural_cat 1.tankrabatt_1#2.hhi_km_quantile 1.tankrabatt_1#3.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile 1.tankrabatt_1#1.d_comp 1.tankrabatt_1#2.d_comp _cons) nomtitles addnotes(Notes follow)

*###############################################################################
*E5
clear all
cd "C:\Users\frede\Desktop\Bachelor-Thesis\Daten\tankstellenpreise"
use "roadside_price_2022_0309_openandassume_revised.dta", clear

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
drop mean_diesel mean_e10 /// 
     d_gborder_10_BEL d_gborder_10_CHE d_gborder_10_CZE d_gborder_10_DNK d_gborder_10_FRA d_gborder_10_LUX d_gborder_10_NLD d_gborder_10_POL d_gborder_10_AUT /// 
	 own_share_3km totalSize_4km own_share_4km own_share_5km /// 
	 month land station_uuid NUTS

xtile hhi_km_quantile = hhi_4km, nquantiles(4)	 
	 
*################### HHI and VERT-INT (one regression) ########################*
cd "C:\Users\frede\Desktop\RDiT JTEP Results\full interaction HHI VERT"
est clear	
est clear
gen d_comp =0
replace d_comp = 1 if inrange(d_next_competitor,1,3)
replace d_comp = 2 if d_next_competitor >3
	* Region Type & Vert Int:
	capture program drop _all
		program first_e5_reg_typ
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.brand_vert_int /// 
			                1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.urban_rural_cat [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e5_reg_typ_file, replace): first_e5_reg_typ
	estimates save first_e5_reg_typ , replace
	* Distance next competitor & Vert Int:
	capture program drop _all
		program first_e5_dist_comp
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.brand_vert_int /// 
			                1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.d_comp [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e5_dist_comp_file, replace): first_e5_dist_comp
	estimates save first_e5_dist_comp , replace
	* HHI Quartiles (4km) & Vert Int:
		capture program drop _all
		program first_e5_hhi_4km
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.brand_vert_int /// 
			                1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e5_hhi_4km_file, replace): first_e5_hhi_4km
	estimates save E4_first_e5_hhi_rad_4 , replace
	esttab using "Het_inter_e5_first_boot.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Heterogeneous Effects First Intervention e5) keep(1.tankrabatt_1 1.tankrabatt_1#2.brand_vert_int 1.tankrabatt_1#3.brand_vert_int 1.tankrabatt_1#2.urban_rural_cat 1.tankrabatt_1#3.urban_rural_cat 1.tankrabatt_1#2.hhi_km_quantile 1.tankrabatt_1#3.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile 1.tankrabatt_1#1.d_comp 1.tankrabatt_1#2.d_comp _cons) nomtitles addnotes(Notes follow)

*###############################################################################
*E10
clear all
cd "C:\Users\frede\Desktop\Bachelor-Thesis\Daten\tankstellenpreise"
use "roadside_price_2022_0309_openandassume_revised.dta", clear

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

xtile hhi_km_quantile = hhi_4km, nquantiles(4)	 
	 
*################### HHI and VERT-INT (one regression) ########################*
cd "C:\Users\frede\Desktop\RDiT JTEP Results\full interaction HHI VERT"
est clear	
est clear
gen d_comp =0
replace d_comp = 1 if inrange(d_next_competitor,1,3)
replace d_comp = 2 if d_next_competitor >3
	* Region Type & Vert Int:
	capture program drop _all
		program first_e10_reg_typ
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.brand_vert_int /// 
			                1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.urban_rural_cat [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_reg_typ_file, replace): first_e10_reg_typ
	estimates save first_e10_reg_typ , replace
	* Distance next competitor & Vert Int:
	capture program drop _all
		program first_e10_dist_comp
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.brand_vert_int /// 
			                1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.d_comp [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_dist_comp_file, replace): first_e10_dist_comp
	estimates save first_e10_dist_comp , replace
	* HHI Quartiles (4km) & Vert Int:
		capture program drop _all
		program first_e10_hhi_4km
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.brand_vert_int /// 
			                1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_hhi_4km_file, replace): first_e10_hhi_4km
	estimates save E4_first_e10_hhi_rad_4 , replace
	esttab using "Het_inter_e10_first_boot.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Heterogeneous Effects First Intervention e10) keep(1.tankrabatt_1 1.tankrabatt_1#2.brand_vert_int 1.tankrabatt_1#3.brand_vert_int 1.tankrabatt_1#2.urban_rural_cat 1.tankrabatt_1#3.urban_rural_cat 1.tankrabatt_1#2.hhi_km_quantile 1.tankrabatt_1#3.hhi_km_quantile 1.tankrabatt_1#4.hhi_km_quantile 1.tankrabatt_1#1.d_comp 1.tankrabatt_1#2.d_comp _cons) nomtitles addnotes(Notes follow)
		 
*-------------------------------------------------------------------------------
* SECOND INTERVENTION
*-------------------------------------------------------------------------------
clear all
cd "C:\Users\frede\Desktop\Bachelor-Thesis\Daten\tankstellenpreise"
use "roadside_price_2022_0309_openandassume_revised.dta", clear

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
drop mean_e5 mean_e10 /// 
     d_gborder_10_BEL d_gborder_10_CHE d_gborder_10_CZE d_gborder_10_DNK d_gborder_10_FRA d_gborder_10_LUX d_gborder_10_NLD d_gborder_10_POL d_gborder_10_AUT /// 
	 own_share_3km totalSize_4km own_share_4km own_share_5km /// 
	 month land station_uuid NUTS

xtile hhi_km_quantile = hhi_4km, nquantiles(4)	 
	 
*################### HHI and VERT-INT (one regression) ########################*
cd "C:\Users\frede\Desktop\RDiT JTEP Results\full interaction HHI VERT"
est clear	
est clear
gen d_comp =0
replace d_comp = 1 if inrange(d_next_competitor,1,3)
replace d_comp = 2 if d_next_competitor >3
	* Region Type & Vert Int:
	capture program drop _all
		program second_diesel_reg_typ
			reghdfe mean_diesel $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.brand_vert_int /// 
			                1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.urban_rural_cat [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_diesel 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_diesel_reg_typ_file, replace): second_diesel_reg_typ
	estimates save second_diesel_reg_typ , replace
	* Distance next competitor & Vert Int:
	capture program drop _all
		program second_diesel_dist_comp
			reghdfe mean_diesel $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.brand_vert_int /// 
			                1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.d_comp [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_diesel 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_diesel_dist_comp_file, replace): second_diesel_dist_comp
	estimates save second_diesel_dist_comp , replace
	* HHI Quartiles (4km) & Vert Int:
		capture program drop _all
		program second_diesel_hhi_4km
			reghdfe mean_diesel $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.brand_vert_int /// 
			                1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_diesel 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_diesel_hhi_4km_file, replace): second_diesel_hhi_4km
	estimates save E4_second_diesel_hhi_rad_4 , replace
	esttab using "Het_inter_diesel_second_boot.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Heterogeneous Effects second Intervention diesel) keep(1.tankrabatt_2 1.tankrabatt_2#2.brand_vert_int 1.tankrabatt_2#3.brand_vert_int 1.tankrabatt_2#2.urban_rural_cat 1.tankrabatt_2#3.urban_rural_cat 1.tankrabatt_2#2.hhi_km_quantile 1.tankrabatt_2#3.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile 1.tankrabatt_2#1.d_comp 1.tankrabatt_2#2.d_comp _cons) nomtitles addnotes(Notes follow)
		 
*###############################################################################
*E5
clear all
cd "C:\Users\frede\Desktop\Bachelor-Thesis\Daten\tankstellenpreise"
use "roadside_price_2022_0309_openandassume_revised.dta", clear

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

xtile hhi_km_quantile = hhi_4km, nquantiles(4)	 
	 
*################### HHI and VERT-INT (one regression) ########################*
cd "C:\Users\frede\Desktop\RDiT JTEP Results\full interaction HHI VERT"
est clear	
est clear
gen d_comp =0
replace d_comp = 1 if inrange(d_next_competitor,1,3)
replace d_comp = 2 if d_next_competitor >3
	* Region Type & Vert Int:
	capture program drop _all
		program second_e5_reg_typ
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.brand_vert_int /// 
			                1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.urban_rural_cat [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_reg_typ_file, replace): second_e5_reg_typ
	estimates save second_e5_reg_typ , replace
	* Distance next competitor & Vert Int:
	capture program drop _all
		program second_e5_dist_comp
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.brand_vert_int /// 
			                1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.d_comp [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_dist_comp_file, replace): second_e5_dist_comp
	estimates save second_e5_dist_comp , replace
	* HHI Quartiles (4km) & Vert Int:
		capture program drop _all
		program second_e5_hhi_4km
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.brand_vert_int /// 
			                1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_hhi_4km_file, replace): second_e5_hhi_4km
	estimates save E4_second_e5_hhi_rad_4 , replace
	esttab using "Het_inter_e5_second_boot.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Heterogeneous Effects second Intervention e5) keep(1.tankrabatt_2 1.tankrabatt_2#2.brand_vert_int 1.tankrabatt_2#3.brand_vert_int 1.tankrabatt_2#2.urban_rural_cat 1.tankrabatt_2#3.urban_rural_cat 1.tankrabatt_2#2.hhi_km_quantile 1.tankrabatt_2#3.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile 1.tankrabatt_2#1.d_comp 1.tankrabatt_2#2.d_comp _cons) nomtitles addnotes(Notes follow)

*###############################################################################
*E10
clear all
cd "C:\Users\frede\Desktop\Bachelor-Thesis\Daten\tankstellenpreise"
use "roadside_price_2022_0309_openandassume_revised.dta", clear

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
drop mean_e5 mean_diesel /// 
     d_gborder_10_BEL d_gborder_10_CHE d_gborder_10_CZE d_gborder_10_DNK d_gborder_10_FRA d_gborder_10_LUX d_gborder_10_NLD d_gborder_10_POL d_gborder_10_AUT /// 
	 own_share_3km totalSize_4km own_share_4km own_share_5km /// 
	 month land station_uuid NUTS

xtile hhi_km_quantile = hhi_4km, nquantiles(4)	 
	 
*################### HHI and VERT-INT (one regression) ########################*
cd "C:\Users\frede\Desktop\RDiT JTEP Results\full interaction HHI VERT"
est clear	
est clear
gen d_comp =0
replace d_comp = 1 if inrange(d_next_competitor,1,3)
replace d_comp = 2 if d_next_competitor >3
	* Region Type & Vert Int:
	capture program drop _all
		program second_e10_reg_typ
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.brand_vert_int /// 
			                1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.urban_rural_cat [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e10_reg_typ_file, replace): second_e10_reg_typ
	estimates save second_e10_reg_typ , replace
	* Distance next competitor & Vert Int:
	capture program drop _all
		program second_e10_dist_comp
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.brand_vert_int /// 
			                1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.d_comp [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e10_dist_comp_file, replace): second_e10_dist_comp
	estimates save second_e10_dist_comp , replace
	* HHI Quartiles (4km) & Vert Int:
		capture program drop _all
		program second_e10_hhi_4km
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.brand_vert_int /// 
			                1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2)##i.hhi_km_quantile [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e10 
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e10_hhi_4km_file, replace): second_e10_hhi_4km
	estimates save E4_second_e10_hhi_rad_4 , replace
	esttab using "Het_inter_e10_second_boot.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Heterogeneous Effects second Intervention e10) keep(1.tankrabatt_2 1.tankrabatt_2#2.brand_vert_int 1.tankrabatt_2#3.brand_vert_int 1.tankrabatt_2#2.urban_rural_cat 1.tankrabatt_2#3.urban_rural_cat 1.tankrabatt_2#2.hhi_km_quantile 1.tankrabatt_2#3.hhi_km_quantile 1.tankrabatt_2#4.hhi_km_quantile 1.tankrabatt_2#1.d_comp 1.tankrabatt_2#2.d_comp _cons) nomtitles addnotes(Notes follow)