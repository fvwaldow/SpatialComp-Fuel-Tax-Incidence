*###############################################################################
* DO-FILE (HIGHWAY)
* Spatial Competition and Fuel Tax Pass-Through
* Hourly Prices & Two-Stage RDiT Approach (incl. Bootstrapping) - Revision
* Frederik von Waldow
* 15.04.2025
*###############################################################################
clear all
ssc install estout, replace
ssc install outreg2
ssc install reghdfe
ssc install ftools
ssc install coefplot, replace
net from http://www.stata-journal.com/software/sj3-2/
net describe st0039
net install st0039
graph set window fontface "Times New Roman"

*-------------------------------------------------------------------------------
*###############################################################################
* HIGHWAY STATIONS
*###############################################################################
cd "C:\Users\frede\Desktop\Code_TwoStage_RDiT_vWaldow\HIGH"
log using "BOOST_HIGHWAY", replace

*-------------------------------------------------------------------------------
* Bootstrap Iterations:
global boot_rep reps(1000)
*-------------------------------------------------------------------------------

*-------------------------------------------------------------------------------
* FIRST INTERVENTION
*-------------------------------------------------------------------------------
cd "C:\Users\frede\Desktop\Bachelor-Thesis\Daten\tankstellenpreise"
use "roadside_price_2022_0309_openandassume_high_revised.dta", clear
cd "C:\Users\frede\Desktop\Code_TwoStage_RDiT_vWaldow\HIGH"
gen id_new = station_id
xtset station_id time
keep if inrange(time,9961,13321) // 20 Weeks (10 Weeks pre treatment, 10 Week with treatment)
bysort station_id (time): gen L24_BRENT_FOB = L24.BRENT_FOB
bysort station_id (time): gen L24_usd_to_euro = L24.usd_to_euro
bysort station_id (time): gen L1_kfz_menge_A = L1.kfz_menge_A

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
global var_dem kfz_menge_A L1_kfz_menge_A temperature precipitation holiday_state day_before_holiday_state // demand covariates 
/*
* RD PLOTS:
foreach i in mean_e5 mean_e10 mean_diesel {
	reghdfe `i' $var_seas, absorb(station_id) resid(res_`i'_seas)
	rdplot res_`i'_seas hour_to_treat_1 if  (inrange(hour_to_treat_1, -120,-9) | inrange(hour_to_treat_1,9,120)), c(0) h(120) p(1) nbins(28 28) masspoints(adjust) graph_options(ytitle(Residuals `i' (ct/l)) xtitle(Hours to Treatment) legend(off) scale(1.3) name(`i'_first, replace))
	graph save `i'_rdplot_first, replace
	drop res_`i'_seas
}	
*/
timer on 1
*-------------------------------------------------------------------------------
*############################### E5 ###########################################*
*-------------------------------------------------------------------------------
*######################### AVERAGE EFFECTS ####################################*
	* Seasonality
est clear
	capture program drop _all
		program first_e5_season_high
			reghdfe mean_e5 $var_seas , absorb(station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e5
		end 	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e5_season_file_high, replace): first_e5_season_high
	estimates save first_e5_season_high , replace
	* Supply & Seasonality
		program first_e5_sup_high
			reghdfe mean_e5 $var_sup $var_seas , absorb(station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e5
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new) saving(first_e5_sup_file_high, replace): first_e5_sup_high
	estimates save first_e5_sup_high , replace
	* Demand & Seasonality
		program first_e5_dem_high
			reghdfe mean_e5 $var_dem $var_seas , absorb(station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e5
		end
		xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new) saving(first_e5_sem_file_high, replace): first_e5_dem_high
	estimates save first_e5_dem_high , replace
	* Supply & Demand & Seasonality
	capture program drop _all
		program first_e5_allcov_high
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e5 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e5_allcov_file_high, replace): first_e5_allcov_high
	estimates save first_e5_allcov_high , replace
	esttab using "Main_e5_first_boot_high.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Main specification e5 First Intervention) nomtitles addnotes(Notes follow)

* #################### ROBUSTNESS - BANDWIDTH & DONUT #########################*
est clear
	* 240 hour (10 day) bandwidth
	scalar h = 240 // Bandwidth 240 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program first_10days_high
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e5 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_10days_high, replace): first_10days_high
	estimates save first_10days_high, replace
	* 288 hour (12 day) bandwidth
	scalar h = 288 // Bandwidth 288 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program first_12days_high
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e5 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_12days_high, replace): first_12days_high
	estimates save first_12days_high, replace
	* 336 hour (14 day) bandwidth (no donut)
	scalar da = 0 // Donut 0 hours
	scalar db = -da
	scalar h = 336 // Bandwidth 336 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program first_14days_0d_high
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e5 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_14days_0d_high, replace): first_14days_0d_high
	estimates save first_14days_0d_high, replace
	* 336 hour bandwidth (9 hours donut) siehe Average Effects
	estimates use "first_e5_allcov_high"
	eststo: estimates replay "first_e5_allcov_high"
	* 336 hour (14 day) bandwidth (18 hours donut)
	scalar da = 18 // Donut 18 hours
	scalar db = -da
	scalar h = 336 // Bandwidth 336 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program first_14days_18d_high
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e5 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_14days_18d_high, replace): first_14days_18d_high
	estimates save first_14days_18d_high, replace
	* 384 hour (16 day) bandwidth
	scalar da = 9 // Donut 9 hours
	scalar db = -da
	scalar h = 384 // Bandwidth 384 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program first_16days_high
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e5 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_16days_high, replace): first_16days_high
	estimates save first_16days_high, replace
		esttab using "ROBUST_e5_first_boot_high.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Heterogeneous Effects First Intervention) nomtitles addnotes(Notes follow)

*-------------------------------------------------------------------------------
*############################### E10 ##########################################*
*-------------------------------------------------------------------------------
*######################### AVERAGE EFFECTS ####################################*
scalar da = 9 // Donut 9 hours (baseline)
scalar db = -da
scalar h = 336 // Bandwidth 336 hours (baseline)
drop kernel_wgt
gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h)) // triangular kernel (non-parametric estimation)
	replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
	replace kernel_wgt = 0 if mi(kernel_wgt)
	* Seasonality
est clear
	capture program drop _all
		program first_e10_season_high
			reghdfe mean_e10 $var_seas , absorb(station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10
		end 	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_season_file_high, replace): first_e10_season_high
	estimates save first_e10_season_high , replace
	* Supply & Seasonality
		program first_e10_sup_high
			reghdfe mean_e10 $var_sup $var_seas , absorb(station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new) saving(first_e10_sup_file_high, replace): first_e10_sup_high
	estimates save first_e10_sup_high , replace
	* Demand & Seasonality
		program first_e10_dem_high
			reghdfe mean_e10 $var_dem $var_seas , absorb(station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10
		end
		xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new) saving(first_e10_sem_file_high, replace): first_e10_dem_high
	estimates save first_e10_dem_high , replace
	* Supply & Demand & Seasonality
	capture program drop _all
		program first_e10_allcov_high
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_e10_allcov_file_high, replace): first_e10_allcov_high
	estimates save first_e10_allcov_high , replace
	esttab using "Main_e10_first_boot_high.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Main specification e10 First Intervention) nomtitles addnotes(Notes follow)

* #################### ROBUSTNESS - BANDWIDTH & DONUT #########################*
est clear
	* 240 hour (10 day) bandwidth
	scalar h = 240 // Bandwidth 240 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program first_10days_high
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_10days_high, replace): first_10days_high
	estimates save first_10days_high, replace
	* 288 hour (12 day) bandwidth
	scalar h = 288 // Bandwidth 288 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program first_12days_high
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_12days_high, replace): first_12days_high
	estimates save first_12days_high, replace
	* 336 hour (14 day) bandwidth (no donut)
	scalar da = 0 // Donut 0 hours
	scalar db = -da
	scalar h = 336 // Bandwidth 336 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program first_14days_0d_high
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_14days_0d_high, replace): first_14days_0d_high
	estimates save first_14days_0d_high, replace
	* 336 hour bandwidth (9 hours donut) siehe Average Effects
	estimates use "first_e10_allcov_high"
	eststo: estimates replay "first_e10_allcov_high"
	* 336 hour (14 day) bandwidth (18 hours donut)
	scalar da = 18 // Donut 18 hours
	scalar db = -da
	scalar h = 336 // Bandwidth 336 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program first_14days_18d_high
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_14days_18d_high, replace): first_14days_18d_high
	estimates save first_14days_18d_high, replace
	* 384 hour (16 day) bandwidth
	scalar da = 9 // Donut 9 hours
	scalar db = -da
	scalar h = 384 // Bandwidth 384 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program first_16days_high
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_e10 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_16days_high, replace): first_16days_high
	estimates save first_16days_high, replace
		esttab using "ROBUST_e10_first_boot_high.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Heterogeneous Effects First Intervention) nomtitles addnotes(Notes follow)

*-------------------------------------------------------------------------------
*############################# DIESEL #########################################*
*-------------------------------------------------------------------------------
*######################### AVERAGE EFFECTS ####################################*
scalar da = 9 // Donut 9 hours (baseline)
scalar db = -da
scalar h = 336 // Bandwidth 336 hours (baseline)
drop kernel_wgt
gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h)) // triangular kernel (non-parametric estimation)
	replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
	replace kernel_wgt = 0 if mi(kernel_wgt)
		* Seasonality
est clear
	capture program drop _all
		program first_diesel_season_high
			reghdfe mean_diesel $var_seas , absorb(station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_diesel
		end 	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_diesel_season_file_high, replace): first_diesel_season_high
	estimates save first_diesel_season_high , replace
	* Supply & Seasonality
		program first_diesel_sup_high
			reghdfe mean_diesel $var_sup $var_seas , absorb(station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_diesel
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new) saving(first_diesel_sup_file_high, replace): first_diesel_sup_high
	estimates save first_diesel_sup_high , replace
	* Demand & Seasonality
		program first_diesel_dem_high
			reghdfe mean_diesel $var_dem $var_seas , absorb(station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_diesel
		end
		xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new) saving(first_diesel_sem_file_high, replace): first_diesel_dem_high
	estimates save first_diesel_dem_high , replace
	* Supply & Demand & Seasonality
	capture program drop _all
		program first_diesel_allcov_high
			reghdfe mean_diesel $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_diesel 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_diesel_allcov_file_high, replace): first_diesel_allcov_high
	estimates save first_diesel_allcov_high , replace
	esttab using "Main_diesel_first_boot_high.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Main specification diesel First Intervention) nomtitles addnotes(Notes follow)

* #################### ROBUSTNESS - BANDWIDTH & DONUT #########################*
est clear
	* 240 hour (10 day) bandwidth
	scalar h = 240 // Bandwidth 240 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program first_10days_high
			reghdfe mean_diesel $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_diesel 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_10days_high, replace): first_10days_high
	estimates save first_10days_high, replace
	* 288 hour (12 day) bandwidth
	scalar h = 288 // Bandwidth 288 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program first_12days_high
			reghdfe mean_diesel $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_diesel 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_12days_high, replace): first_12days_high
	estimates save first_12days_high, replace
	* 336 hour (14 day) bandwidth (no donut)
	scalar da = 0 // Donut 0 hours
	scalar db = -da
	scalar h = 336 // Bandwidth 336 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program first_14days_0d_high
			reghdfe mean_diesel $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_diesel 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_14days_0d_high, replace): first_14days_0d_high
	estimates save first_14days_0d_high, replace
	* 336 hour bandwidth (9 hours donut) siehe Average Effects
	estimates use "first_diesel_allcov_high"
	eststo: estimates replay "first_diesel_allcov_high"
	* 336 hour (14 day) bandwidth (18 hours donut)
	scalar da = 18 // Donut 18 hours
	scalar db = -da
	scalar h = 336 // Bandwidth 336 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program first_14days_18d_high
			reghdfe mean_diesel $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_diesel 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_14days_18d_high, replace): first_14days_18d_high
	estimates save first_14days_18d_high, replace
	* 384 hour (16 day) bandwidth
	scalar da = 9 // Donut 9 hours
	scalar db = -da
	scalar h = 384 // Bandwidth 384 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program first_16days_high
			reghdfe mean_diesel $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_1##(c.hour_to_treat_1##c.hour_to_treat_1) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
			drop res_mean_diesel 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(first_16days_high, replace): first_16days_high
	estimates save first_16days_high, replace
		esttab using "ROBUST_diesel_first_boot_high.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Heterogeneous Effects First Intervention) nomtitles addnotes(Notes follow)
		
		
		
*-------------------------------------------------------------------------------
* SECOND INTERVENTION
*-------------------------------------------------------------------------------
cd "C:\Users\frede\Desktop\Bachelor-Thesis\Daten\tankstellenpreise"
use "roadside_price_2022_0309_openandassume_high_revised.dta", clear
cd "C:\Users\frede\Desktop\Code_TwoStage_RDiT_vWaldow\HIGH"
gen id_new = station_id
xtset station_id time
keep if inrange(time,12169,15529) // 20 Weeks (10 Weeks pre treatment, 10 Week with treatment)
bysort station_id (time): gen L24_BRENT_FOB = L24.BRENT_FOB
bysort station_id (time): gen L24_usd_to_euro = L24.usd_to_euro
bysort station_id (time): gen L1_kfz_menge_A = L1.kfz_menge_A

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
global var_dem kfz_menge_A L1_kfz_menge_A temperature precipitation holiday_state day_before_holiday_state // demand covariates
/*
* RD PLOTS:
foreach i in mean_e5 mean_e10 mean_diesel {
	reghdfe `i' $var_seas, absorb(station_id) resid(res_`i'_seas)
	rdplot res_`i'_seas hour_to_treat_2 if  (inrange(hour_to_treat_2, -120,-9) | inrange(hour_to_treat_2,9,120)), c(0) h(120) p(1) nbins(28 28) masspoints(adjust) graph_options(ytitle(Residuals `i' (ct/l)) xtitle(Hours to Treatment) legend(off) scale(1.3) name(`i'_second, replace))
	graph save `i'_rdplot_second, replace
	drop res_`i'_seas
}	
*/

*-------------------------------------------------------------------------------
*############################### E5 ###########################################*
*-------------------------------------------------------------------------------
*######################### AVERAGE EFFECTS ####################################*
	* Seasonality
est clear
	capture program drop _all
		program second_e5_season_high
			reghdfe mean_e5 $var_seas , absorb(station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5
		end 	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_season_file_high, replace): second_e5_season_high
	estimates save second_e5_season_high , replace
	* Supply & Seasonality
		program second_e5_sup_high
			reghdfe mean_e5 $var_sup $var_seas , absorb(station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new) saving(second_e5_sup_file_high, replace): second_e5_sup_high
	estimates save second_e5_sup_high , replace
	* Demand & Seasonality
		program second_e5_dem_high
			reghdfe mean_e5 $var_dem $var_seas , absorb(station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5
		end
		xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new) saving(second_e5_sem_file_high, replace): second_e5_dem_high
	estimates save second_e5_dem_high , replace
	* Supply & Demand & Seasonality
	capture program drop _all
		program second_e5_allcov_high
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e5_allcov_file_high, replace): second_e5_allcov_high
	estimates save second_e5_allcov_high , replace
	esttab using "Main_e5_second_boot_high.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Main specification e5 second Intervention) nomtitles addnotes(Notes follow)

* #################### ROBUSTNESS - BANDWIDTH & DONUT #########################*
est clear
	* 240 hour (10 day) bandwidth
	scalar h = 240 // Bandwidth 240 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program second_10days_high
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_10days_high, replace): second_10days_high
	estimates save second_10days_high, replace
	* 288 hour (12 day) bandwidth
	scalar h = 288 // Bandwidth 288 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program second_12days_high
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_12days_high, replace): second_12days_high
	estimates save second_12days_high, replace
	* 336 hour (14 day) bandwidth (no donut)
	scalar da = 0 // Donut 0 hours
	scalar db = -da
	scalar h = 336 // Bandwidth 336 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program second_14days_0d_high
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_14days_0d_high, replace): second_14days_0d_high
	estimates save second_14days_0d_high, replace
	* 336 hour bandwidth (9 hours donut) siehe Average Effects
	estimates use "second_e5_allcov_high"
	eststo: estimates replay "second_e5_allcov_high"
	* 336 hour (14 day) bandwidth (18 hours donut)
	scalar da = 18 // Donut 18 hours
	scalar db = -da
	scalar h = 336 // Bandwidth 336 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program second_14days_18d_high
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_14days_18d_high, replace): second_14days_18d_high
	estimates save second_14days_18d_high, replace
	* 384 hour (16 day) bandwidth
	scalar da = 9 // Donut 9 hours
	scalar db = -da
	scalar h = 384 // Bandwidth 384 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program second_16days_high
			reghdfe mean_e5 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e5)
			reg res_mean_e5 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e5 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_16days_high, replace): second_16days_high
	estimates save second_16days_high, replace
		esttab using "ROBUST_e5_second_boot_high.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Heterogeneous Effects second Intervention) nomtitles addnotes(Notes follow)

*-------------------------------------------------------------------------------
*############################### E10 ##########################################*
*-------------------------------------------------------------------------------
*######################### AVERAGE EFFECTS ####################################*
scalar da = 9 // Donut 9 hours (baseline)
scalar db = -da
scalar h = 336 // Bandwidth 336 hours (baseline)
drop kernel_wgt
gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h)) // triangular kernel (non-parametric estimation)
	replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
	replace kernel_wgt = 0 if mi(kernel_wgt)
	* Seasonality
est clear
	capture program drop _all
		program second_e10_season_high
			reghdfe mean_e10 $var_seas , absorb(station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e10
		end 	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e10_season_file_high, replace): second_e10_season_high
	estimates save second_e10_season_high , replace
	* Supply & Seasonality
		program second_e10_sup_high
			reghdfe mean_e10 $var_sup $var_seas , absorb(station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e10
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new) saving(second_e10_sup_file_high, replace): second_e10_sup_high
	estimates save second_e10_sup_high , replace
	* Demand & Seasonality
		program second_e10_dem_high
			reghdfe mean_e10 $var_dem $var_seas , absorb(station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e10
		end
		xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new) saving(second_e10_sem_file_high, replace): second_e10_dem_high
	estimates save second_e10_dem_high , replace
	* Supply & Demand & Seasonality
	capture program drop _all
		program second_e10_allcov_high
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e10 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_e10_allcov_file_high, replace): second_e10_allcov_high
	estimates save second_e10_allcov_high , replace
	esttab using "Main_e10_second_boot_high.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Main specification e10 second Intervention) nomtitles addnotes(Notes follow)

* #################### ROBUSTNESS - BANDWIDTH & DONUT #########################*
est clear
	* 240 hour (10 day) bandwidth
	scalar h = 240 // Bandwidth 240 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program second_10days_high
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e10 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_10days_high, replace): second_10days_high
	estimates save second_10days_high, replace
	* 288 hour (12 day) bandwidth
	scalar h = 288 // Bandwidth 288 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program second_12days_high
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e10 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_12days_high, replace): second_12days_high
	estimates save second_12days_high, replace
	* 336 hour (14 day) bandwidth (no donut)
	scalar da = 0 // Donut 0 hours
	scalar db = -da
	scalar h = 336 // Bandwidth 336 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program second_14days_0d_high
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e10 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_14days_0d_high, replace): second_14days_0d_high
	estimates save second_14days_0d_high, replace
	* 336 hour bandwidth (9 hours donut) siehe Average Effects
	estimates use "second_e10_allcov_high"
	eststo: estimates replay "second_e10_allcov_high"
	* 336 hour (14 day) bandwidth (18 hours donut)
	scalar da = 18 // Donut 18 hours
	scalar db = -da
	scalar h = 336 // Bandwidth 336 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program second_14days_18d_high
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e10 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_14days_18d_high, replace): second_14days_18d_high
	estimates save second_14days_18d_high, replace
	* 384 hour (16 day) bandwidth
	scalar da = 9 // Donut 9 hours
	scalar db = -da
	scalar h = 384 // Bandwidth 384 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program second_16days_high
			reghdfe mean_e10 $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_e10)
			reg res_mean_e10 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_e10 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_16days_high, replace): second_16days_high
	estimates save second_16days_high, replace
		esttab using "ROBUST_e10_second_boot_high.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Heterogeneous Effects second Intervention) nomtitles addnotes(Notes follow)

*-------------------------------------------------------------------------------
*############################# DIESEL #########################################*
*-------------------------------------------------------------------------------
*######################### AVERAGE EFFECTS ####################################*
scalar da = 9 // Donut 9 hours (baseline)
scalar db = -da
scalar h = 336 // Bandwidth 336 hours (baseline)
drop kernel_wgt
gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h)) // triangular kernel (non-parametric estimation)
	replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
	replace kernel_wgt = 0 if mi(kernel_wgt)
		* Seasonality
est clear
	capture program drop _all
		program second_diesel_season_high
			reghdfe mean_diesel $var_seas , absorb(station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_diesel
		end 	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_diesel_season_file_high, replace): second_diesel_season_high
	estimates save second_diesel_season_high , replace
	* Supply & Seasonality
		program second_diesel_sup_high
			reghdfe mean_diesel $var_sup $var_seas , absorb(station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_diesel
		end	
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new) saving(second_diesel_sup_file_high, replace): second_diesel_sup_high
	estimates save second_diesel_sup_high , replace
	* Demand & Seasonality
		program second_diesel_dem_high
			reghdfe mean_diesel $var_dem $var_seas , absorb(station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_diesel
		end
		xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new) saving(second_diesel_sem_file_high, replace): second_diesel_dem_high
	estimates save second_diesel_dem_high , replace
	* Supply & Demand & Seasonality
	capture program drop _all
		program second_diesel_allcov_high
			reghdfe mean_diesel $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_diesel 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_diesel_allcov_file_high, replace): second_diesel_allcov_high
	estimates save second_diesel_allcov_high , replace
	esttab using "Main_diesel_second_boot_high.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Main specification diesel second Intervention) nomtitles addnotes(Notes follow)

* #################### ROBUSTNESS - BANDWIDTH & DONUT #########################*
est clear
	* 240 hour (10 day) bandwidth
	scalar h = 240 // Bandwidth 240 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program second_10days_high
			reghdfe mean_diesel $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_diesel 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_10days_high, replace): second_10days_high
	estimates save second_10days_high, replace
	* 288 hour (12 day) bandwidth
	scalar h = 288 // Bandwidth 288 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program second_12days_high
			reghdfe mean_diesel $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_diesel 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_12days_high, replace): second_12days_high
	estimates save second_12days_high, replace
	* 336 hour (14 day) bandwidth (no donut)
	scalar da = 0 // Donut 0 hours
	scalar db = -da
	scalar h = 336 // Bandwidth 336 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program second_14days_0d_high
			reghdfe mean_diesel $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_diesel 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_14days_0d_high, replace): second_14days_0d_high
	estimates save second_14days_0d_high, replace
	* 336 hour bandwidth (9 hours donut) siehe Average Effects
	estimates use "second_diesel_allcov_high"
	eststo: estimates replay "second_diesel_allcov_high"
	* 336 hour (14 day) bandwidth (18 hours donut)
	scalar da = 18 // Donut 18 hours
	scalar db = -da
	scalar h = 336 // Bandwidth 336 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program second_14days_18d_high
			reghdfe mean_diesel $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_diesel 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_14days_18d_high, replace): second_14days_18d_high
	estimates save second_14days_18d_high, replace
	* 384 hour (16 day) bandwidth
	scalar da = 9 // Donut 9 hours
	scalar db = -da
	scalar h = 384 // Bandwidth 384 hours
	drop kernel_wgt
	gen kernel_wgt = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	capture program drop _all
		program second_16days_high
			reghdfe mean_diesel $var_sup $var_dem $var_seas , absorb( station_id) resid(res_mean_diesel)
			reg res_mean_diesel 1.tankrabatt_2##(c.hour_to_treat_2##c.hour_to_treat_2) [aw=kernel_wgt] ///
				if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
			drop res_mean_diesel 
		end			
	xtset, clear
	eststo: bootstrap, $boot_rep cluster(station_id) idcluster(id_new)  saving(second_16days_high, replace): second_16days_high
	estimates save second_16days_high, replace
		esttab using "ROBUST_diesel_second_boot_high.rtf", replace b(2) se(2) ar2 label star(* 0.10 ** 0.05 *** 0.01) title(Heterogeneous Effects second Intervention) nomtitles addnotes(Notes follow)
		
		
timer off 1
timer list
log close