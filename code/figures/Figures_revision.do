*** ######################################################################## ***
** FIGURES publication ready **
		* Install ftools (remove program if it existed previously)
		cap ado uninstall ftools
		net install ftools, from("https://raw.githubusercontent.com/sergiocorreia/ftools/master/src/")

		* Install reghdfe 6.x
		cap ado uninstall reghdfe
		net install reghdfe, from("https://raw.githubusercontent.com/sergiocorreia/reghdfe/master/src/")

		* Install parallel, if using the parallel() option; don't install from SSC
		cap ado uninstall parallel
		net install parallel, from(https://raw.github.com/gvegayon/parallel/stable/) replace
		mata mata mlib index

*** ######################################################################## ***
* MAP GERMANY GASOLINE STATIONS - HHI (4km)
grmap, activate
cd "C:\Users\frede\Desktop\Bachelor-Thesis\Daten\vg2500_geo84 (2)"
use vg2500_bld.dta, clear
	grmap USE, fcolor(gs16) ocolor(gs1) osize(vthin) legenda(off) legend(pos(6) col(3) colgap(5) size(small)) plotregion(margin(0 0 0 0)) point(data("C:\Users\frede\Desktop\Bachelor-Thesis\Daten\tankstellenpreise\stations_with_pricechanges_2022_HHI.dta") xcoord(longitude) ycoord(latitude) by(localHHId_4) size(1.1pt 1.1pt 1.1pt 1.1pt 1.1pt 1.1pt 1.1pt 1.1pt 1.1pt 1.1pt) osize(none none none none none none none none none none) fcolor(BuYlRd) legenda(off) legcount) saving("C:\Users\frede\Desktop\Bachelor-Thesis\Two_Stage_RDiT_JTEP\revised figures_hquality\Map1.gph", replace)
	graph export "C:\Users\frede\Desktop\Bachelor-Thesis\Two_Stage_RDiT_JTEP\revised figures_hquality\Map1.pdf", replace

* MAP GERMANY GASOLINE STATIONS - ROADSIDE/HIGHWAY
grmap USE, fcolor(gs16) ocolor(gs1) osize(vthin) legenda(off) legend(pos(6) col(2) colgap(5) size(small)) plotregion(margin(0 0 8 0)) point(data("C:\Users\frede\Desktop\Bachelor-Thesis\Daten\tankstellenpreise\stations_with_pricechanges_2022.dta") xcoord(longitude) ycoord(latitude) by(autobahn) size(.8pt 1.27pt) osize(none none) fcolor("0 167 153*0.7" "233 77 61") shape(O D) legenda(on) legcount) saving("C:\Users\frede\Desktop\Bachelor-Thesis\Two_Stage_RDiT_JTEP\revised figures_hquality\Map2.gph", replace)
	graph export "C:\Users\frede\Desktop\Bachelor-Thesis\Two_Stage_RDiT_JTEP\revised figures_hquality\Map2.pdf", replace

*** ######################################################################## ***
* PRICE CYCLES:
cd "C:\Users\frede\Desktop\Bachelor-Thesis\Daten\tankstellenpreise"
use "roadside_price_2022_0309_openandassume_revised.dta", clear
append using "roadside_price_2022_0309_openandassume_high_revised.dta"
replace autobahn = 0 if mi(autobahn)
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
keep if inrange(month,5,9)

xtile hhi_km_quantile = hhi_4km, nquantiles(4)
foreach i in mean_diesel mean_e5 {
bysort autobahn hhi_km_quantile date station_id: egen day_`i' = mean(`i') 
gen demean_`i' = `i' - day_`i'
bysort autobahn hhi_km_quantile hour: egen mean_dem_`i' = mean(demean_`i')
egen q1_`i' = pctile(demean_`i'), p(25) by(hour autobahn hhi_km_quantile)
egen q3_`i' = pctile(demean_`i'), p(75) by(hour autobahn hhi_km_quantile)
}

drop day_mean_diesel demean_mean_diesel mean_dem_mean_diesel q1_mean_diesel q3_mean_diesel day_mean_e5 demean_mean_e5 mean_dem_mean_e5 q1_mean_e5 q3_mean_e5
foreach i in mean_diesel mean_e5 {
bysort autobahn date station_id: egen day_`i' = mean(`i') 
gen demean_`i' = `i' - day_`i'
bysort autobahn hour: egen mean_dem_`i' = mean(demean_`i')
egen q1_`i' = pctile(demean_`i'), p(25) by(hour autobahn)
egen q3_`i' = pctile(demean_`i'), p(75) by(hour autobahn)
}
sort station_id time

twoway line q1_mean_e5 q3_mean_e5 hour if station_id == 2 & date == 22772 & autobahn == 0, lcolor(stc1%60 stc1%60) lwidth(medthick) lpattern(dash dash) lwidth(medthick medthick) || ///
	   line q1_mean_diesel q3_mean_diesel hour if station_id == 2 & date == 22772 & autobahn == 0, lcolor(stc2%60 stc2%60) lwidth(medthick) lpattern(dash dash) lwidth(medthick medthick) || /// 
	   line mean_dem_mean_e5 hour if station_id == 2 & date == 22772 & autobahn == 0, lcolor(stc1) lwidth(medthick) || /// 
	   line mean_dem_mean_diesel hour if station_id == 2 & date == 22772 & autobahn == 0, lcolor(stc2) lwidth(medthick) ylabel(-4(4)11) ytitle(Demeaned Price (ct/l)) xtitle(Hour of the Day) xlabel(0 "00:00" 3 "03:00" 6 "06.00" 9 "09:00" 12 "12:00" 15 "15:00"  18 "18:00" 21 "21:00", angle(-45)) legend(row(1) pos(6) label(5 "Super e5") label(6 "Diesel") order(5 6)) name(hourofday_road, replace) xsize(8.0cm) ysize(6cm) scale(1.3) saving("C:\Users\frede\Desktop\Bachelor-Thesis\Two_Stage_RDiT_JTEP\revised figures_hquality\price_cycle_road_e5_diesel.gph", replace)
graph export "C:\Users\frede\Desktop\Bachelor-Thesis\Two_Stage_RDiT_JTEP\revised figures_hquality\price_cycle_road_e5_diesel.emf", replace
	   
twoway line q1_mean_e5 q3_mean_e5 hour if station_id == 1413 & date == 22772 & autobahn == 1, lcolor(stc1%60 stc1%60) lwidth(medthick) lpattern(dash dash) lwidth(medthick medthick) || ///
	   line q1_mean_diesel q3_mean_diesel hour if station_id == 1413 & date == 22772 & autobahn == 1, lcolor(stc2%60 stc2%60) lwidth(medthick) lpattern(dash dash) lwidth(medthick medthick) || /// 
	   line mean_dem_mean_e5 hour if station_id == 1413 & date == 22772 & autobahn == 1, lcolor(stc1) lwidth(medthick) || /// 
	   line mean_dem_mean_diesel hour if station_id == 1413 & date == 22772 & autobahn == 1, lcolor(stc2) lwidth(medthick) ytitle(Demeaned Price (ct/l)) xtitle(Hour of the Day) xlabel(0 "00:00" 3 "03:00" 6 "06.00" 9 "09:00" 12 "12:00" 15 "15:00"  18 "18:00" 21 "21:00", angle(-45)) legend(row(1) pos(6) label(5 "Super e5") label(6 "Diesel") order(5 6)) name(hourofday_high, replace) xsize(8.0cm) ysize(6cm) scale(1.3) saving("C:\Users\frede\Desktop\Bachelor-Thesis\Two_Stage_RDiT_JTEP\revised figures_hquality\price_cycle_high_e5_diesel.gph", replace)
graph export "C:\Users\frede\Desktop\Bachelor-Thesis\Two_Stage_RDiT_JTEP\revised figures_hquality\price_cycle_high_e5_diesel.emf", replace

*** ######################################################################## ***
* RDiT Plots:
clear all
cd "C:\Users\frede\Desktop\Bachelor-Thesis\Daten\tankstellenpreise"
use "roadside_price_2022_0309_openandassume_revised.dta", clear
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
bysort station_id (time): gen L24_BRENT_FOB = L24.BRENT_FOB
bysort station_id (time): gen L24_usd_to_euro = L24.usd_to_euro
bysort station_id (time): gen L1_kfz_menge_B = L1.kfz_menge_B
gen hour_to_treat_1 = time
replace hour_to_treat_1 = (time- 11641) // Treatment June 1st 2022
gen tankrabatt_1 = 0
replace tankrabatt_1 = 1 if inrange(month,6,8)
gen hour_to_treat_2 = time
replace hour_to_treat_2 = (time- 13849) // Treatment September 1st 2022
gen tankrabatt_2 = 1
replace tankrabatt_2 = 0 if inrange(month,6,8)
scalar h = 336 // Bandwidth 336 hours (baseline)
gen kernel_wgt_1 = 1- (hour_to_treat_1)/h if inrange(hour_to_treat_1,0,scalar(h)) // triangular kernel (non-parametric estimation) FIRST
	replace kernel_wgt_1 = 1- (hour_to_treat_1)/-h if inrange(hour_to_treat_1,scalar(-h),0)
	replace kernel_wgt_1 = 0 if mi(kernel_wgt_1)
gen kernel_wgt_2 = 1- (hour_to_treat_2)/h if inrange(hour_to_treat_2,0,scalar(h)) // triangular kernel (non-parametric estimation) SECOND
	replace kernel_wgt_2 = 1- (hour_to_treat_2)/-h if inrange(hour_to_treat_2,scalar(-h),0)
	replace kernel_wgt_2 = 0 if mi(kernel_wgt_2)
	
global var_seas i.hour i.week_days // seasonality and price cycles covariates
global var_sup BRENT_FOB L24_BRENT_FOB usd_to_euro L24_usd_to_euro // supply covariates
global var_dem kfz_menge_B L1_kfz_menge_B temperature precipitation holiday_state day_before_holiday_state // demand covariates

*FIRST INTERVENTION
	reghdfe mean_e5 $var_seas if inrange(time,9961,13321), absorb(station_id) resid(res_mean_e5_first)
	rdplot res_mean_e5_first hour_to_treat_1 if (inrange(hour_to_treat_1, -336,-9) | inrange(hour_to_treat_1,9,336)), c(0) h(336) p(2) kernel(triangular) nbins(27 27) masspoints(adjust) graph_options(ytitle(Residuals E5 (ct/l)) xtitle(Hours to Treatment) xlabel(-336(168)336) legend(off) scale(1.4) xsize(8.0cm) ysize(5cm) name(mean_e5_first, replace))
	graph export "C:\Users\frede\Desktop\Bachelor-Thesis\Two_Stage_RDiT_JTEP\revised figures_hquality\rdplot_e5_first_2.emf", replace
	reghdfe mean_e10 $var_seas if inrange(time,9961,13321), absorb(station_id) resid(res_mean_e10_first)
	rdplot res_mean_e10_first hour_to_treat_1 if (inrange(hour_to_treat_1, -336,-9) | inrange(hour_to_treat_1,9,336)), c(0) h(336) p(2) kernel(triangular) nbins(27 27) masspoints(adjust) graph_options(ytitle(Residuals E10 (ct/l)) xtitle(Hours to Treatment) xlabel(-336(168)336) legend(off) scale(1.4) xsize(8.0cm) ysize(5cm) name(mean_e10_first, replace))
	graph export "C:\Users\frede\Desktop\Bachelor-Thesis\Two_Stage_RDiT_JTEP\revised figures_hquality\rdplot_e10_first_2.emf", replace
	reghdfe mean_diesel $var_seas if inrange(time,9961,13321), absorb(station_id) resid(res_mean_diesel_first)
	rdplot res_mean_diesel_first hour_to_treat_1 if (inrange(hour_to_treat_1, -336,-9) | inrange(hour_to_treat_1,9,336)), c(0) h(336) p(2) kernel(triangular) nbins(27 27) masspoints(adjust) graph_options(ytitle(Residuals Diesel (ct/l)) xtitle(Hours to Treatment) xlabel(-336(168)336) legend(off) scale(1.4) xsize(8.0cm) ysize(5cm) name(mean_diesel_first, replace))
	graph export "C:\Users\frede\Desktop\Bachelor-Thesis\Two_Stage_RDiT_JTEP\revised figures_hquality\rdplot_diesel_first_2.emf", replace	
*SECOND INTERVENTION
	reghdfe mean_e5 $var_seas if inrange(time,12169,15529), absorb(station_id) resid(res_mean_e5_second)
	rdplot res_mean_e5_second hour_to_treat_2 if (inrange(hour_to_treat_2, -336,-9) | inrange(hour_to_treat_2,9,336)), c(0) h(336) p(2) kernel(triangular) nbins(27 27) masspoints(adjust) graph_options(ytitle(Residuals E5 (ct/l)) xtitle(Hours to Treatment) xlabel(-336(168)336) legend(off) scale(1.4) xsize(8.0cm) ysize(5cm) name(mean_e5_second, replace))
	graph export "C:\Users\frede\Desktop\Bachelor-Thesis\Two_Stage_RDiT_JTEP\revised figures_hquality\rdplot_e5_second_2.emf", replace
	reghdfe mean_e10 $var_seas if inrange(time,12169,15529), absorb(station_id) resid(res_mean_e10_second)
	rdplot res_mean_e10_second hour_to_treat_2 if (inrange(hour_to_treat_2, -336,-9) | inrange(hour_to_treat_2,9,336)), c(0) h(336) p(2) kernel(triangular) nbins(27 27) masspoints(adjust) graph_options(ytitle(Residuals E10 (ct/l)) xtitle(Hours to Treatment) xlabel(-336(168)336) legend(off) scale(1.4) xsize(8.0cm) ysize(5cm) name(mean_e10_second, replace))
	graph export "C:\Users\frede\Desktop\Bachelor-Thesis\Two_Stage_RDiT_JTEP\revised figures_hquality\rdplot_e10_second_2.emf", replace
	reghdfe mean_diesel $var_seas if inrange(time,12169,15529), absorb(station_id) resid(res_mean_diesel_second)
	rdplot res_mean_diesel_second hour_to_treat_2 if (inrange(hour_to_treat_2, -336,-9) | inrange(hour_to_treat_2,9,336)), c(0) h(336) p(2) kernel(triangular) nbins(27 27) masspoints(adjust) graph_options(ytitle(Residuals Diesel (ct/l)) xtitle(Hours to Treatment) xlabel(-336(168)336) legend(off) scale(1.4) xsize(8.0cm) ysize(5cm) name(mean_diesel_second, replace))
	graph export "C:\Users\frede\Desktop\Bachelor-Thesis\Two_Stage_RDiT_JTEP\revised figures_hquality\rdplot_diesel_second_2.emf", replace	

*** ######################################################################## ***
* ROBUST - POLYNOMIAL RDPLOT
scalar h = 336 // Bandwidth 336 hours (baseline)
scalar da = 9 // Donut 9 hours
scalar db = -da
foreach i in mean_diesel mean_e5 mean_e10 {
*FIRST INTERVENTION
	* Extracting bins from RDPLOT:
	rdplot res_`i'_first hour_to_treat_1 if (inrange(hour_to_treat_1, -336,-9) | inrange(hour_to_treat_1,9,336)), c(0) h(336) p(1) kernel(triangular) nbins(27 27) masspoints(adjust) hide genvars
	rename rdplot_mean_y rd_`i'_y_first
	rename rdplot_mean_x rd_`i'_x_first
	drop rdplot_id rdplot_N rdplot_min_bin rdplot_max_bin rdplot_mean_bin rdplot_se_y rdplot_ci_l rdplot_ci_r rdplot_hat_y
	* Linear Fit
	reg res_`i'_first 1.tankrabatt_1##c.hour_to_treat_1 [aw=kernel_wgt_1] ///
										if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
	predict res_`i'_lin_first, xb
	* Quadratic Fit
	reg res_`i'_first 1.tankrabatt_1##c.hour_to_treat_1##c.hour_to_treat_1 [aw=kernel_wgt_1] ///
										if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
	predict res_`i'_quad_first, xb
	* Cubic Fit
	reg res_`i'_first 1.tankrabatt_1##c.hour_to_treat_1##c.hour_to_treat_1##c.hour_to_treat_1 [aw=kernel_wgt_1] ///
										if (inrange(hour_to_treat_1,scalar(-h),scalar(db)) | inrange(hour_to_treat_1,scalar(da),scalar(h)))
	predict res_`i'_cub_first, xb
*SECOND INTERVENTION
	* Extracting bins from RDPLOT:
	rdplot res_`i'_second hour_to_treat_2 if (inrange(hour_to_treat_2, -336,-9) | inrange(hour_to_treat_2,9,336)), c(0) h(336) p(1) kernel(triangular) nbins(27 27) masspoints(adjust) hide genvars
	rename rdplot_mean_y rd_`i'_y_second
	rename rdplot_mean_x rd_`i'_x_second
	drop rdplot_id rdplot_N rdplot_min_bin rdplot_max_bin rdplot_mean_bin rdplot_se_y rdplot_ci_l rdplot_ci_r rdplot_hat_y
	* Linear Fit
	reg res_`i'_second 1.tankrabatt_2##c.hour_to_treat_2 [aw=kernel_wgt_2] ///
										if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
	predict res_`i'_lin_second, xb
	* Quadratic Fit
	reg res_`i'_second 1.tankrabatt_2##c.hour_to_treat_2##c.hour_to_treat_2 [aw=kernel_wgt_2] ///
										if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
	predict res_`i'_quad_second, xb
	* Cubic Fit
	reg res_`i'_second 1.tankrabatt_2##c.hour_to_treat_2##c.hour_to_treat_2##c.hour_to_treat_2 [aw=kernel_wgt_2] ///
										if (inrange(hour_to_treat_2,scalar(-h),scalar(db)) | inrange(hour_to_treat_2,scalar(da),scalar(h)))
	predict res_`i'_cub_second, xb
}
* POLYNOMIAL RDPLOT:
*Diesel
use "C:\Users\frede\Desktop\Bachelor-Thesis\Daten\tankstellenpreise\residual first stage_polynm_rob.dta", clear

xtset station_id time
twoway (scatter rd_mean_diesel_y_first rd_mean_diesel_x_first if station_id == 2, color(gs10) xline(0, lpattern(dash))) || /// 
	   (line res_mean_diesel_lin_first hour_to_treat_1 if inrange(hour_to_treat_1,scalar(-h),scalar(-1)) & station_id == 2, color(stc1)) || /// 
	   (line res_mean_diesel_lin_first hour_to_treat_1 if inrange(hour_to_treat_1,scalar(0),scalar(h)) & station_id == 2, color(stc1)) || ///
	   (line res_mean_diesel_cub_first hour_to_treat_1 if inrange(hour_to_treat_1,scalar(-h),scalar(-1)) & station_id == 2, color(stc2)) || /// 
	   (line res_mean_diesel_cub_first hour_to_treat_1 if inrange(hour_to_treat_1,scalar(0),scalar(h)) & station_id == 2, color(stc2)) || ///
	   (line res_mean_diesel_quad_first hour_to_treat_1 if inrange(hour_to_treat_1,scalar(-h),scalar(-1)) & station_id == 2, color(stc3)) || /// 
	   (line res_mean_diesel_quad_first hour_to_treat_1 if inrange(hour_to_treat_1,scalar(0),scalar(h)) & station_id == 2, color(stc3) legend(pos(6) col(4) off) xtitle("Hours to Treatment") ytitle("Residuals Diesel (ct/l)") ysize(5cm) xsize(8cm) scale(1.3) xlabel(-336(168)336) name(rd_poly_diesel_first, replace))
	   	graph export "C:\Users\frede\Desktop\Bachelor-Thesis\Two_Stage_RDiT_JTEP\revised figures_hquality\poly_rdplot_diesel_first.emf", replace	
twoway (scatter rd_mean_diesel_y_second rd_mean_diesel_x_second if station_id == 2, color(gs10) xline(0, lpattern(dash))) || /// 
	   (line res_mean_diesel_lin_second hour_to_treat_2 if inrange(hour_to_treat_2,scalar(-h),scalar(-1)) & station_id == 2, color(stc1)) || /// 
	   (line res_mean_diesel_lin_second hour_to_treat_2 if inrange(hour_to_treat_2,scalar(0),scalar(h)) & station_id == 2, color(stc1)) || ///
	   (line res_mean_diesel_cub_second hour_to_treat_2 if inrange(hour_to_treat_2,scalar(-h),scalar(-1)) & station_id == 2, color(stc2)) || /// 
	   (line res_mean_diesel_cub_second hour_to_treat_2 if inrange(hour_to_treat_2,scalar(0),scalar(h)) & station_id == 2, color(stc2)) || ///
	   (line res_mean_diesel_quad_second hour_to_treat_2 if inrange(hour_to_treat_2,scalar(-h),scalar(-1)) & station_id == 2, color(stc3)) || /// 
	   (line res_mean_diesel_quad_second hour_to_treat_2 if inrange(hour_to_treat_2,scalar(0),scalar(h)) & station_id == 2, color(stc3) legend(pos(6) col(4) off) xtitle("Hours to Treatment") ytitle("Residuals Diesel (ct/l)") ysize(5cm) xsize(8cm) scale(1.3) xlabel(-336(168)336) name(rd_poly_diesel_second, replace))
	   	graph export "C:\Users\frede\Desktop\Bachelor-Thesis\Two_Stage_RDiT_JTEP\revised figures_hquality\poly_rdplot_diesel_second.emf", replace
*E5   
twoway (scatter rd_mean_e5_y_first rd_mean_e5_x_first if station_id == 2, color(gs10) xline(0, lpattern(dash))) || /// 
	   (line res_mean_e5_lin_first hour_to_treat_1 if inrange(hour_to_treat_1,scalar(-h),scalar(-1)) & station_id == 2, color(stc1)) || /// 
	   (line res_mean_e5_lin_first hour_to_treat_1 if inrange(hour_to_treat_1,scalar(0),scalar(h)) & station_id == 2, color(stc1)) || ///
	   (line res_mean_e5_cub_first hour_to_treat_1 if inrange(hour_to_treat_1,scalar(-h),scalar(-1)) & station_id == 2, color(stc2)) || /// 
	   (line res_mean_e5_cub_first hour_to_treat_1 if inrange(hour_to_treat_1,scalar(0),scalar(h)) & station_id == 2, color(stc2)) || ///
	   (line res_mean_e5_quad_first hour_to_treat_1 if inrange(hour_to_treat_1,scalar(-h),scalar(-1)) & station_id == 2, color(stc3)) || /// 
	   (line res_mean_e5_quad_first hour_to_treat_1 if inrange(hour_to_treat_1,scalar(0),scalar(h)) & station_id == 2, color(stc3) legend(pos(6) col(4) off) xtitle("Hours to Treatment") ytitle("Residuals E5(ct/l)") ysize(5cm) xsize(8cm) scale(1.3) xlabel(-336(168)336) name(rd_poly_e5_first, replace))
	   graph export "C:\Users\frede\Desktop\Bachelor-Thesis\Two_Stage_RDiT_JTEP\revised figures_hquality\poly_rdplot_e5_first.emf", replace
twoway (scatter rd_mean_e5_y_second rd_mean_e5_x_second if station_id == 2, color(gs10) xline(0, lpattern(dash))) || /// 
	   (line res_mean_e5_lin_second hour_to_treat_2 if inrange(hour_to_treat_2,scalar(-h),scalar(-1)) & station_id == 2, color(stc1)) || /// 
	   (line res_mean_e5_lin_second hour_to_treat_2 if inrange(hour_to_treat_2,scalar(0),scalar(h)) & station_id == 2, color(stc1)) || ///
	   (line res_mean_e5_cub_second hour_to_treat_2 if inrange(hour_to_treat_2,scalar(-h),scalar(-1)) & station_id == 2, color(stc2)) || /// 
	   (line res_mean_e5_cub_second hour_to_treat_2 if inrange(hour_to_treat_2,scalar(0),scalar(h)) & station_id == 2, color(stc2)) || ///
	   (line res_mean_e5_quad_second hour_to_treat_2 if inrange(hour_to_treat_2,scalar(-h),scalar(-1)) & station_id == 2, color(stc3)) || /// 
	   (line res_mean_e5_quad_second hour_to_treat_2 if inrange(hour_to_treat_2,scalar(0),scalar(h)) & station_id == 2, color(stc3) legend(pos(6) col(4) off) xtitle("Hours to Treatment") ytitle("Residuals E5 (ct/l)") ysize(5cm) xsize(8cm) scale(1.3) xlabel(-336(168)336) name(rd_poly_e5_second, replace))
	   graph export "C:\Users\frede\Desktop\Bachelor-Thesis\Two_Stage_RDiT_JTEP\revised figures_hquality\poly_rdplot_e5_second.emf", replace
*E10
twoway (scatter rd_mean_e10_y_first rd_mean_e10_x_first if station_id == 2, color(gs10) xline(0, lpattern(dash)) ) || /// 
	   (line res_mean_e10_lin_first hour_to_treat_1 if inrange(hour_to_treat_1,scalar(-h),scalar(-1)) & station_id == 2, color(stc1)) || /// 
	   (line res_mean_e10_lin_first hour_to_treat_1 if inrange(hour_to_treat_1,scalar(0),scalar(h)) & station_id == 2, color(stc1)) || ///
	   (line res_mean_e10_cub_first hour_to_treat_1 if inrange(hour_to_treat_1,scalar(-h),scalar(-1)) & station_id == 2, color(stc2)) || /// 
	   (line res_mean_e10_cub_first hour_to_treat_1 if inrange(hour_to_treat_1,scalar(0),scalar(h)) & station_id == 2, color(stc2)) || ///
	   (line res_mean_e10_quad_first hour_to_treat_1 if inrange(hour_to_treat_1,scalar(-h),scalar(-1)) & station_id == 2, color(stc3)) || /// 
	   (line res_mean_e10_quad_first hour_to_treat_1 if inrange(hour_to_treat_1,scalar(0),scalar(h)) & station_id == 2, color(stc3) legend(pos(6) col(4) off) xtitle("Hours to Treatment") ytitle("Residuals E10(ct/l)") ysize(5cm) xsize(8cm) scale(1.3) xlabel(-336(168)336) name(rd_poly_e10_first, replace))
	   graph export "C:\Users\frede\Desktop\Bachelor-Thesis\Two_Stage_RDiT_JTEP\revised figures_hquality\poly_rdplot_e10_first.emf", replace
twoway (scatter rd_mean_e10_y_second rd_mean_e10_x_second if station_id == 2, color(gs10) xline(0, lpattern(dash))) || /// 
	   (line res_mean_e10_lin_second hour_to_treat_2 if inrange(hour_to_treat_2,scalar(-h),scalar(-1)) & station_id == 2, color(stc1)) || /// 
	   (line res_mean_e10_lin_second hour_to_treat_2 if inrange(hour_to_treat_2,scalar(0),scalar(h)) & station_id == 2, color(stc1)) || ///
	   (line res_mean_e10_cub_second hour_to_treat_2 if inrange(hour_to_treat_2,scalar(-h),scalar(-1)) & station_id == 2, color(stc2)) || /// 
	   (line res_mean_e10_cub_second hour_to_treat_2 if inrange(hour_to_treat_2,scalar(0),scalar(h)) & station_id == 2, color(stc2)) || ///
	   (line res_mean_e10_quad_second hour_to_treat_2 if inrange(hour_to_treat_2,scalar(-h),scalar(-1)) & station_id == 2, color(stc3)) || /// 
	   (line res_mean_e10_quad_second hour_to_treat_2 if inrange(hour_to_treat_2,scalar(0),scalar(h)) & station_id == 2, color(stc3) legend(pos(6) col(4) off) xtitle("Hours to Treatment") ytitle("Residuals E10 (ct/l)") ysize(5cm) xsize(8cm) scale(1.3) xlabel(-336(168)336) name(rd_poly_e10_second, replace))
	   graph export "C:\Users\frede\Desktop\Bachelor-Thesis\Two_Stage_RDiT_JTEP\revised figures_hquality\poly_rdplot_e10_second.emf", replace
	
save "residual first stage_polynm_rob.dta"

*** ######################################################################## ***
* RDiT Plot - Alternative Polynomials
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
keep if inrange(time,9625,15865) // 11 Weeks (12 Weeks pre treatment 9625 (10 weeks after tax cut 13345), 1 Week with treatment 11809) (12 weeks after 2nd treatment 15865 
gen lnBRENT_FOB = ln(BRENT_FOB)
bysort station_id (time): gen L24_BRENT_FOB = L24.BRENT_FOB
bysort station_id (time): gen L24_lnBRENT_FOB = L24.lnBRENT_FOB
bysort station_id (time): gen L24_usd_to_euro = L24.usd_to_euro
bysort station_id (time): gen L1_kfz_menge_B = L1.kfz_menge_B
gen hour_to_treat_1 = time
replace hour_to_treat_1 = (time- 11641) // Treatment June 1st 2022
gen tankrabatt_1 = 0
replace tankrabatt_1 = 1 if inrange(month,6,8)
gen hour_to_treat_2 = time
replace hour_to_treat_2 = (time- 13849) // Treatment September 1st 2022
gen tankrabatt_2 = 1
replace tankrabatt_2 = 0 if inrange(month,6,8)
global var_seas i.hour i.week_days
global var_sup lnBRENT_FOB L24_lnBRENT_FOB usd_to_euro L24_usd_to_euro
global var_dem kfz_menge_B L1_kfz_menge_B temperature precipitation holiday_state day_before_holiday_state
* tax cut
*pre treatment
forval b = 1/12 {
	scalar minus`b'week_1= 11641 - (24*7*`b')
	di minus`b'week_1
}
*after treatment:
forval a = 1/12 {
	scalar plus`a'week_1= 11641 + (24*7*`a')
	di plus`a'week_1
}
* tax increase
*pre treatment
forval b = 1/12 {
	scalar minus`b'week_2= 13849 - (24*7*`b')
	di minus`b'week_2
}
*after treatment:
forval a = 1/12 {
	scalar plus`a'week_2= 13849 + (24*7*`a')
	di plus`a'week_2
}

scalar da = 9 // Donut 9 hours
scalar h = 336 // Bandwidth 120 hours
scalar db = -da
foreach inter in 1 2 {
	gen kernel_wgt = 1- (hour_to_treat_`inter')/h if inrange(hour_to_treat_`inter',0,scalar(h))
		replace kernel_wgt = 1- (hour_to_treat_`inter')/-h if inrange(hour_to_treat_`inter',scalar(-h),0)
		replace kernel_wgt = 0 if mi(kernel_wgt)
	foreach type in e5 e10 diesel { 
		reghdfe mean_`type' $var_sup $var_dem $var_seas if inrange(time,minus10week_`inter', plus10week_`inter'), absorb(station_id) resid(res_mean_`type')
		* Linear Fit
		reg res_mean_`type' 1.tankrabatt_`inter'##c.hour_to_treat_`inter' [aw=kernel_wgt] ///
											if (inrange(hour_to_treat_`inter',scalar(-h),scalar(db)) | inrange(hour_to_treat_`inter',scalar(da),scalar(h)))
		predict res_mean_`type'_lin_`inter', xb
		* Quadratic Fit
		reg res_mean_`type' 1.tankrabatt_`inter'##c.hour_to_treat_`inter'##c.hour_to_treat_`inter' [aw=kernel_wgt] ///
											if (inrange(hour_to_treat_`inter',scalar(-h),scalar(db)) | inrange(hour_to_treat_`inter',scalar(da),scalar(h)))
		predict res_mean_`type'_quad_`inter', xb
		* Cubic Fit
		reg res_mean_`type' 1.tankrabatt_`inter'##c.hour_to_treat_`inter'##c.hour_to_treat_`inter'##c.hour_to_treat_`inter' [aw=kernel_wgt] ///
											if (inrange(hour_to_treat_`inter',scalar(-h),scalar(db)) | inrange(hour_to_treat_`inter',scalar(da),scalar(h)))
		predict res_mean_`type'_cub_`inter', xb
		gen bin = ceil(hour_to_treat_`inter'/4)*4
		bysort bin: egen res_mean_`type'_bin = mean(res_mean_`type')
		bysort station_id: replace res_mean_`type'_bin = . if res_mean_`type'_bin[_n-1] == res_mean_`type'_bin
		bysort station_id: replace res_mean_`type'_bin = . if res_mean_`type'_bin[_n-2] == res_mean_`type'_bin
		bysort station_id: replace res_mean_`type'_bin = . if res_mean_`type'_bin[_n-3] == res_mean_`type'_bin
		xtset station_id time
		twoway (scatter res_mean_`type'_bin hour_to_treat_`inter' if (inrange(hour_to_treat_`inter',-h,0) | inrange(hour_to_treat_`inter',0,h)) & station_id == 2, color(gs10) xline(0, lpattern(dash)) xline(-9, lpattern(dash) lcolor(gs15)) xline(9, lpattern(dash) lcolor(gs15))) || /// 
			   (line res_mean_`type'_lin_`inter' hour_to_treat_`inter' if inrange(hour_to_treat_`inter',scalar(-h),scalar(-1)) & station_id == 2, color(stc1)) || /// 
			   (line res_mean_`type'_lin_`inter' hour_to_treat_`inter' if inrange(hour_to_treat_`inter',scalar(0),scalar(h)) & station_id == 2, color(stc1)) || ///
			   (line res_mean_`type'_cub_`inter' hour_to_treat_`inter' if inrange(hour_to_treat_`inter',scalar(-h),scalar(-1)) & station_id == 2, color(stc2)) || /// 
			   (line res_mean_`type'_cub_`inter' hour_to_treat_`inter' if inrange(hour_to_treat_`inter',scalar(0),scalar(h)) & station_id == 2, color(stc2)) || ///
			   (line res_mean_`type'_quad_`inter' hour_to_treat_`inter' if inrange(hour_to_treat_`inter',scalar(-h),scalar(-1)) & station_id == 2, color(stc3)) || /// 
			   (line res_mean_`type'_quad_`inter' hour_to_treat_`inter' if inrange(hour_to_treat_`inter',scalar(0),scalar(h)) & station_id == 2, color(stc3) legend(pos(6) col(4) off) xtitle("Hours to Treatment") ytitle("Residual Price `type' (ct/l)") ysize(5cm) xsize(8cm) scale(1.1))
			   graph save "C:\Users\frede\Desktop\Bachelor-Thesis\Two_Stage_RDiT_JTEP\revised figures_hquality\rob_poly_`type'_`inter'.gph", replace
			   graph export "C:\Users\frede\Desktop\Bachelor-Thesis\Two_Stage_RDiT_JTEP\revised figures_hquality\rob_poly_`type'_`inter'.emf", replace
		drop bin res_mean_`type' res_mean_`type'_bin
	}
	drop kernel_wgt
}

*** ######################################################################## ***
* Heterogeneous Pass-Through HHI and VERT INTERVENTION
graph set window fontface "Times New Roman"
scalar scal_tax_mean_diesel = 1/16.71 * 100 // (1/+-16.71)*100 = +-5.984440455
scalar scal_tax_mean_e5 = 1/35.16 * 100 // (1/+-35.16)*100 = +-2.844141069
scalar scal_tax_mean_e10 = 1/35.16 * 100 // (1/+-35.16)*100 = +-2.844141069
* tax reduction:
foreach i in e5 e10 diesel {
		cd "C:\Users\frede\Desktop\RDiT JTEP Results\FIRST_`i'"
		estimates use "first_`i'_hhi_vert_1"
			scalar scal_tankrabatt_1 = _b[1.tankrabatt_1]
			scalar scal_tankrabatt_1_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_`i'
		estimates use "first_`i'_hhi_vert_2"
			scalar scal_tankrabatt_2 = _b[1.tankrabatt_1]
			scalar scal_tankrabatt_2_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_`i'
		estimates use "first_`i'_hhi_vert_3"
			scalar scal_tankrabatt_3 = _b[1.tankrabatt_1]
			scalar scal_tankrabatt_3_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_`i'
			graph set window fontface "Times New Roman"
		graph use "HHI_VERT_first_`i'_per", scheme(stgcolor)
		graph export "C:\Users\frede\Desktop\RDiT JTEP Results\revised figures_hquality\HHI_VERT_first_`i'_per_final.svg", replace
}
*tax increase:
foreach i in e5 e10 diesel {
		cd "C:\Users\frede\Desktop\RDiT JTEP Results\SECOND_`i'"
		estimates use "second_`i'_hhi_vert_1"
			scalar scal_tankrabatt_1 = _b[1.tankrabatt_2]
			scalar scal_tankrabatt_1_per = abs(_b[1.tankrabatt_2])*scal_tax_mean_`i'
		estimates use "second_`i'_hhi_vert_2"
			scalar scal_tankrabatt_2 = _b[1.tankrabatt_2]
			scalar scal_tankrabatt_2_per = abs(_b[1.tankrabatt_2])*scal_tax_mean_`i'
		estimates use "second_`i'_hhi_vert_3"
			scalar scal_tankrabatt_3 = _b[1.tankrabatt_2]
			scalar scal_tankrabatt_3_per = abs(_b[1.tankrabatt_2])*scal_tax_mean_`i'
			graph set window fontface "Times New Roman"
		graph use "HHI_VERT_second_`i'_per", scheme(stgcolor)
		graph export "C:\Users\frede\Desktop\RDiT JTEP Results\revised figures_hquality\HHI_VERT_second_`i'_per_final.svg", replace
}

*** ######################################################################## ***
* Heterogeneous Pass-Through HHI Sensitivity
scalar scal_tax_mean_diesel = 1/16.71 * 100 // (1/+-16.71)*100 = +-5.984440455
scalar scal_tax_mean_e5 = 1/35.16 * 100 // (1/+-35.16)*100 = +-2.844141069
scalar scal_tax_mean_e10 = 1/35.16 * 100 // (1/+-35.16)*100 = +-2.844141069
* tax reduction:
foreach i in e5 e10 diesel {
		cd "C:\Users\frede\Desktop\RDiT JTEP Results\FIRST_`i'"
		  estimates use "E3_first_`i'_hhi_rad_3"
			scalar scal_tankrabatt_3 = _b[1.tankrabatt_1]
			scalar scal_tankrabatt_3_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_`i'
		estimates use "E4_first_`i'_hhi_rad_4" 
			scalar scal_tankrabatt_4 = _b[1.tankrabatt_1]
			scalar scal_tankrabatt_4_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_`i'
		estimates use "E5_first_`i'_hhi_rad_5"
			scalar scal_tankrabatt_5 = _b[1.tankrabatt_1]
			scalar scal_tankrabatt_5_per = abs(_b[1.tankrabatt_1])*scal_tax_mean_`i'
			graph set window fontface "Times New Roman"
		graph use "HHI_SENS_first_`i'_per", scheme(stgcolor)
		graph export "C:\Users\frede\Desktop\RDiT JTEP Results\revised figures_hquality\HHI_SENS_first_`i'_per.svg", replace
}
*tax increase:
foreach i in e5 e10 diesel  {
		cd "C:\Users\frede\Desktop\RDiT JTEP Results\SECOND_`i'"
		  estimates use "E3_second_`i'_hhi_rad_3"
			scalar scal_tankrabatt_3 = _b[1.tankrabatt_2]
			scalar scal_tankrabatt_3_per = abs(_b[1.tankrabatt_2])*scal_tax_mean_`i'
		estimates use "E4_second_`i'_hhi_rad_4" 
			scalar scal_tankrabatt_4 = _b[1.tankrabatt_2]
			scalar scal_tankrabatt_4_per = abs(_b[1.tankrabatt_2])*scal_tax_mean_`i'
		estimates use "E5_second_`i'_hhi_rad_5"
			scalar scal_tankrabatt_5 = _b[1.tankrabatt_2]
			scalar scal_tankrabatt_5_per = abs(_b[1.tankrabatt_2])*scal_tax_mean_`i'
			graph set window fontface "Times New Roman"
		graph use "HHI_SENS_second_`i'_per", scheme(stgcolor)
		graph export "C:\Users\frede\Desktop\RDiT JTEP Results\revised figures_hquality\HHI_SENS_second_`i'_per.svg", replace
}


*** ######################################################################## ***
* Margins Plot for three-way interaction HHI & Vert
* FIRST INTERVENTION
* E5
	cd "C:\Users\frede\Desktop\RDiT JTEP Results\full interaction HHI VERT"
	estimates use "first_e5_het_included"
	regress
forvalues k = 1/3 {
    forvalues j = 1/4 {
        matrix I_`k'_`j' = J(1,3,.)
        lincom 1.tankrabatt_1 + ///
               1.tankrabatt_1#`k'.brand_vert_int + ///
               1.tankrabatt_1#`j'.hhi_km_quantile + ///
               1.tankrabatt_1#`j'.hhi_km_quantile#`k'.brand_vert_int
		matrix I_`k'_`j'[1,1] = r(estimate), r(lb), r(ub)
	}
	matrix Vert_`k' = I_`k'_1 \ I_`k'_2 \ I_`k'_3 \ I_`k'_4
	matrix rown Vert_`k' = hhi1 hhi2 hhi3 hhi4
	matrix coln Vert_`k' = b ll95 ul95
	matlist Vert_`k'
}
scalar scal_tax_mean_e10 = 1/35.16 * 100 // (1/+-35.16)*100 = +-2.844141069
coefplot ///
    (matrix(Vert_1[,1]), rescale(-2.844141069) recast(connected) lpattern(dash) offset(-0.12) ci((2 3)) label(oligopolistic) color(stc1) ///
        msymbol(circle) mlwidth(medthick) msize(1.4) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop) ///
    (matrix(Vert_2[,1]), rescale(-2.844141069) recast(connected) lpattern(dash) ci((2 3)) label(other integrated) color(stc2) ///
        msymbol(triangle) mlwidth(medthick) msize(1.4) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop)  ///
    (matrix(Vert_3[,1]), rescale(-2.844141069) recast(connected) lpattern(dash) offset(0.12) ci((2 3)) label(independent) color(stc3) ///
        msymbol(square) mlwidth(medthick) msize(1.4) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop), ///
    vertical ///
    xtitle("HHI Quantile") ///
	ytitle("Tax Pass-Through E5 (%)") xtitle("HHI Quartiles") /// 
    legend(pos(6) col(3)) xlabel(1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4") /// 
	scale(1.3) ysize(5cm) xsize(8cm)
	graph export "Interact_HHI_VERT_E5_first_per.emf", replace

*E10
	cd "C:\Users\frede\Desktop\RDiT JTEP Results\full interaction HHI VERT"
	estimates use "first_e10_het_included"
	regress
forvalues k = 1/3 {
    forvalues j = 1/4 {
        matrix I_`k'_`j' = J(1,3,.)
        lincom 1.tankrabatt_1 + ///
               1.tankrabatt_1#`k'.brand_vert_int + ///
               1.tankrabatt_1#`j'.hhi_km_quantile + ///
               1.tankrabatt_1#`j'.hhi_km_quantile#`k'.brand_vert_int
		matrix I_`k'_`j'[1,1] = r(estimate), r(lb), r(ub)
	}
	matrix Vert_`k' = I_`k'_1 \ I_`k'_2 \ I_`k'_3 \ I_`k'_4
	matrix rown Vert_`k' = hhi1 hhi2 hhi3 hhi4
	matrix coln Vert_`k' = b ll95 ul95
	matlist Vert_`k'
}
scalar scal_tax_mean_e10 = 1/35.16 * 100 // (1/+-35.16)*100 = +-2.844141069
coefplot ///
    (matrix(Vert_1[,1]), rescale(-2.844141069) recast(connected) lpattern(dash) offset(-0.12) ci((2 3)) label(oligopolistic) color(stc1) ///
        msymbol(circle) mlwidth(medthick) msize(1.4) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop) ///
    (matrix(Vert_2[,1]), rescale(-2.844141069) recast(connected) lpattern(dash) ci((2 3)) label(other integrated) color(stc2) ///
        msymbol(triangle) mlwidth(medthick) msize(1.4) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop)  ///
    (matrix(Vert_3[,1]), rescale(-2.844141069) recast(connected) lpattern(dash) offset(0.12) ci((2 3)) label(independent) color(stc3) ///
        msymbol(square) mlwidth(medthick) msize(1.4) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop), ///
    vertical ///
    xtitle("HHI Quantile") ///
	ytitle("Tax Pass-Through E10 (%)") xtitle("HHI Quartiles") /// 
    legend(pos(6) col(3)) xlabel(1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4") /// 
	scale(1.3) ysize(5cm) xsize(8cm)
	graph export "Interact_HHI_VERT_E10_first_per.emf", replace

	
* DIESEL
	cd "C:\Users\frede\Desktop\RDiT JTEP Results\full interaction HHI VERT"
	estimates use "first_diesel_het_included"
	regress
forvalues k = 1/3 {
    forvalues j = 1/4 {
        matrix I_`k'_`j' = J(1,3,.)
        lincom 1.tankrabatt_1 + ///
               1.tankrabatt_1#`k'.brand_vert_int + ///
               1.tankrabatt_1#`j'.hhi_km_quantile + ///
               1.tankrabatt_1#`j'.hhi_km_quantile#`k'.brand_vert_int
		matrix I_`k'_`j'[1,1] = r(estimate), r(lb), r(ub)
	}
	matrix Vert_`k' = I_`k'_1 \ I_`k'_2 \ I_`k'_3 \ I_`k'_4
	matrix rown Vert_`k' = hhi1 hhi2 hhi3 hhi4
	matrix coln Vert_`k' = b ll95 ul95
	matlist Vert_`k'
}
scalar scal_tax_mean_diesel = 1/16.71 * 100 // (1/+-16.71)*100 = +-5.984440455
coefplot ///
    (matrix(Vert_1[,1]), rescale(-5.984440455) recast(connected) lpattern(dash) offset(-0.12) ci((2 3)) label(oligopolistic) color(stc1) ///
        msymbol(circle) mlwidth(medthick) msize(1.4) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop) ///
    (matrix(Vert_2[,1]), rescale(-5.984440455) recast(connected) lpattern(dash) ci((2 3)) label(other integrated) color(stc2) ///
        msymbol(triangle) mlwidth(medthick) msize(1.4) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop)  ///
    (matrix(Vert_3[,1]), rescale(-5.984440455) recast(connected) lpattern(dash) offset(0.12) ci((2 3)) label(independent) color(stc3) ///
        msymbol(square) mlwidth(medthick) msize(1.4) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop), ///
    vertical ///
    xtitle("HHI Quantile") ///
	ytitle("Tax Pass-Through Diesel (%)") xtitle("HHI Quartiles") /// 
    legend(pos(6) col(3)) xlabel(1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4") /// 
	scale(1.3) ysize(5cm) xsize(8cm)
	graph export "Interact_HHI_VERT_Diesel_first_per.emf", replace


	
* SECOND INTERVENTION
* E5 
cd "C:\Users\frede\Desktop\RDiT JTEP Results\full interaction HHI VERT"
	estimates use "second_e5_het_included"
	regress
forvalues k = 1/3 {
    forvalues j = 1/4 {
        matrix I_`k'_`j' = J(1,3,.)
        lincom 1.tankrabatt_2 + ///
               1.tankrabatt_2#`k'.brand_vert_int + ///
               1.tankrabatt_2#`j'.hhi_km_quantile + ///
               1.tankrabatt_2#`j'.hhi_km_quantile#`k'.brand_vert_int
		matrix I_`k'_`j'[1,1] = r(estimate), r(lb), r(ub)
	}
	matrix Vert_`k' = I_`k'_1 \ I_`k'_2 \ I_`k'_3 \ I_`k'_4
	matrix rown Vert_`k' = hhi1 hhi2 hhi3 hhi4
	matrix coln Vert_`k' = b ll95 ul95
	matlist Vert_`k'
}
scalar scal_tax_mean_e10 = 1/35.16 * 100 // (1/+-35.16)*100 = +-2.844141069
coefplot ///
    (matrix(Vert_1[,1]), rescale(2.844141069) recast(connected) lpattern(dash) offset(-0.12) ci((2 3)) label(oligopolistic) color(stc1) ///
        msymbol(circle) mlwidth(medthick) msize(1.4) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop) ///
    (matrix(Vert_2[,1]), rescale(2.844141069) recast(connected) lpattern(dash) ci((2 3)) label(other integrated) color(stc2) ///
        msymbol(triangle) mlwidth(medthick) msize(1.4) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop)  ///
    (matrix(Vert_3[,1]), rescale(2.844141069) recast(connected) lpattern(dash) offset(0.12) ci((2 3)) label(independent) color(stc3) ///
        msymbol(square) mlwidth(medthick) msize(1.4) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop), ///
    vertical ///
    xtitle("HHI Quantile") ///
	ytitle("Tax Pass-Through E5 (%)") xtitle("HHI Quartiles") /// 
    legend(pos(6) col(3)) xlabel(1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4") /// 
	scale(1.3) ysize(5cm) xsize(8cm)
	graph export "Interact_HHI_VERT_E5_second_per.emf", replace
* E10
cd "C:\Users\frede\Desktop\RDiT JTEP Results\full interaction HHI VERT"
	estimates use "second_e10_het_included"
	regress
forvalues k = 1/3 {
    forvalues j = 1/4 {
        matrix I_`k'_`j' = J(1,3,.)
        lincom 1.tankrabatt_2 + ///
               1.tankrabatt_2#`k'.brand_vert_int + ///
               1.tankrabatt_2#`j'.hhi_km_quantile + ///
               1.tankrabatt_2#`j'.hhi_km_quantile#`k'.brand_vert_int
		matrix I_`k'_`j'[1,1] = r(estimate), r(lb), r(ub)
	}
	matrix Vert_`k' = I_`k'_1 \ I_`k'_2 \ I_`k'_3 \ I_`k'_4
	matrix rown Vert_`k' = hhi1 hhi2 hhi3 hhi4
	matrix coln Vert_`k' = b ll95 ul95
	matlist Vert_`k'
}
scalar scal_tax_mean_e10 = 1/35.16 * 100 // (1/+-35.16)*100 = +-2.844141069
coefplot ///
    (matrix(Vert_1[,1]), rescale(2.844141069) recast(connected) lpattern(dash) offset(-0.12) ci((2 3)) label(oligopolistic) color(stc1) ///
        msymbol(circle) mlwidth(medthick) msize(1.4) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop) ///
    (matrix(Vert_2[,1]), rescale(2.844141069) recast(connected) lpattern(dash) ci((2 3)) label(other integrated) color(stc2) ///
        msymbol(triangle) mlwidth(medthick) msize(1.4) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop)  ///
    (matrix(Vert_3[,1]), rescale(2.844141069) recast(connected) lpattern(dash) offset(0.12) ci((2 3)) label(independent) color(stc3) ///
        msymbol(square) mlwidth(medthick) msize(1.4) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop), ///
    vertical ///
    xtitle("HHI Quantile") ///
	ytitle("Tax Pass-Through E10 (%)") xtitle("HHI Quartiles") /// 
    legend(pos(6) col(3)) xlabel(1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4") /// 
	scale(1.3) ysize(5cm) xsize(8cm)
	graph export "Interact_HHI_VERT_E10_second_per.emf", replace
* DIESEL
cd "C:\Users\frede\Desktop\RDiT JTEP Results\full interaction HHI VERT"
	estimates use "second_diesel_het_included"
	regress
forvalues k = 1/3 {
    forvalues j = 1/4 {
        matrix I_`k'_`j' = J(1,3,.)
        lincom 1.tankrabatt_2 + ///
               1.tankrabatt_2#`k'.brand_vert_int + ///
               1.tankrabatt_2#`j'.hhi_km_quantile + ///
               1.tankrabatt_2#`j'.hhi_km_quantile#`k'.brand_vert_int
		matrix I_`k'_`j'[1,1] = r(estimate), r(lb), r(ub)
	}
	matrix Vert_`k' = I_`k'_1 \ I_`k'_2 \ I_`k'_3 \ I_`k'_4
	matrix rown Vert_`k' = hhi1 hhi2 hhi3 hhi4
	matrix coln Vert_`k' = b ll95 ul95
	matlist Vert_`k'
}
scalar scal_tax_mean_diesel = 1/16.71 * 100 // (1/+-16.71)*100 = +-5.984440455
coefplot ///
    (matrix(Vert_1[,1]), rescale(5.984440455) recast(connected) lpattern(dash) offset(-0.12) ci((2 3)) label(oligopolistic) color(stc1) ///
        msymbol(circle) mlwidth(medthick) msize(1.4) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop) ///
    (matrix(Vert_2[,1]), rescale(5.984440455) recast(connected) lpattern(dash) ci((2 3)) label(other integrated) color(stc2) ///
        msymbol(triangle) mlwidth(medthick) msize(1.4) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop)  ///
    (matrix(Vert_3[,1]), rescale(5.984440455) recast(connected) lpattern(dash) offset(0.12) ci((2 3)) label(independent) color(stc3) ///
        msymbol(square) mlwidth(medthick) msize(1.4) finten(40) ciopts(recast(rcap) lwidth(medthick)) citop), ///
    vertical ///
    xtitle("HHI Quantile") ///
	ytitle("Tax Pass-Through Diesel (%)") xtitle("HHI Quartiles") /// 
    legend(pos(6) col(3)) xlabel(1 "Q1" 2 "Q2" 3 "Q3" 4 "Q4") /// 
	scale(1.3) ysize(5cm) xsize(8cm)
	graph export "Interact_HHI_VERT_Diesel_second_per.emf", replace






*###############################################################################
* Table - Descriptive Statistics before and after 1. an 2. Intervention
*###############################################################################


*###############################################################################
* First Intervention - ROAD
*###############################################################################
cd "C:\Users\frede\Desktop\Bachelor-Thesis\Daten\tankstellenpreise"
use "roadside_price_2022_0309_openandassume_revised.dta", clear

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
scalar ba = 336 // Bandwidth 120 hours
scalar db = -da
scalar bb = -ba

*-------------------------------------------------------------------------------
* descriptives BEFORE and AFTER
est clear 
eststo: estpost gtabstat mean_e5 ///
				 mean_e10 /// 
				 mean_diesel /// 
				 BRENT_FOB /// 
				 usd_to_euro /// 
				 kfz_menge_B ///
				 temperature /// 
				 precipitation if inrange(hour_to_treat_1, scalar(bb),0), statistics(mean sd) columns(statistics)
eststo: estpost gtabstat mean_e5 ///
				 mean_e10 /// 
				 mean_diesel /// 
				 BRENT_FOB /// 
				 usd_to_euro /// 
				 kfz_menge_B ///
				 temperature /// 
				 precipitation if inrange(hour_to_treat_1, 0,scalar(ba)), statistics(mean sd) columns(statistics) 
distinct station_id if inrange(hour_to_treat_1, scalar(bb),0)
distinct station_id if inrange(hour_to_treat_1, 0,scalar(ba))


*###############################################################################
* First Intervention - HIGHWAY
*###############################################################################
cd "C:\Users\frede\Desktop\Bachelor-Thesis\Daten\tankstellenpreise"
use "roadside_price_2022_0309_openandassume_high_revised.dta", clear
gen id_new = station_id
xtset station_id time
keep if inrange(time,9961,11809) // 11 Weeks (10 Weeks pre treatment, 1 Week with treatment)
gen lnBRENT_FOB = ln(BRENT_FOB)
bysort station_id (time): gen L24_lnBRENT_FOB = L24.lnBRENT_FOB
bysort station_id (time): gen L24_usd_to_euro = L24.usd_to_euro
bysort station_id (time): gen L1_kfz_menge_A = L1.kfz_menge_A
*encode NUTS, gen(NUTS3)

*** RDiT Indicators ***
gen hour_to_treat_1 = time
replace hour_to_treat_1 = (time- 11641) // Treatment June 1st 2022
gen tankrabatt_1 = 0
replace tankrabatt_1 = 1 if inrange(month,6,8)
scalar da = 9 // Donut 9 hours
scalar ba = 336 // Bandwidth 120 hours
scalar db = -da
scalar bb = -ba
global var_seas i.hour i.week_days
global var_sup lnBRENT_FOB L24_lnBRENT_FOB usd_to_euro L24_usd_to_euro	 
global var_dem kfz_menge_A L1_kfz_menge_A temperature precipitation holiday_state day_before_holiday_state

*-------------------------------------------------------------------------------
* descriptives BEFORE and AFTER
eststo: estpost gtabstat mean_e5 ///
				 mean_e10 /// 
				 mean_diesel /// 
				 BRENT_FOB /// 
				 usd_to_euro /// 
				 kfz_menge_A ///
				 temperature /// 
				 precipitation if inrange(hour_to_treat_1, scalar(bb),0), statistics(mean sd) columns(statistics)
eststo: estpost gtabstat mean_e5 ///
				 mean_e10 /// 
				 mean_diesel /// 
				 BRENT_FOB /// 
				 usd_to_euro /// 
				 kfz_menge_A ///
				 temperature /// 
				 precipitation if inrange(hour_to_treat_1, 0,scalar(ba)), statistics(mean sd) columns(statistics) 
distinct station_id if inrange(hour_to_treat_1, scalar(bb),0)
distinct station_id if inrange(hour_to_treat_1, 0,scalar(ba))
				 
*###############################################################################
* Second Intervention - ROAD
*###############################################################################
			 
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
keep if inrange(time,12169,14017) // 11 Weeks (10 Weeks pre treatment, 1 Week with treatment)
gen lnBRENT_FOB = ln(BRENT_FOB)
bysort station_id (time): gen L24_lnBRENT_FOB = L24.lnBRENT_FOB
bysort station_id (time): gen L24_usd_to_euro = L24.usd_to_euro
bysort station_id (time): gen L1_kfz_menge_B = L1.kfz_menge_B
*encode NUTS, gen(NUTS3)

*** RDiT Indicators ***
gen hour_to_treat_1 = time
replace hour_to_treat_1 = (time- 13849) // Treatment September 1st 2022
gen tankrabatt_1 = 1
replace tankrabatt_1 = 0 if inrange(month,6,8)
scalar da = 9 // Donut 9 hours
scalar ba = 336 // Bandwidth 120 hours
scalar db = -da
scalar bb = -ba
global var_seas i.hour i.week_days
global var_sup lnBRENT_FOB L24_lnBRENT_FOB usd_to_euro L24_usd_to_euro	 
global var_dem kfz_menge_B L1_kfz_menge_B temperature precipitation holiday_state day_before_holiday_state

*-------------------------------------------------------------------------------
* descriptives BEFORE and AFTER
eststo: estpost gtabstat mean_e5 ///
				 mean_e10 /// 
				 mean_diesel /// 
				 BRENT_FOB /// 
				 usd_to_euro /// 
				 kfz_menge_B ///
				 temperature /// 
				 precipitation if inrange(hour_to_treat_1, scalar(bb),0), statistics(mean sd) columns(statistics)
eststo: estpost gtabstat mean_e5 ///
				 mean_e10 /// 
				 mean_diesel /// 
				 BRENT_FOB /// 
				 usd_to_euro /// 
				 kfz_menge_B ///
				 temperature /// 
				 precipitation if inrange(hour_to_treat_1, 0,scalar(ba)), statistics(mean sd) columns(statistics) 
distinct station_id if inrange(hour_to_treat_1, scalar(bb),0)
distinct station_id if inrange(hour_to_treat_1, 0,scalar(ba))

*###############################################################################
* Second Intervention - HIGHWAY
*###############################################################################

cd "C:\Users\frede\Desktop\Bachelor-Thesis\Daten\tankstellenpreise"
use "roadside_price_2022_0309_openandassume_high_revised.dta", clear
gen id_new = station_id
xtset station_id time
keep if inrange(time,12169,14017) // 11 Weeks (10 Weeks pre treatment, 1 Week with treatment)
gen lnBRENT_FOB = ln(BRENT_FOB)
bysort station_id (time): gen L24_lnBRENT_FOB = L24.lnBRENT_FOB
bysort station_id (time): gen L24_usd_to_euro = L24.usd_to_euro
bysort station_id (time): gen L1_kfz_menge_A = L1.kfz_menge_A
*encode NUTS, gen(NUTS3)

*** RDiT Indicators ***
gen hour_to_treat_1 = time
replace hour_to_treat_1 = (time- 13849) // Treatment September 1st 2022
gen tankrabatt_1 = 1
replace tankrabatt_1 = 0 if inrange(month,6,8)
scalar da = 9 // Donut 9 hours
scalar ba = 336 // Bandwidth 120 hours
scalar db = -da
scalar bb = -ba
global var_seas i.hour i.week_days
global var_sup lnBRENT_FOB L24_lnBRENT_FOB usd_to_euro L24_usd_to_euro	 
global var_dem kfz_menge_A L1_kfz_menge_A temperature precipitation holiday_state day_before_holiday_state

*-------------------------------------------------------------------------------
* descriptives BEFORE and AFTER
eststo: estpost gtabstat mean_e5 ///
				 mean_e10 /// 
				 mean_diesel /// 
				 BRENT_FOB /// 
				 usd_to_euro /// 
				 kfz_menge_A ///
				 temperature /// 
				 precipitation if inrange(hour_to_treat_1, scalar(bb),0), statistics(mean sd) columns(statistics)
eststo: estpost gtabstat mean_e5 ///
				 mean_e10 /// 
				 mean_diesel /// 
				 BRENT_FOB /// 
				 usd_to_euro /// 
				 kfz_menge_A ///
				 temperature /// 
				 precipitation if inrange(hour_to_treat_1, 0,scalar(ba)), statistics(mean sd) columns(statistics) 
distinct station_id if inrange(hour_to_treat_1, scalar(bb),0)
distinct station_id if inrange(hour_to_treat_1, 0,scalar(ba))

esttab est* using "C:\Users\frede\Desktop\Bachelor-Thesis\Daten\tankstellenpreise\descriptives_before_after_hour_rev.rtf", replace  main(mean 2) aux(sd 2) nostar unstack nonote label 



