library(lubridate)

flux_path <- './R/FLUXNET2015/'
files <- list.files(flux_path)
sites <- read.csv('./R/site.info.csv')



transfer <- function(flux_HH_yearly){
  # "DTIME", 
  old_list <- c("Year", "DoY", "Time", "NEE_f", "NEE_f_delta", "GPP_f", "GPP_f_delta", "Reco", 
                "NEE_GPP_qc", "NEE_GPP_qcOK", "LE_f", "LE_fqc", "LE_fqcOK", "H_f", "H_fqc", "H_fqcOK", 
                "G_f", "G_fqc", "G_fqcOK", "Ta_f", "Ta_fqc", "Ta_fqcOK", "Ts1_f", "Ts1_fqc", "Ts1_fqcOK", 
                "Ts2_f", "Ts2_fqc", "Ts2_fqcOK", "VPD_f", "VPD_fqc", "VPD_fqcOK", "Precip_f", "Precip_fqc", 
                "Precip_fqcOK", "SWC1_f", "SWC1_fqc", "SWC1_fqcOK", "SWC2_f", "SWC2_fqc", "SWC2_fqcOK", "WS_f", 
                "WS_fqc", "WS_fqcOK", "Rg_f", "Rg_fqc", "Rg_fqcOK", "PPFD_f", "PPFD_fqc", "PPFD_fqcOK", "Rn_f", 
                "Rn_fqc", "Rn_fqcOK", "Rg_pot", "Rd", "Rr", "PPFDbc", "PPFDd", "PPFDr", "FAPAR", "LWin", "LWout", 
                "SWin", "SWout", "WD", "ustar", "ZL", "Rh", "CO2", "H2O", "H2Ostor1", "H2Ostor2", "NEE_f_unc", 
                "NEE_fqc_unc", "NEE_fqcOK_unc", "LE_f_unc", "LE_fqc_unc", "LE_fqcOK_unc", "H_f_unc", "H_fqc_unc", 
                "H_fqcOK_unc", "Reco_E0_100", "Reco_E0_200", "Reco_E0_300", "wbal_clim", "wbal_act", "wdef_cum", "Epot_viaLE_H", 
                "EpotPT_viaLE_H", "Epot_viaRg", "Epot_viaRn", "Epot_f", "Epot_flag", "gsurf_viaRg", "gsurf_viaRn", "gsurf_viaLE_H", 
                "gsurf_f", "gsurf_flag", "Drain", "EpotPT_viaRn", "H2Ostor1_hWHC", "H2Ostor2_hWHC", "Drain_hWHC")
  
  
  
  # Year,DoY,Time,DTIME,
  flux_HH_yearly$Year = year(ymd_hm(as.character(flux_HH_yearly$TIMESTAMP_START)))
  flux_HH_yearly$DoY = yday(ymd_hm(as.character(flux_HH_yearly$TIMESTAMP_START)))
  flux_HH_yearly$Time = minute(ymd_hm(as.character(flux_HH_yearly$TIMESTAMP_START)))/60 + hour(ymd_hm(as.character(flux_HH_yearly$TIMESTAMP_START)))
  flux_HH_yearly$DTIME = -9999
  
  # NEE_f,NEE_f_delta,GPP_f,GPP_f_delta,Reco,NEE_GPP_qc,NEE_GPP_qcOK,
  
  flux_HH_yearly$NEE_f = flux_HH_yearly$NEE_VUT_REF
  flux_HH_yearly$NEE_f_delta = -9999
  flux_HH_yearly$GPP_f = flux_HH_yearly$GPP_DT_VUT_REF # NT is not right
  flux_HH_yearly$GPP_f_delta = -9999
  flux_HH_yearly$Reco = flux_HH_yearly$RECO_DT_VUT_REF
  flux_HH_yearly$NEE_GPP_qc = -9999
  flux_HH_yearly$NEE_GPP_qcOK = -9999
  
  # LE_f,LE_fqc,LE_fqcOK,H_f,H_fqc,H_fqcOK,G_f,G_fqc,G_fqcOK,
  
  flux_HH_yearly$LE_f = flux_HH_yearly$LE_F_MDS
  flux_HH_yearly$LE_fqc = -9999 #flux_HH_yearly$LE_F_MDS_QC
  flux_HH_yearly$LE_fqcOK = -9999
  flux_HH_yearly$H_f = flux_HH_yearly$H_F_MDS
  flux_HH_yearly$H_fqc = -9999 #flux_HH_yearly$H_F_MDS_QC
  flux_HH_yearly$H_fqcOK = -9999
  flux_HH_yearly$G_f = flux_HH_yearly$G_F_MDS
  flux_HH_yearly$G_fqc = -9999 #flux_HH_yearly$G_F_MDS_QC
  flux_HH_yearly$G_fqcOK = -9999
  
  # Ta_f,Ta_fqc,Ta_fqcOK,Ts1_f,Ts1_fqc,Ts1_fqcOK,Ts2_f,Ts2_fqc,Ts2_fqcOK
  
  flux_HH_yearly$Ta_f = flux_HH_yearly$TA_F
  flux_HH_yearly$Ta_fqc = -9999 #flux_HH_yearly$TA_F_QC
  flux_HH_yearly$Ta_fqcOK = -9999
  flux_HH_yearly$Ts1_f = flux_HH_yearly$TS_F_MDS_1
  flux_HH_yearly$Ts1_fqc = -9999 #fflux_HH_yearly$TS_F_MDS_1_QC
  flux_HH_yearly$Ts1_fqcOK = -9999
  flux_HH_yearly$Ts2_f = flux_HH_yearly$TS_F_MDS_2
  flux_HH_yearly$Ts2_fqc = -9999 #flux_HH_yearly$TS_F_MDS_2_QC
  flux_HH_yearly$Ts2_fqcOK = -9999
  
  
  # VPD_f,VPD_fqc,VPD_fqcOK,Precip_f,Precip_fqc,Precip_fqcOK,
  
  flux_HH_yearly$VPD_f = flux_HH_yearly$VPD_F
  flux_HH_yearly$VPD_fqc = -9999 #fflux_HH_yearly$VPD_F_QC
  flux_HH_yearly$VPD_fqcOK = -9999
  flux_HH_yearly$Precip_f = flux_HH_yearly$P_F
  flux_HH_yearly$Precip_fqc = -9999 #fflux_HH_yearly$P_F_QC
  flux_HH_yearly$Precip_fqcOK = -9999
  
  
  # SWC1_f,SWC1_fqc,SWC1_fqcOK,SWC2_f,SWC2_fqc,SWC2_fqcOK,
  flux_HH_yearly$SWC1_f = flux_HH_yearly$SWC_F_MDS_1
  flux_HH_yearly$SWC1_fqc = -9999 #fflux_HH_yearly$SWC_F_MDS_1_QC
  flux_HH_yearly$SWC1_fqcOK = -9999
  flux_HH_yearly$SWC2_f = flux_HH_yearly$SWC_F_MDS_2
  flux_HH_yearly$SWC2_fqc = -9999 #fflux_HH_yearly$SWC_F_MDS_2_QC
  flux_HH_yearly$SWC2_fqcOK = -9999 
  
  
  # WS_f,WS_fqc,WS_fqcOK
  flux_HH_yearly$WS_f = flux_HH_yearly$WS_F
  flux_HH_yearly$WS_fqc = -9999 #fflux_HH_yearly$WS_F_QC
  flux_HH_yearly$WS_fqcOK = -9999
  
  
  # Rg_f,Rg_fqc,Rg_fqcOK,PPFD_f,PPFD_fqc,PPFD_fqcOK,Rn_f,Rn_fqc,Rn_fqcOK,Rg_pot,Rd,Rr,
  flux_HH_yearly$Rg_f = flux_HH_yearly$LW_IN_F + flux_HH_yearly$SW_IN_F
  flux_HH_yearly$Rg_fqc = -9999
  flux_HH_yearly$Rg_fqcOK = -9999
  flux_HH_yearly$PPFD_f = flux_HH_yearly$PPFD_IN
  flux_HH_yearly$PPFD_fqc = -9999 #flux_HH_yearly$PPFD_IN_QC
  flux_HH_yearly$PPFD_fqcOK = -9999
  flux_HH_yearly$Rn_f = flux_HH_yearly$NETRAD # ? Global Radiation == Net radiation ??
  flux_HH_yearly$Rn_fqc = -9999# flux_HH_yearly$NETRAD_QC
  flux_HH_yearly$Rn_fqcOK = -9999
  flux_HH_yearly$Rg_pot = -9999
  flux_HH_yearly$Rd = -9999
  flux_HH_yearly$Rr = -9999
  
  # PPFDbc,PPFDd,PPFDr,FAPAR,LWin,LWout,SWin,SWout,WD,
  # ,,,,,,,,,
  flux_HH_yearly$PPFDbc = -9999
  flux_HH_yearly$PPFDd = flux_HH_yearly$PPFD_DIF
  flux_HH_yearly$PPFDr = -9999
  flux_HH_yearly$FAPAR = -9999
  flux_HH_yearly$LWin = flux_HH_yearly$LW_IN_F
  flux_HH_yearly$LWout = flux_HH_yearly$LW_OUT
  flux_HH_yearly$SWin = flux_HH_yearly$SW_IN_F
  flux_HH_yearly$SWout = flux_HH_yearly$SW_OUT
  flux_HH_yearly$WD = flux_HH_yearly$WD
  
  # ustar,ZL,Rh,CO2,H2O,H2Ostor1,H2Ostor2,NEE_f_unc,NEE_fqc_unc,NEE_fqcOK_unc,
  # ,,,,,,,,,,
  flux_HH_yearly$ustar = flux_HH_yearly$USTAR
  flux_HH_yearly$ZL = -9999
  flux_HH_yearly$Rh = flux_HH_yearly$RH
  flux_HH_yearly$CO2 = flux_HH_yearly$CO2_F_MDS
  flux_HH_yearly$H2O = -9999
  flux_HH_yearly$H2Ostor1 = -9999
  flux_HH_yearly$H2Ostor2 = -9999
  
  flux_HH_yearly$NEE_f_unc = -9999
  flux_HH_yearly$NEE_fqc_unc = -9999
  flux_HH_yearly$NEE_fqcOK_unc = -9999
  
  
  # LE_f_unc,LE_fqc_unc,LE_fqcOK_unc,H_f_unc,H_fqc_unc,H_fqcOK_unc, Reco_E0_100,Reco_E0_200,Reco_E0_300,
  flux_HH_yearly$LE_f_unc = -9999
  flux_HH_yearly$LE_fqc_unc = -9999
  flux_HH_yearly$LE_fqcOK_unc = -9999
  flux_HH_yearly$H_f_unc = -9999
  flux_HH_yearly$H_fqc_unc = -9999
  flux_HH_yearly$H_fqcOK_unc = -9999
  
  flux_HH_yearly$Reco_E0_100 = -9999
  flux_HH_yearly$Reco_E0_200 = -9999
  flux_HH_yearly$Reco_E0_300 = -9999
  
  # wbal_clim,wbal_act,wdef_cum,Epot_viaLE_H,EpotPT_viaLE_H,Epot_viaRg,Epot_viaRn,Epot_f,Epot_flag,
  flux_HH_yearly$wbal_clim = -9999
  flux_HH_yearly$wbal_act = -9999
  flux_HH_yearly$wdef_cum = -9999
  flux_HH_yearly$Epot_viaLE_H = -9999
  flux_HH_yearly$EpotPT_viaLE_H = -9999
  flux_HH_yearly$Epot_viaRg = -9999
  flux_HH_yearly$Epot_viaRn = -9999
  flux_HH_yearly$Epot_f = -9999
  flux_HH_yearly$Epot_flag = -9999
  
  # gsurf_viaRg,gsurf_viaRn,gsurf_viaLE_H,gsurf_f,gsurf_flag,Drain,EpotPT_viaRn,H2Ostor1_hWHC,H2Ostor2_hWHC,Drain_hWHC
  
  flux_HH_yearly$gsurf_viaRg = -9999
  flux_HH_yearly$gsurf_viaRn = -9999
  flux_HH_yearly$gsurf_viaLE_H = -9999
  flux_HH_yearly$gsurf_f = -9999
  flux_HH_yearly$gsurf_flag = -9999
  flux_HH_yearly$Drain = -9999
  flux_HH_yearly$EpotPT_viaRn = -9999
  flux_HH_yearly$H2Ostor1_hWHC = -9999
  flux_HH_yearly$H2Ostor2_hWHC = -9999
  flux_HH_yearly$Drain_hWHC = -9999
  
  old_list <- old_list[old_list %in% colnames(flux_HH_yearly)]
  flux_HH_old <- flux_HH_yearly[old_list]
  colnames(flux_HH_yearly)
  return(flux_HH_old)
}


site.info = data.frame()
site.info.full = data.frame()

for (filename in files) {
  
  flux_HH = read.csv(paste0(flux_path,filename))
  flux_HH$Year_filter <- year(ymd_hm(as.character(flux_HH$TIMESTAMP_START)))
  site_ID <- gsub("FLX_(.*)_FLUXNET.*", "\\1", filename)
  year_start <- min(flux_HH$Year_filter,na.rm = TRUE)
  year_end <-  max(flux_HH$Year_filter,na.rm = TRUE)
  # year_i =2018
  year_i = year_start
  
  while (year_i <= year_end) {
    print(paste0(site_ID,year_i))
    flux_HH_yearly <- subset(flux_HH,Year_filter == year_i)
    if (nrow(flux_HH_yearly) == 17520 | nrow(flux_HH_yearly) == 17568){
      yearly_output_old <- transfer(flux_HH_yearly)
      
      # for one year simulation
      write.csv(yearly_output_old,paste0('./Inputs/',site_ID,year_i,".",year_i,".synth.hourly.allvars.csv"), row.names=F)
      
      # for all years simulation
      write.csv(yearly_output_old,paste0('./Inputs/Inputs/',site_ID,".",year_i,".synth.hourly.allvars.csv"), row.names=F)
      
      
      # save yearly sites information
      one_site_info <- data.frame(
        Site = paste0(site_ID,year_i),
        FirstY = year_i,
        LastY = year_i,
        Lat = sites$Lat[which(sites$Site==site_ID)],
        Lon = sites$Lon[which(sites$Site==site_ID)],
        UTCtime = sites$UTCtime[which(sites$Site==site_ID)],
        timeres = 0.5
      )
      
      site.info <- rbind(site.info,one_site_info)
      
    }
    year_i = year_i +1
  }
  
  # save sites information at once

  one_site_info.full <- data.frame(
    Site = site_ID,
    FirstY = year_start,
    LastY = year_end,
    Lat = sites$Lat[which(sites$Site==site_ID)],
    Lon = sites$Lon[which(sites$Site==site_ID)],
    UTCtime = sites$UTCtime[which(sites$Site==site_ID)],
    timeres = 0.5
  )
  site.info.full <- rbind(site.info.full,one_site_info.full)
}

write.table(site.info, "./Inputs/siteinfo.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(site.info.full, "./Inputs/Inputs/siteinfo.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


# read climate =====================================
#---------------------------------------------------
# Tair    : Kelvin
# PSurf   : Pa
# Qair    : kg/kg
# Wind_E  : m s-1
# Wind_N  : m s-1
# Rainf   : kg m-2 s-1
# Snowf   : kg m-2 s-1
# SWdown  : W m-2
# LWdown  : W m-2
#---------------------------------------------------

# read observed =====================================
#---------------------------------------------------
# NEE_f	umolCO2 m-2 s-1		Net Ecosystem Exchange
# GPP_f	umolCO2 m-2 s-1		Gross ecosystem production
# Reco	umolCO2 m-2 s-1		Ecosystem Respiration calculated after Reichstein et al. (2005)
# LE_f	W m-2		        Latent heat flux
# H_f	W m-2		        Sensible heat flux
# G_f	W m-2		        Soil heat flux
# Ts1_f	deC		        Soil temperature (upper layer)
# Ts2_f	degC		        Soil temperature (lower layer)
# SWC1_f	%		        Soil water content (upper layer)
# SWC2_f	%		        Soil water content (lower layer)
# PPFD_f	umol m-2 s-1		Photosynthetic Photon Flux Density
# Rn_f	W m-2		        Net radiation
# FAPAR	%		        Fraction of Absorbed Photosynthetically Active Radiation 
# Epot_f	mm h-1		        Potential evapotranspiration after Penman according to priority: take Epot_viaRn, if not available, take Epot_viaLE_H, else Epot_viaRg 
# gsurf_f	mmol m-2 s-1		canopy conductance via inversion of Penman-Monteith,  according to priority: take gsurf_viaRn, if not available, take gsurf_viaLE_H, else gsurf_viaRg; filtered on hourly basis: take only if Rg > 150 W m-2 and gsurf_f within [0,2000 mmol m-2 s-1]
#---------------------------------------------------

# read weather =====================================
#---------------------------------------------------
# Ta_f     : Air Temperature (C)
# Precip_f : Precipitation (mm)
# Rg_f     : Global Radiation (W m-2)
# VPD_f    : Vapour Pressure Deficit (hPa)
# WS_f     : Wind horizontal speed (m s-1)
# LWin     : Incomming Longwave Radiation (W m-2)
#---------------------------------------------------