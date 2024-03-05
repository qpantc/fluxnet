# Important meteorological drivers for model
## ED2 【温度，降水，湿度，风速，风向，短波辐射】
需要用小时的数据作为驱动
tmp    -  Air temperature                  [        K] !") fluxnet_data$TA_F + 273.15 # unit: K
pres   -  Pressure                         [       Pa] !") fluxnet_data$PA_F*1000 #unit: Pa
sh     -  Specific humidity                [    kg/kg] !") absolute_humidity(fluxnet_data\$PA_F,fluxnet_data\$TA_F,fluxnet_data$RH) #unit: "kgH2O kgAir⁻¹"
ugrd   -  Zonal wind                       [      m/s] !") calculated by fluxnet_data\$WS_F and fluxnet_data\$WD
vgrd   -  Zonal wind                       [      m/s] !") calculated by fluxnet_data\$WS_F and fluxnet_data\$WD
prate  -  Precipitation rate               [  kg/m2/s] !") rainfall_kg_m2_s1(fluxnet_data$P_F) # unit: kg/m2*s1
dlwrf  -  Downward long wave radiation     [     W/m2] !") rad.in=datum\$rshort.in <- fluxnet_data\$LW_IN_F 
nbdsf  -  Near-IR beam radiation           [     W/m2] !") rshort.bdown()根据短波辐射，气压和角度计算而来
nddsf  -  Near-IR diffuse radiation        [     W/m2] !") rshort.bdown()根据短波辐射，气压和角度计算而来
vbdsf  -  Visible beam radiation           [     W/m2] !") rshort.bdown()根据短波辐射，气压和角度计算而来
vddsf  -  Visible beam radiation           [     W/m2] !") rshort.bdown()根据短波辐射，气压和角度计算而来
co2**  -  CO2 mixing ratio                 [ umol/mol] !")


## P-model 【温度，海拔，光子通量，VPD，光合有效辐射】
Tc.deg_C <- 25 # temperature in degreeC #均温
elv.m <- seq(-5,3000,1) # elevation in meter #海拔
PPFD.mol_m2 <- 1000 # PAR in mol photon m-2 month-1 每月光量子的通量
VPD.kPa <- 1 # Vapour pressure deficit in kPa
fAPAR <- 1 # factor absorbed PAR # 表示光合有效辐射，与LAI有关，当为1时全部被吸收
CO2.ppm <- 400 #ppm 二氧化碳浓度


## Model drivers: ORCHIDEE