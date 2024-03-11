import os
import numpy as np
import pandas as pd
from datetime import datetime
from dlm_functions import computeAnormaly,forwardFilteringM, Model,PlotEWS

import matplotlib.pyplot as plt

fill_value = -9999


Path = './Inputs/Fluxnet2015'
files = os.listdir(Path)




file = '/Users/quan/projects/gitprojects/fluxnet/Data/FLX_GF-Guy_FLUXNET2015_FULLSET_2004-2014_2-4/FLX_GF-Guy_FLUXNET2015_FULLSET_DD_2004-2014_2-4.csv'

variable_list = ['TA_F', 'VPD_F', 'P_F', 'NETRAD']
GPP = pd.read_csv(file,sep = ',')[['TIMESTAMP','GPP_DT_VUT_REF', 'RECO_DT_VUT_REF','NEE_VUT_REF','NEE_VUT_REF_QC'] + variable_list]
GPP[GPP == -9999] = np.nan
GPP['datetime'] = pd.to_datetime(GPP['TIMESTAMP'], format='%Y%m%d')
GPP = GPP[GPP['datetime']<= '2015-01-01']
GPP_monthly = GPP.groupby([GPP['datetime'].dt.year, GPP['datetime'].dt.month]).mean()
N = GPP_monthly['NEE_VUT_REF'].to_numpy();N

CLM = GPP_monthly[variable_list].to_numpy();CLM

# daily averages of climate conditions in the same order
AvgCLM = GPP.groupby([GPP['datetime'].dt.month,GPP['datetime'].dt.day]).mean()[variable_list].to_numpy();AvgCLM
date0 = GPP['datetime'][0] # the data of first NDVI obervation
GPP 
# GPP_DT_VUT_REF	TA_F	VPD_F	P_F	NETRAD	datetime
# fluxnet2015 11.27800	25.648	3.872	9.8	105.013750	2004-01-01
# icos: 11.26490	25.648	3.232	9.8	105.013750	2004-01-01


# compute climate anomaly within each interval of two NDVI observations
anCLM = computeAnormaly(CLM,AvgCLM,date0) # 计算气象因素的异常值,16天的数据与对应期间年平均数据的差值


Y = N[1:]-np.nanmean(N) # 不理解为什么不是从N[0:]开始, 好像是要用当期的环境因素预测下一期的NDVI
# use two seasonal harmonic components # 使用两个季节性谐波分量
rseas = [1,2] 

# include lag-1 centerred NDVI and precipitation in the regression module #只包含上起的NDVI和降水
X = np.column_stack((N[:-1]-np.nanmean(N),anCLM[:-1,0]))

# set up model and run forward filtering
delta = 0.98

M = Model(Y,X,rseas,delta) # 构建模型
FF = forwardFilteringM(M,period=12) # 求解模型

