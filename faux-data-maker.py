import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import math
import datetime
import time
import sys
import seaborn as sns
import os 
from scipy import interpolate
from scipy.interpolate import UnivariateSpline, interp1d,LSQUnivariateSpline, CubicSpline
import scipy.optimize as optimization
from numpy import linspace,exp
import matplotlib.colors as mcolors
from scipy.optimize import least_squares, leastsq
#from scipy import signal
import scipy.signal as signal
import statistics
from photutils.utils import ShepardIDWInterpolator as idw


"""
INPUT PARAMETERS
"""

# WHAT IS THE NAME OF THE FILE YOU ARE OVERWRITING?
filename = 'fakeday-3.h5'

# WHAT KIND OF IONOSPHERE WOULD YOU LIKE?
#   1 = CHAPMAN PROFILE, 2 = LATITUDINALLY ELONGATED ENHANCEMENT, 
#   3 = LONGITUDINALLY ELONGATED ENHANCEMENT
ionospheretype = 1

#hmf2 LOCATION [km]
hmf2 = 350

#nmf2 IN THE ABSENCE OF A PATCH [m-3]
nmf2 = 1e11
zm0 = 250000
nem0 = 1e12 # 5e11

# HOW ENHANCED IS THE PATCH RELATIVE TO THE BACKGROUND?
factor = 2

# WOULD YOU LIKE TO USE THE SOLAR ZENITH ANGLE ALREADY IN THE DATA FILE? "yes" 
#   OR "no"
szatype = 'no'
szavalue = 60

# WHAT KIND OF ERRORS WOULD YOU LIKE? 1 = UNIFORM, 2 = A CERTAIN PERCENTAGE OF 
#   THE PLASMA DENSITY
errortype = 1

# IF YOU ANSWERED 1 TO THE PREVIOUS QUESTION, ENTER THE NUMBER YOU WOULD LIKE 
#   THE BACKGROUND TO BE
# IF YOU ANSWERED 2 TO THE PREVIOUS QUESTION, ENTER THE PERCENTAGE YOU WOULD 
#   LIKE THE BACKGROUND TO BE
errorvalue =  1e11 # 1e11, 5e10, 1E12


####################################
####################################

# THIS OPENS THE DATA FILE OF INTEREST
f  = h5.File(filename, 'r+')

# EPOCH TIME
epoch = f['Time/UnixTime']
epochData = epoch[:,:]
epochData = epochData.astype(float)    
    
# ALTITUDE
alt = f['Geomag/Altitude']
altdata = alt[:,:]
altdata = altdata.astype(float)    

# CORRECTED LATITUDE [deg]
GLAT = f['Geomag/Latitude']
GLATData = GLAT[:,:]
GLATData = GLATData.astype(float)    
    
# CORRECTED LONGITUDE [deg]
GLON = f['Geomag/Longitude']
GLONData = GLON[:,:]
GLONData = GLONData.astype(float)    

# SOLAR ZENITH ANGLE [deg]
SolarZen = f['MSIS/SolarZen']
SolarZenData = SolarZen[:,:,:]
SolarZenData = SolarZenData.astype(float)   

if szatype == 'no':
  SolarZenData[:,:,:] = szavalue

# PLASMA DENSITY [m-3]
ne = f['FittedParams/Ne']
nedata = ne[:,:,:]
nedata = nedata.astype(float)    


# PLASMA DENSITY [m-3]
dne = f['FittedParams/dNe']
dnedata = dne[:,:,:]
dnedata = dnedata.astype(float)    


# fitcode
fitcode = f['FittedParams/FitInfo/fitcode']
fitcodedata = fitcode[:,:,:]
fitcodedata = fitcodedata.astype(float)    

# chi2
chi2 = f['FittedParams/FitInfo/chi2']
chi2data = chi2[:,:,:]
chi2data = chi2data.astype(float)    




kb = 1.38*10**(-23)
tr = 3000
g = 9.81
mi = 16*(1.66*10**(-27))
scaleheight = kb*tr/(mi*g)


for t in range(len(nedata)):
  for i in range(len(nedata[0])):
    for j in range(len(nedata[0,0])):
      zprime = ((altdata[i,j]-zm0)/scaleheight)
      y0 =nem0*math.exp(0.5*(1-zprime-abs(1/math.cos(SolarZenData[t,i,j]*math.pi/180))*math.exp(-1*zprime)))

      if ionospheretype == 1:
        nedata[t,i,j] = y0
        
      if ionospheretype == 2:
        ym = y0*factor
        xm = np.nanmin(GLONData)+abs((np.nanmax(GLONData)-np.nanmin(GLONData))/2)
        v0y = math.sqrt((y0 - ym)*2*(-g))
        x0 = np.nanmin(GLONData)+0.5*abs((np.nanmax(GLONData)-np.nanmin(GLONData))/2)
        v0x =(x0 - xm)*(-g)/v0y
        nedata[t,i,j] = -g*((((GLONData[i,j]-x0)/v0x)**2)/2) + v0y*((GLONData[i,j]-x0)/v0x) + y0

      if ionospheretype == 3:
        ym = y0*factor
        xm = np.nanmin(GLATData)+abs((np.nanmax(GLATData)-np.nanmin(GLATData))/2)
        v0y = math.sqrt((y0 - ym)*2*(-g))
        x0 = np.nanmin(GLATData)+0.5*abs((np.nanmax(GLATData)-np.nanmin(GLATData))/2)
        v0x =(x0 - xm)*(-g)/v0y
        nedata[t,i,j] = -g*((((GLATData[i,j]-x0)/v0x)**2)/2) + v0y*((GLATData[i,j]-x0)/v0x) + y0

      if nedata[t,i,j] < y0:
        nedata[t,i,j] = y0

      if np.isnan(altdata[i,j]) == True:
        nedata[t,i,j] = math.nan

ne[...] = nedata

dnedata = np.empty((nedata.shape[0],nedata.shape[1],nedata.shape[2]))#time, beams, ranges
for t in range(len(nedata)):
  for i in range(len(nedata[0])):
    for j in range(len(nedata[0,0])):
      if errortype == 1:
        dnedata[t,i,j] = errorvalue
      if errortype == 2:
        dnedata[t,i,j] = errorvalue*nedata[t,i,j]
      if np.isnan(altdata[i,j]) == True:
        dnedata[t,i,j] = math.nan
dne[...] = dnedata



dnedata = np.empty((nedata.shape[0],nedata.shape[1],nedata.shape[2]))#time, beams, ranges
for t in range(len(nedata)):
  for i in range(len(nedata[0])):
    for j in range(len(nedata[0,0])):
      chi2data[t,i,j] = 1
      fitcodedata[t,i,j] = 1
      if np.isnan(altdata[i,j]) == True:
        chi2data[t,i,j] = 0
        fitcodedata[t,i,j] = 0
chi2[...] = chi2data
fitcode[...] = fitcodedata


f.close() 