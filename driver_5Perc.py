import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import datetime as dt

import nsidc_extent
import seaIceMasks
import helpers

timeStringFormat = "%Y-%m-%d-%H"

def calc_eventLocations_concLoc(perc=1,showFigs=False,validMonth=range(6,8+1),dtUniqueCase=dt.timedelta(days=5), yrMax=2015):
  #Within the largest loss object, take the location w/ largest concentration loss
  
  caseDates = nsidc_extent.driver_getCaseDates(perc=perc,showFigs=showFigs,validMonth=validMonth,dtUniqueCase=dtUniqueCase)
  print 'Keeping years .le. ', yrMax; caseDates = [t0 for t0 in caseDates if t0.year<=yrMax]
  
  infoCases = []
  for t0 in caseDates:
    #ice loss mask
    fIceStart,fIceEnd = helpers.get_ice_filenames(t0);
    latIce, lonIce, dConcExtr = seaIceMasks.calc_iceLossMask_refPoint(fIceStart, fIceEnd,returnMask=False,dxiceThresh=-.1)
    
    print 'Event location: ', t0, latIce, lonIce
    infoCases.append([t0, latIce, lonIce])
  
  if (showFigs):
    iceVals = 0
    for t0 in caseDates:
      fIceStart,fIceEnd = helpers.get_ice_filenames(t0);
      lat,lon, mask = seaIceMasks.calc_iceLossMask_refPoint(fIceStart, fIceEnd,returnMask=True, dxiceThresh=-.1)
      
      iceVals += mask
    seaIceMasks.plot_2d_ll(lat,lon, iceVals, cmap=plt.cm.RdBu_r, showFig=True)
  
  print infoCases
  return infoCases

def calc_loc_maxWind_iceLossObject(u,v,mask):
  #Input u,v[times,lats,lons]. mask[lats,lons].
  #return iTime,iLat,iLon, and wind angle (math angle wrt x-axis) for max wind over valid mask
  
  speed = np.sqrt(u*u+v*v)
  speed = speed*mask[None,:,:]
  ind = np.argmax(speed)
  iTime,iLat,iLon = np.unravel_index(ind, speed.shape)
  
  angle = np.arctan2(v[iTime,iLat,iLon],u[iTime,iLat,iLon])
  
  return (iTime,iLat,iLon,angle)

def calc_eventLocations_windMax(perc=1,showFigs=False,validMonth=range(6,8+1),dtUniqueCase=dt.timedelta(days=5), dtWindMax=dt.timedelta(days=5),yrMax=2015):
  #Within the largest loss object interpolated to atmospheric mesh, take the location of maximum wind over the previous dtWindMax.
  #The resulting date and location should correspond more closely to driving atmospheric conditions.
  
  caseDates = nsidc_extent.driver_getCaseDates(perc=perc,showFigs=showFigs,validMonth=validMonth,dtUniqueCase=dtUniqueCase)
  print 'Keeping years .le. ', yrMax; caseDates = [t0 for t0 in caseDates if t0.year<=yrMax]
  
  #get atmosphere grid
  fAtmo = helpers.get_atmo_filename_year(caseDates[0],info='sfc')
  data = netCDF4.Dataset(fAtmo,'r')
  lon1d = data.variables['longitude'][:]
  lat1d = data.variables['latitude'][:]
  data.close()
  
  infoCases = []
  for t0 in caseDates:
    #ice loss mask
    fIceStart,fIceEnd = helpers.get_ice_filenames(t0);
    atmoLat2d, atmoLon2d, iceMask_atmo = seaIceMasks.make_ice_mask(fIceStart,fIceEnd,lat1d,lon1d)
    
    fAtmo = helpers.get_atmo_filename_year(t0, info='sfc')
    data = netCDF4.Dataset(fAtmo,'r')
    fileTimes = data.variables['time'][:]
    iTime1 = helpers.eraI_time2TimeInd(fileTimes,t0); iTime0 = helpers.eraI_time2TimeInd(fileTimes,t0-dtWindMax)
    u = data.variables['u10'][iTime0:iTime1+1,:,:]
    v = data.variables['v10'][iTime0:iTime1+1,:,:]
    data.close()
    
    iTime,iLat,iLon,angle = calc_loc_maxWind_iceLossObject(u,v,iceMask_atmo)
    tWind = helpers.eraI_timeInd2Time(fileTimes[iTime0+iTime])
    print 'Event location: ', t0, tWind,lat1d[iLat],lon1d[iLon],angle*180./np.pi
    
    infoCases.append((t0,tWind,lat1d[iLat],lon1d[iLon],angle))
  
  return infoCases

if __name__=='__main__':
  #infoCases = calc_eventLocations_concLoc(perc=5); print infoCases
  infoCases = calc_eventLocations_windMax(perc=5); print infoCases
  


