import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import scipy.ndimage
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata

import helpers
import regrid_map

def plot_2d_ll(lat, lon, var, cmap=plt.cm.RdBu_r, showFig=True, vMinMax=None):
  #Input lat and lon in degrees!!!
  plt.figure()
  
  m = Basemap(projection='npstere',boundinglat=50,lon_0=0,resolution='l')
  #m = Basemap(projection='ortho',lon_0=0,lat_0=90, resolution='l')
  x,y = m(lon, lat)
  
  m.drawcoastlines()
  m.drawmapboundary()
  
  if (vMinMax != None):
    m.pcolor(x,y,var,shading='flat',edgecolors='none',cmap=cmap,vmin=vMinMax[0],vmax=vMinMax[1])
  else:
    m.pcolor(x,y,var,shading='flat',edgecolors='none',cmap=cmap) #cmap=plt.cm.hot_r) #,vmin=100,vmax=1000)
  
  plt.colorbar()
  if (showFig==True):
    plt.show()
  else:
    print 'Saving to: ', showFig
    plt.savefig(showFig); plt.close()

def calc_objects(candidates, szFilter=None, takeBiggest=True):
  candidates = scipy.ndimage.morphology.binary_closing(candidates, structure=np.ones((3,3),dtype=int), iterations=1)
  objs, nObjs = scipy.ndimage.measurements.label(candidates)
  
  labelForTooSmall = -1
  if (szFilter != None):
    #eliminate objects that are too small
    #labelForTooSmall = -1
    for objLabel in xrange(nObjs+1):
      inObj = objs==objLabel
      nInObj = np.sum(inObj)
      if (nInObj<=szFilter):
        objs[inObj] = labelForTooSmall
  if (takeBiggest):
    #set every region that's not the largest to tooSmall
    objSizes = []; objLabels = []
    for objLabel in xrange(1,nObjs+1):
      inObj = objs==objLabel
      nInObj = np.sum(inObj)
      objLabels.append(objLabel); objSizes.append(nInObj)
    iBig = np.argmax(objSizes); objLabel = objLabels[iBig]
    objs[objs!=objLabel] = labelForTooSmall

  if (szFilter != None or takeBiggest):      
    #reorder the object labels so can loop through    
    uniqueObjs = np.unique(objs); uniqueObjs.sort(); #so -1 for objects that were too small to be first element
    nUniqueObjs = len(uniqueObjs)
    newLabels = objs.copy() #to not overwrite
    for iObj in xrange(nUniqueObjs):
      oldLabel = uniqueObjs[iObj]
      #newLabels[objs==oldLabel] = iObj-1
      newLabels[objs==oldLabel] = iObj
    nObjs = nUniqueObjs -1 #-1 since we toss objects that were too small
    objs = newLabels
    
  return (objs, nObjs)

def regrid_nearest(lat_src, lon_src, vals, lat_dest, lon_dest):
  #given 2d lat/lons for src and dest, interpolate vals src->dest grids
  
  #distance works naturally on a map projection, so we'll use that
  m = Basemap(projection='npstere',boundinglat=50,lon_0=0,resolution='l')
  x_src, y_src = m(lon_src, lat_src)
  x_dest, y_dest = m(lon_dest, lat_dest)
  
  vals_dest = regrid_map.interpolate_nearest(x_src,y_src,x_dest,y_dest,vals)
  
  return vals_dest
  
def demo():
  #user input
  #ice='/data01/class/reu/seaice_conc_daily_nh_20050721_v02r00.nc'
  #iceRef='/data01/class/reu/seaice_conc_daily_nh_20050716_v02r00.nc'
  #iceRef = '/arctic3/datasets/seaIce/concentration/seaice_conc_daily_nh_19790713_v02r00.nc'
  #ice = '/arctic3/datasets/seaIce/concentration/seaice_conc_daily_nh_19790719_v02r00.nc'
  ice = '/raid3/datasets/seaIce/concentration/seaice_conc_daily_nh_n07_19840722_v02r00.nc'
  iceRef = '/raid3/datasets/seaIce/concentration/seaice_conc_daily_nh_n07_19840728_v02r00.nc'

  dLat = .5; dLon = dLat
  atmoLat = np.arange(90,-90-dLat/2,-dLat); atmoLon = np.arange(0,360,dLon)
  atmoLat, atmoLon = np.meshgrid(atmoLat, atmoLon)
  
  #get ice change for event
  lat,lon, dxice = helpers.get_iceConcChange(ice, iceRef)
  plot_2d_ll(lat, lon, dxice, showFig=False)
  
  #group into "significant" regions
  candidates = dxice<-0.1
  #objs, nObjs = calc_objects(candidates, szFilter=4,takeBiggest=False)
  objs, nObjs = calc_objects(candidates, szFilter=4,takeBiggest=True)
  plot_2d_ll(lat, lon, objs, cmap=plt.cm.Set1,showFig=False)
  iceLoss = objs>0

  #interpolate candidate regions onto atmo grid
  iceLoss_atmo = regrid_nearest(lat, lon, iceLoss, atmoLat, atmoLon)
  plot_2d_ll(atmoLat, atmoLon, iceLoss_atmo, showFig=False)
  
  plt.show()

def make_ice_mask(fname1,fname2,atmoLat,atmoLon):
  ice=fname2
  iceRef=fname1
  #dLat = .5; dLon = dLat
  #atmoLat = np.arange(90,-90-dLat/2,-dLat); atmoLon = np.arange(0,360,dLon)
  atmoLat, atmoLon = np.meshgrid(atmoLat, atmoLon)
  atmoLat = atmoLat.T; atmoLon = atmoLon.T; #transpose for lat,lon since meshgrid is silly

  #get ice change for event
  lat,lon, dxice = helpers.get_iceConcChange(ice, iceRef)
  
  #group into "significant" regions
  candidates = dxice<-0.1
  objs, nObjs = calc_objects(candidates, szFilter=4)
  #plot_2d_ll(lat, lon, objs, cmap=plt.cm.Set1,showFig=False)
  iceLoss = objs>0

  #interpolate candidate regions onto atmo grid
  iceLoss_atmo = regrid_nearest(lat, lon, iceLoss, atmoLat, atmoLon)

  return (atmoLat, atmoLon, iceLoss_atmo)

def calc_iceLossMask_refPoint(fname1,fname2,returnMask=False, dxiceThresh=-.1):
  ice=fname2
  iceRef=fname1

  #get ice change for event
  lat,lon, dxice = helpers.get_iceConcChange(ice, iceRef)
  
  #group into "significant" regions
  #candidates = dxice<-0.1
  candidates = dxice<dxiceThresh
  objs, nObjs = calc_objects(candidates, szFilter=4)
  #plot_2d_ll(lat, lon, objs, cmap=plt.cm.Set1,showFig=False)
  iceLoss = objs>0
  
  if (returnMask):
    return (lat,lon,iceLoss)
  
  cellInd = np.unravel_index(np.argmin(iceLoss*dxice), dxice.shape)
  latPt = lat[cellInd]; lonPt = lon[cellInd]
  
  return (latPt, lonPt, dxice[cellInd])

if __name__=='__main__':
  demo()

