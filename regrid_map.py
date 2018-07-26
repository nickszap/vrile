import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata
from netCDF4 import Dataset

lon0Map = -97.44; lat0Map = 35.22 #Norman OK

def make_referenceGrid(lat0 = lat0Map, lon0 = lon0Map,dxRef=30.e3, refGridRadius=3000.e3):
  #make a computational grid that's evenly spaced centered on our map location
  #specs for map and grid
  wdth = 10.e6;
  
  #make map
  m0 = Basemap(width=wdth,height=wdth,resolution='l',projection='stere',lat_ts=lat0,lat_0=lat0,lon_0=lon0)
  
  #make grid
  x0,y0 = m0(lon0,lat0)
  nGrid = int(refGridRadius/dxRef)
  refX = x0+dxRef*np.arange(-nGrid, nGrid); refY = y0+dxRef*np.arange(-nGrid, nGrid); nx = len(refX); ny = len(refY)
  gridX, gridY = np.meshgrid(refX, refY); #print gridX; print gridY
  
  return (gridX,gridY, m0)
  
def driver(latRefCase, lonRefCase, lat2d, lon2d, vals2d):
  #for a given case, map case values onto reference grid.
  
  #make reference grid common to all cases
  gridX,gridY,mRef = make_referenceGrid()
  
  #make a case grid
  gridX_case,gridY_case,mCase = make_referenceGrid(lat0 = latRefCase, lon0 = lonRefCase)
  
  #interpolate values of case to reference grid
  xCase, yCase = mCase(lon2d,lat2d)
  vals_dest = interpolate_nearest(xCase,yCase,gridX,gridY,vals2d)
  
  return (gridX,gridY, mRef,vals_dest)

def interpolate_nearest(xsource,ysource,xdest,ydest,values_source):
  vals = values_source.flatten()
  
  pts_src = np.column_stack( (xsource.flatten(), ysource.flatten()) )
  pts_dest = np.column_stack( (xdest.flatten(), ydest.flatten()) )
  vals_dest = griddata(pts_src, vals, pts_dest, method='nearest')
  
  vals_dest = np.reshape(vals_dest, xdest.shape)
  return (vals_dest)
  
if __name__=='__main__':
  pass


