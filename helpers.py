import numpy as np
import netCDF4
import datetime as dt
from time import strftime
import glob
import os
from ftplib import FTP

def eraI_timeInd2Time(timeInd):
  #return the datetime of a given timeInd=hours since t0
  t0 = dt.datetime(1900,1,1,0)
  deltaT = dt.timedelta(hours=1)
  timeVal = t0+deltaT*int(timeInd)
  return timeVal

def eraI_time2TimeInd(fileTimes,t0):
  #return index of fileTimes closest to t0
  tRef = dt.datetime(1900,1,1,0)
  timeVal = (t0-tRef).total_seconds()/3600.; timeVal = int(timeVal)
  return np.argmin(np.absolute(timeVal-fileTimes))

def datetime_to_string(t0):
  #return yr-m-d
  tStruct = dt.datetime.timetuple(t0)
  ymd0 = strftime('%Y-%m-%d',tStruct); #print t0,ymd0
  return ymd0

def get_atmo_filename(event_date,days=5):
  #given datetime of case, return atmo filename.
  #example file is /arctic3/datasets/ERAI_mslp/2DayChange/ERAI_sfcVars_2007-06-26To2007-07-01.nc
  fDir = '/raid3/datasets/ERAI_mslp/2DayChange/'
  if (days>5):
    fDir = fDir+'forLags/'
  tStart=event_date-dt.timedelta(days=days)
  s0 = datetime_to_string(tStart); s1 = datetime_to_string(event_date)
  dateInfo = s0+'To'+s1
  fPath = fDir+'ERAI_sfcVars_'+dateInfo+'.nc'
  #fPath = fDir+'ERAI_DT_'+dateInfo+'.nc'
  print 'Atmo file: ', event_date, fPath
  
  return fPath

def get_atmo_filename_year(event_date,info='mslp'):
  #given datetime of case, return atmo filename.
  #use info to switch between mslp, sfcVars, DT,...
  fDir = '/raid3/datasets/ERAI_mslp/anl/'
  yr = event_date.year
  dateInfo = '{0}-01-01-to-{0}-12-31'.format(yr)
  f = 'ERAI_{0}_{1}.nc'.format(info, dateInfo)
  fPath = fDir+f
  if (True):
    print 'Atmo file: ', event_date, fPath
  return fPath

def get_ice_filenames(event_date):
  #example file is /data01/seaIce/concentration/seaice_conc_daily_nh_20140701_v02r00.nc
  #fields are daily after 1987 and every 2 days before
  #fDir = '/data01/seaIce/concentration/'
  fDir = '/raid3/datasets/seaIce/concentration/'
  deltaT = dt.timedelta(days=5)
  if (event_date.year<1988):
    deltaT = dt.timedelta(days=6)
  
  s0 = datetime_to_string(event_date-deltaT); s1 = datetime_to_string(event_date)
  s0 = s0.replace('-',''); s1 = s1.replace('-','')
  
  fStart = fDir+'seaice_conc_daily_nh_'+s0+'_v02r00.nc'
  fEnd = fDir+'seaice_conc_daily_nh_'+s1+'_v02r00.nc'
  print 'Ice files: ', event_date, fStart,fEnd
  
  return (fStart,fEnd)

def get_iceConc(f,iceKey='goddard_merged_seaice_conc',canSwitch=True):
  #return lat,lon,xice for a given case
  #the tricky bit is that not all files have valid Goddard SIC.
  
  data = netCDF4.Dataset(f,'r')
  xice = data.variables[iceKey][0,:,:]
  if (canSwitch):
    try:
      xice[xice.mask] = 0
    except:
      pass
    if (np.amax(xice)<1.e-6):
      #switch datasets
      print 'Switching ice for ',f
      iceKey = 'seaice_conc_cdr'
      xice = data.variables[iceKey][0,:,:]
  
  lat = data.variables['latitude'][:,:]
  lon = data.variables['longitude'][:,:] #[-180,180]
  
  data.close()
  return (lat,lon, xice,iceKey)
  
def get_iceConcChange(f,f0):
  #return lat,lon,change in xice for a given case
  
  #seaice data is sometimes missing between datasets. default key will be switched and when that doesn't work, hope the other is there.
  
  data = netCDF4.Dataset(f,'r')
  data0 = netCDF4.Dataset(f0,'r')
  iceKey = 'goddard_merged_seaice_conc'
  xice = data.variables[iceKey][0,:,:]
  xice0 = data0.variables[iceKey][0,:,:]
  
  try:
    xice[xice.mask] = 0
  except:
    pass
  try:
    xice0[xice0.mask] = 0
  except:
    pass

  if (np.amax(xice)<1.e-6 or np.amax(xice0)<1.e-6):
    #switch datasets
    print 'Switching ice for ',f,f0
    iceKey = 'seaice_conc_cdr'
    xice = data.variables[iceKey][0,:,:]
    xice0 = data0.variables[iceKey][0,:,:]
    
  dxice = xice-xice0 #negative if lost ice
  
  lat = data.variables['latitude'][:,:]
  lon = data.variables['longitude'][:,:] #[-180,180]
  
  data.close(); data0.close()
  return (lat,lon, dxice)

def get_monthlyMean(key='msl',scaleFac=1.e-2,iMonth=6):
  #using monthly file from ERA-I, calc mean over years and return values scaled to proper units.
  #iMonth=7 is August since 0-based file indexing.
  if (key=='msl'):
    f = '/raid3/datasets/ERAI_mslp/mslp_monthMean.nc'
  elif (key=='pt'):
    f = '/raid3/datasets/ERAI_mslp/DT_monthMean.nc'
  elif (key=='t2m' or key=='d2m'):
    f = '/raid3/datasets/ERAI_mslp/t2m_monthMean.nc'
  elif (key=='vort_500'):
    f = '/raid3/datasets/ERAI_mslp/pLev_monthMean.nc'
  else:
    print 'Unknown key: ', key
    
  data = netCDF4.Dataset(f,'r')
  vals = data.variables[key][:] #time,lat,lon
  
  valsMonth = np.mean(vals[iMonth::12,:,:],axis=0)*scaleFac
  
  return valsMonth

def calc_distSphere_multiple(r, lat1, lon1, lat2, lon2):
  '''
  #return the distance between 1 ll1 point and >=1 ll2 points.
  on a sphere.
  input lat/lon in radians!!
  '''
  dlat = lat2-lat1
  dlon = lon2-lon1
  latTerm = np.sin(.5*dlat); latTerm = latTerm*latTerm;
  lonTerm = np.sin(.5*dlon); lonTerm = lonTerm*lonTerm*np.cos(lat1)*np.cos(lat2);
  dAngle = np.sqrt(latTerm+lonTerm)
  dist = 2.*r*np.arcsin(dAngle)
  
  return dist

def formURL_iceConcentration(t0):
  #NSIDC CDR
  #example URL: ftp://sidads.colorado.edu/pub/DATASETS/NOAA/G02202_v2/north/daily/2000/seaice_conc_daily_nh_f13_20000115_v02r00.nc
  #Return -1 if date does not exist
  
  yr = t0.strftime('%Y')
  ymd = t0.strftime('%Y%m%d')
  
  ftp = FTP('sidads.colorado.edu')
  ftp.login()
  #ftp.cwd('pub/DATASETS/NOAA/G02202_v2/north/daily/{0}'.format(yr))
  ftp.cwd('pub/DATASETS/NOAA/G02202_V3/north/daily/{0}'.format(yr))
  fNames = ftp.nlst()
  
  matchName = [a for a in fNames if ymd in a]
  if (len(matchName)<1):
    f = -1
  else:
    f = 'ftp://sidads.colorado.edu/pub/DATASETS/NOAA/G02202_V3/north/daily/{0}/{1}'.format(yr, matchName[0])
  return f
  
def download_iceConcentration_forChange(t0,fDirSave = '/raid3/datasets/seaIce/concentration/'):
  #for algorithm, need a spatial change in ice concentration. we settled on ~5 days between 2 ice states. but, fields are daily only after 1987.
  #fDirSave = '/raid3/datasets/seaIce/concentration/'
  
  f = formURL_iceConcentration(t0)
  if (f==-1): #no matching time was found...assume before 1988, which was once every 2 days
    t0 = t0+dt.timedelta(days=1)
  caseTimes = []
  if (t0.year<1988):
    deltaT = dt.timedelta(days=6)
    caseTimes = [t0, t0-deltaT]
  else:
    deltaT = dt.timedelta(days=5)
    caseTimes = [t0, t0-deltaT]
    
  cmdFmt = 'wget --directory-prefix={0} {1}'
  for caseTime in caseTimes:
    f = formURL_iceConcentration(caseTime)
    cmd = cmdFmt.format(fDirSave,f); print cmd; os.system(cmd)

def softLink_iceFileNames():
  #there's weird f08, n07, f13,f11,f17,... descriptors on the filenames from nsidc.
  #make symbolic links to unified filenames
  
  specList = ['n07_','f08_','f11_','f13_','f17_']; versionReplace = ['v03r01','v02r00']
  for spec in specList:
    fList = glob.glob('*{0}*.nc'.format(spec))
    for f in fList:
      svName = f.replace(spec,''); svName = svName.replace(versionReplace[0], versionReplace[1])
      cmd = 'ln -s {0} {1}'.format(f, svName)
      print cmd; os.system(cmd)
      
if __name__=='__main__':
  softLink_iceFileNames()


