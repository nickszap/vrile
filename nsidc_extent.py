import numpy as np
from scipy import signal, stats
import matplotlib as mlab
import matplotlib.pyplot as plt
import datetime as dt

import helpers

def read_ExtentCSV(f):
  data = np.loadtxt(f,skiprows=2, delimiter=',',usecols=(0,1,2,3))
  yr = data[:,0].astype(int) #e.g., 1981
  month = data[:,1].astype(int) #e.g., 01
  day=data[:,2].astype(int) #e.g., 01
  extent = data[:,3] #e.g., 11.345 [10^6km^2]
  return (yr,month,day,extent)
  
#Butterworth filtering from scipy-cookbook: http://scipy-cookbook.readthedocs.io/items/ButterworthBandpass.html
def butter_highpass(fCut, fs, order=9):
  nyq = 0.5 * fs
  low = fCut / nyq
  b, a = signal.butter(order, low, btype='highpass')
  return b, a

def butter_highpass_filter(data, fCut, fs, order=9):
  #answers are different if use lfilter vs filtfilt...i'll just be consistent w/ what uriel did
  b, a = butter_highpass(fCut, fs, order=order)
  #y = signal.lfilter(b, a, data)
  y = signal.filtfilt(b,a,data)
  return y

#def driver_getCaseDates(perc=1,showFigs=False,monthMin=6,monthMax=8,dtUniqueCase=dt.timedelta(days=5)):
def driver_getCaseDates(perc=1,showFigs=False,validMonth=range(6,8+1),dtUniqueCase=dt.timedelta(days=5)):
  #read sie time series ------------------------
  #f = '/data01/class/reu/NH_seaice_extent_final.csv'
  #f = '/data01/class/reu/N_seaice_extent_daily_v2.1.1979-2016.csv'
  f = '/raid3/datasets/seaIce/N_seaice_extent_daily_v2.1.1979-2016.csv'
  yr,month,day,extent = read_ExtentCSV(f)
  nVals = len(extent)
  timeList = [dt.datetime(yr[i],month[i],day[i]) for i in xrange(nVals)]; timeList = np.array(timeList)
  
  #filter sie ------------------------
  #Interpolate to uniform grid
  tStart = timeList[0]; tEnd = timeList[-1]
  nDays = (tEnd-tStart).total_seconds()/(24.*3600.); nDays = int(nDays); print '{0} days between {1} {2}'.format(nDays, tStart, tEnd)
  times = [tStart+ dt.timedelta(days=i) for i in range(0, nDays)]; times = np.array(times)
  
  dTime = [(a-timeList[0]).total_seconds() for a in timeList]
  dTimeResample = [(a-timeList[0]).total_seconds() for a in times]
  extent = np.interp(dTimeResample, dTime, extent)
  
  #impose filter - fourier bandpass, butterworth,...
  #filter setup
  fs = 1./(24.*3600.) #1/day in samples/s
  fCut = 1./(18.*24.*3600.)
  #windowing for finite length time series really changes time series. seems better to just ignore values near boundaries ("boundary distance" depends on order of filter)
  extent -= np.mean(extent)
  wts = np.kaiser(extent.size,14)
  extentF = extent*wts
  if (showFigs):
    plt.figure()
    plt.plot(times,extent,label='interp')
    plt.plot(times,extentF,label='window')
    plt.legend()
    plt.show()
  
  #extent_highPass = butter_highpass_filter(extent, fCut, fs, order=9)
  
  dExtent = extent.copy()
  dExtent[2:] = extent[2:]-extent[0:-2]
  dExtent[0:2] = dExtent[2] #constant extrapolation
  dExtent_highPass = butter_highpass_filter(dExtent, fCut, fs, order=9)  
  
  #toss values near endpoints
  #timeCutoff = dt.datetime(1980,1,1); isValid = times>timeCutoff; timeCutoff = dt.datetime(2014,1,1); isValid = isValid*(times<timeCutoff)
  #times = times[isValid]; extent = extent[isValid]; dExtent_highPass = dExtent_highPass[isValid]; #extent_highPass = extent_highPass[isValid]
  
  if (showFigs):
    fig, ax1 = plt.subplots()
    ax1.plot(times,extent,'k-',label='extent')
    ax2 = ax1.twinx()
    ax2.plot(times,dExtent_highPass,'b-',label='Filtered dExtent')
    #make common legend: http://stackoverflow.com/questions/5484922/secondary-axis-with-twinx-how-to-add-to-legend
    # ask matplotlib for the plotted objects and their labels
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2)
    plt.show()
  
  #candidate cases ------------------------
  vals = dExtent_highPass
  #vals = extent_highPass
  #vals = dExtent
  thresh = np.percentile(vals, perc)
  print '{0}% threshold for dExtent: {1}'.format(perc,thresh)
  
  caseDates = times[vals<=thresh]
  caseVals = vals[vals<=thresh]
  print caseDates.tolist()
  
  if True:
    iDate = times.tolist().index(dt.datetime(2006,8,20))
    valDate = vals[iDate]
    print times[iDate], valDate, np.sum(vals<valDate), np.sum(vals<valDate)/float(len(vals))
  
  if (showFigs):
    plt.figure()
    plt.plot(caseDates, caseVals,'bo')
    plt.show()
  
    months = [i.month for i in caseDates]
    plt.hist(months,bins=12)
    plt.show()
  
  #filter candidate cases - season, too close in time,... ------------------------
  #get only w/in specified month range
  print 'Limiting to months in: ', validMonth
  #caseDates = [t0 for t0 in caseDates if ((t0.month>=monthMin) and (t0.month<=monthMax)) ]
  caseDates = [t0 for t0 in caseDates if (t0.month in validMonth) ]
  
  #eliminate dates close in time
  print 'Keeping later of dates for duplicate cases if w/in ', dtUniqueCase
  nCases = len(caseDates)
  validDates = []
  for iCase in xrange(nCases-1):
    if (caseDates[iCase+1]-caseDates[iCase] > dtUniqueCase):
      validDates.append(caseDates[iCase])
  validDates.append(caseDates[-1])
  print validDates
  
  if (True):
    #take out dates w/ missing values
    caseDates = []
    for t0 in validDates:
      #19870709 is first file in 1987
      if (t0<dt.datetime(1987,1,1) or t0>dt.datetime(1987,7,10) ):
        caseDates.append(t0)
    print 'Kept {0}/{1} event dates not known to have missing sea ice values'.format(len(caseDates),len(validDates))
    validDates = caseDates
  
  return validDates

def get_seaIceCases():
  #validDates = driver_getCaseDates(perc=5)
  validDates = driver_getCaseDates(perc=5,showFigs=True,validMonth=range(0,13),dtUniqueCase=dt.timedelta(days=5))
  for t0 in validDates:
    helpers.download_iceConcentration_forChange(t0)
  
if __name__=='__main__':
  driver_getCaseDates(perc=1)


