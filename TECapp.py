import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import datetime
import aacgmv2
import magtec
import pyIGRF
import os
from PIL import Image
from numpy import genfromtxt
import spacepy.coordinates as coord
from spacepy.time import Ticktock

def main():
      

    function = st.sidebar.selectbox('Plot Type', ('Interactive Plotter', '3 Hour Global Plots', 'Daily Movie'))
    date = st.sidebar.date_input('Date', datetime.date(2017, 9, 1))
    st.sidebar.markdown('Only dates in 2017 are available, more data coming soon')
    
    if function == 'Interactive Plotter':
        st.title('Total Electron Content Interactive Plotting Tool')
        coordsystem = st.sidebar.selectbox('Coordinate System', ('Geographic', 'AACGMv2', 'IGRF'))
        supermag = st.sidebar.checkbox('Overlay magnetometers?')
        maxtec = st.sidebar.slider('Max TEC Value to Plot, (TECU)', 0, 50, 20)
        st.sidebar.markdown('TEC units are vertical column density measurments.')
        st.sidebar.markdown('1 TECU = 10^16 electrons/m^2')
        intfactor = st.sidebar.slider('Interpolation Level', 0, 26, 4)
        st.sidebar.markdown('Interpolation level is the minimum number of datapoints in the 3x3x3 matrix of latitude, longitude, and time that are needed to "fill in" the point in question')
        medfilter = st.sidebar.checkbox('Apply median filtering?')
        st.sidebar.markdown('Median filtering sets the TEC value of all datapoints to the median of the 3x3x3 matrix of latitude, longitude, and time surrounding it')
        time = st.sidebar.time_input('Time (UT)', datetime.time(12, 00))
        plottype = st.sidebar.selectbox('Plot Layout', ('Northern Hemisphere', 'Southern Hemisphere', 'Global'))
        extent = [0, 0, 0, 0]
        extent[0]= st.sidebar.number_input('Minimum Longitude', min_value=-180, max_value=180, value=-180)
        extent[1]= st.sidebar.number_input('Maximum Longitude', min_value=-180, max_value=180, value=180)
        extent[2]= st.sidebar.number_input('Minimum Latitude', min_value=-90, max_value=90, value=0)
        extent[3]= st.sidebar.number_input('Maximum Latitude', min_value=-90, max_value=90, value=90)
        lat, lon, Z = convert(date, time, intfactor, coordsystem, medfilter)
        plot(lat, lon, Z, date, time, intfactor, maxtec, coordsystem, plottype, extent, supermag)
    elif function == '3 Hour Global Plots':
        st.title('Three Hour Global Total Electron Content Plots')
        showdailyplots(date)
    elif function == 'Daily Movie':  
        st.title('Global Total Electron Content Animation')
        showmovie(date)
        st.write('Animation resolution has been reduced for display purposes.  Higher resolution versions are availabe upon request.')
    st.write('All TEC data is from http://cedar.openmadrigal.org/,')
    st.write('Total electron content is a measurment of the vertical column density of electrons in the atmosphere.  It is usually found using phase delays from GPS signals.  TEC data used in this webtool is binned into single degree longitude and latitude estimates and available in 15 minute intervals.')
def showdailyplots(date):
    os.chdir(str(date.year))
    os.chdir(str(date.month))
    os.chdir(str(date.day))
    st.image('merged_images.png')
    os.chdir('..')
    os.chdir('..')
    os.chdir('..')

def showmovie(date):
    os.chdir(str(date.year))
    os.chdir(str(date.month))
    os.chdir(str(date.day))
    st.image('movie.gif')
    os.chdir('..')
    os.chdir('..')
    os.chdir('..')
    
def timeconvert(time):
    hours = time.hour * 60
    totalmins = hours + time.minute
    totalmins = totalmins / 5
    if(totalmins == 0):
        totalmins = 5
    return int(totalmins)

def interpolate(X, tec, time, factor):
    for index in np.ndindex(X.shape):
        if np.isnan(X[index]):
            Y = tec[(index[0] - 1):(index[0] + 2), (index[1] - 1) : (index[1] + 2), time - 1: time + 2]
            if np.count_nonzero(~np.isnan(Y)) > factor:
                X[index] = np.nanmedian(Y)
    return X     

def medianfilter(X, tec, time):
    for index in np.ndindex(X.shape):
        if  not (np.isnan(X[index])):
            Y = tec[(index[0] - 1):(index[0] + 2), (index[1] - 1) : (index[1] + 2), time - 1: time + 2]
            X[index] = np.nanmedian(Y)
    return X 

def geotomag(lat,lon,alt,plot_date):
    #call with altitude in kilometers and lat/lon in degrees 
    Re=6371.0 #mean Earth radius in kilometers
    #setup the geographic coordinate object with altitude in earth radii 
    cvals = coord.Coords([np.float(alt+Re)/Re, np.float(lat), np.float(lon)], 'GEO', 'sph',['Re','deg','deg'])
    #set time epoch for coordinates:
    cvals.ticks=Ticktock([plot_date.isoformat()], 'ISO')
    #return the magnetic coords in the same units as the geographic:
    return cvals.convert('MAG','sph')

def convert(date, time, intfactor, coordsystem, medfilter):
    directory = str(os.getcwd())
    directory = directory + '\\' + str(date.year) + '\\' + str(date.month) + '\\' + str(date.day)
    with np.load(directory + '\\data.npz') as data:
        tec = data['tec']        
        timeint = timeconvert(time)        
        Z = tec[:,:,timeint]
        
        Z = interpolate(Z, tec, timeint, intfactor) 
        if medfilter:
            Z = medianfilter(Z, tec, timeint)
                       
        lon = np.linspace(-180, 180, 361)
        lat = np.linspace(-90, 90, 181)
        lon2d = np.empty((180, 360))
        lat2d = np.empty((180, 360))
        lon2d[:] = np.nan
        lat2d[:] = np.nan
        if(coordsystem == 'AACGMv2'):
            for i in range(0,179):
                for j in range(0, 359):
                    lat2d[i,j], lon2d[i,j], x = aacgmv2.get_aacgm_coord(i - 90, j - 180, 350, date)
                    
            #remove jumps in longitude by rearranging matrix of lon, lat, and Z
            for row in range(0, 179):
                for element in range(0,358):
                    if (lon2d[row, element] > 0) and (lon2d[row, element+1] < 0) and (180 < (lon2d[row, element] - lon2d[row, element+1])):
                        first = lon2d[row,:element+1]
                        second = lon2d[row,element+1:]
                        lon2d[row,:] = np.append(second, first)
                        first = lat2d[row,:element+1]
                        second = lat2d[row,element+1:]
                        lat2d[row,:] = np.append(second, first)
                        first = Z[row,:element+1]
                        second = Z[row,element+1:]
                        Z[row,:] = np.append(second, first)
            return lat2d, lon2d, Z
        elif(coordsystem == 'IGRF'):
            lon = np.linspace(-180, 180, 361)
            lat = np.linspace(-90, 90, 181)
            lon2d = np.empty((180, 360))
            lat2d = np.empty((180, 360))
            lon2d[:] = np.nan
            lat2d[:] = np.nan
            for i in range(0,179):
                for j in range(0, 359):
                    _mc = geotomag(i - 90, j - 180, 350, date)
                    lat2d[i,j] = float(_mc.data[:,1])
                    lon2d[i,j] = float(_mc.data[:,2])
            return lat2d, lon2d, Z
        else:
            return lat, lon, Z#Note that return could be 1d or 2d array



    

def plot(lat, lon, Z, date, time, intfactor, maxtec, coordsystem, plottype, extent, supermag):

        
    tecmap = plt.figure(figsize=(11,8))
    magstations = genfromtxt('20200515-23-45-supermag-stations.csv', delimiter=',')
    if plottype ==  'Northern Hemisphere':
        map_proj = ccrs.NorthPolarStereo()
    elif plottype == 'Southern Hemisphere':
        map_proj = ccrs.SouthPolarStereo()
    elif plottype == 'Global':
        map_proj = ccrs.PlateCarree()
    if coordsystem=='AACGMv2' or coordsystem=='IGRF':
        ax = tecmap.add_subplot(projection='aacgmv2',map_projection = map_proj)
        ax.overaly_coast_lakes(coords="aacgmv2", plot_date=date)
                
        
        if map_proj == ccrs.PlateCarree():
            mesh = ax.pcolor(lon, lat, Z, cmap='jet', vmax=maxtec, transform=ccrs.PlateCarree()) 
        else:
            mesh = ax.scatter(lon, lat, c=Z, cmap='jet', vmax=maxtec, transform=ccrs.PlateCarree())
    
    else:
        #usepcolormesh in 
        ax = plt.axes(projection=map_proj)
        ax.add_feature(cfeature.COASTLINE)
        #ax.add_feature(cfeature.LAKES)
        mesh = ax.pcolormesh(lon, lat, Z, cmap='jet', vmax=maxtec, transform=ccrs.PlateCarree())
        if(supermag):
            magx = magstations[1:, 1]
            magx = magx
            magy = magstations[1:, 2]
            ax.scatter(magx, magy, c='black', transform=ccrs.PlateCarree())            
    
    ax.set_extent(extent, ccrs.PlateCarree())
    clrbar = plt.colorbar(mesh, shrink=0.5)
    clrbar.set_label('Total Electron Content (TECU)')
    ax.gridlines(linewidth=0.5)
    plt.title('Total Electron Content for ' + str(date.month) + '/' + str(date.day) + '/' + str(date.year) + ' at ' + str(time.hour) + ':' + ('%02d' % time.minute) + ' UT')
    st.pyplot(tecmap)


main()
    

    
