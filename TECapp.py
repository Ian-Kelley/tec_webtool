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


def main():
      
    st.title('TEC Online Plotter')
    st.write('This tool plots binned TEC data from Madrigal.')
    st.write('Daily overview GIF and 3 hour plots are also available.')
    aacgm = False#st.sidebar.checkbox('Use AACGM?')
    date = st.sidebar.date_input('Date', datetime.date(2017, 9, 1))
    function = st.sidebar.selectbox('Plot Type', ('Interactive Plotter', '3 Hour Global Plots', 'Daily Movie'))
    if function == 'Interactive Plotter':
        
        maxtec = st.sidebar.slider('Max TEC Value to Plot', 0, 50, 20)
        intfactor = st.sidebar.slider('Interpolation Level', 0, 26, 4)
        time = st.sidebar.time_input('Time', datetime.time(12, 00))
        plottype = st.sidebar.selectbox('Plot Layout', ('Northern Hemisphere', 'Southern Hemisphere', 'Global'))
        extent = [0, 0, 0, 0]
        extent[0]= st.sidebar.number_input('Xmin', min_value=-180, max_value=180, value=-180)
        extent[1]= st.sidebar.number_input('Xmax', min_value=-180, max_value=180, value=180)
        extent[2]= st.sidebar.number_input('Ymin', min_value=-90, max_value=90, value=0)
        extent[3]= st.sidebar.number_input('Ymax', min_value=-90, max_value=90, value=90)
        lat, lon, Z = convert(date, time, intfactor, aacgm)
        plot(lat, lon, Z, date, time, intfactor, maxtec, aacgm, plottype, extent)
    elif function == '3 Hour Global Plots':
        st.title('Three Hour Global TEC Plots')
        showdailyplots(date)
    elif function == 'Daily Movie':  
        st.title('Global TEC Animation')
        showmovie(date)
        st.write('Animation resolution has been reduced for display purposes.  Higher resolution versions are availabe upon request.')

def showdailyplots(date):
    os.chdir(str(date.year))
    os.chdir(str(date.month))
    os.chdir(str(date.day))
    st.image('merged_images.png', width = 800)
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

def convert(date, time, intfactor, aacgm):
    directory = str(os.getcwd())
    directory = directory + '\\' + str(date.year) + '\\' + str(date.month) + '\\' + str(date.day)
    with np.load(directory + '\\data.npz') as data:
        tec = data['tec']        
        timeint = timeconvert(time)        
        Z = tec[:,:,timeint]
        Z = interpolate(Z, tec, timeint, intfactor)                
        lon = np.linspace(-180, 180, 361)
        lat = np.linspace(-90, 90, 181)
        lon2d = np.empty((180, 360))
        lat2d = np.empty((180, 360))
        lon2d[:] = np.nan
        lat2d[:] = np.nan
        if(aacgm):
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
    
        else:
            return lat, lon, Z#Note that return could be 1d or 2d array


    

def plot(lat, lon, Z, date, time, intfactor, maxtec, aacgm, plottype, extent):

        
    tecmap = plt.figure(figsize=(11,8))
        
    if plottype ==  'Northern Hemisphere':
        map_proj = ccrs.NorthPolarStereo()
    elif plottype == 'Southern Hemisphere':
        map_proj = ccrs.SouthPolarStereo()
    elif plottype == 'Global':
        map_proj = ccrs.PlateCarree()
    if aacgm:
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
        ax.add_feature(cfeature.LAKES)
        mesh = ax.pcolormesh(lon, lat, Z, cmap='jet', vmax=maxtec, transform=ccrs.PlateCarree())  
            
    
    ax.set_extent(extent, ccrs.PlateCarree())
    clrbar = plt.colorbar(mesh, shrink=0.5)
    clrbar.set_label('Total Electron Content (TECU)')
    ax.gridlines(linewidth=0.5)
    plt.title('Total Electron Content for ' + str(date.month) + '/' + str(date.day) + '/' + str(date.year) + ' at ' + str(time.hour) + ':' + ('%02d' % time.minute) + ' UT')
    st.pyplot(tecmap)


main()
    
