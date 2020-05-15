import os
import numpy as np
import datetime
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import PIL
from PIL import Image
from madrigalWeb import madrigalWeb
import h5py


def main():

    day = 1
    date = datetime.date(2017, 1, day)
    os.chdir(str(date.year))
    os.chdir(str(date.month))
    while( day < 31):
        
        date = datetime.date(2017, 1, day)
        #change working directory to date of intrest
        if(os.path.exists(str(date.day))):
            os.chdir(str(date.day))
        else:
            os.mkdir(str(date.day))
            os.chdir(str(date.day))
   
        getdata(date)
        pnggenerate(date)
    
        os.chdir('..')

        day = day + 1


def getdata(date):
    url = "http://cedar.openmadrigal.org/"
    username = "Ian Kelley"
    email = "ikelley@vt.edu"
    affiliation = "Virginia Tech"
    instrument = 8000 #TEC instruments

    startyear = str(date.year)
    startmonth = str(date.month)
    startday = str(date.day)
    starthour = '0'
    startmin = '0'
    startsec = '0'
        



    MadrigalData = madrigalWeb.MadrigalData(url)
    experimentarray = MadrigalData.getExperiments(instrument, startyear, startmonth, startday, starthour, startmin, startsec, startyear, startmonth, startday, starthour, startmin, startsec, 1)

    #save hdf5 format straight from madrigalweb (no parameters)
    print("Downloading: ", startyear, startmonth, startday)
    filearray = MadrigalData.getExperimentFiles(experimentarray[1].id, getNonDefault=False) #need to use index 1 since madrigal returns the previos day at index 0 
    tempfilename = 'tempfile.hdf5'
    MadrigalData.downloadFile(filearray[0].name, tempfilename, username, email, affiliation, format='hdf5')
    print("raw hdf5 file downloaded successfully")

    
    #save only 3d TEC matrix binned by lat, lon, and time
    with h5py.File( tempfilename, 'r') as f:
        data = f['Data']['Array Layout']['2D Parameters']
        tec = data['tec']
        dtec = data['dtec']
        print("saving to compressed .npz file format")
        np.savez_compressed('data.npz', tec=tec, dtec=dtec)


    os.remove(tempfilename)

def interpolate(X, tec, time, factor):
    if time == 0:
        time = 1
    for index in np.ndindex(X.shape):
        if np.isnan(X[index]):
            Y = tec[(index[0] - 1):(index[0] + 2), (index[1] - 1) : (index[1] + 2), time - 1: time + 2]
            if np.count_nonzero(~np.isnan(Y)) > factor:
                X[index] = np.nanmedian(Y)
    return X  

def plot(Z, time, date):
    lon = np.linspace(-180, 180, 361)
    lat = np.linspace(-90, 90, 181)
    plt.figure(figsize=(11,8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines()                
    mesh = ax.pcolormesh(lon, lat, Z, cmap='jet', vmax=25, transform=ccrs.PlateCarree())
    clrbar = plt.colorbar(mesh, shrink=0.5)
    clrbar.set_label('Total Electron Content (TECU)')
    ax.gridlines(linewidth=0.5)
    plt.title('Total Electron Content for ' + str(date.month) + '/' + str(date.day) + '/' + str(date.year) + ' at ' + str(time.hour) + ':' + ('%02d' % time.minute) + ' UT')
    plt.savefig(str(time.hour) + ('%02d' % time.minute) + '.png', format="png", bbox_inches = 'tight')#saves figure as png)
    plt.close()
    
def timeconvert(time):
    hours = time.hour * 60
    totalmins = hours + time.minute
    totalmins = totalmins / 5
    if(totalmins == 0):
        totalmins = 5
    return int(totalmins)

def gifmake(date):
    #save 3 hour plots
    images = []
    i = 0
    while i < 1440:
        images.append(Image.open(str(int(i/60)) + ('%02d' % (i %60)) + '.png'))
        i = i + 15
        

    images[0].save('movie.gif',save_all=True, append_images=images[1:], optimize=True, loop=0, duration=200)

def threehourplot(date):
    images = []
    i = 0
    
    #select all images of 3 hour incriments
    while i < 23:
        images.append(Image.open(str(i) + '00.png'))
        i = i + 3
    size = images[0].size
    new_im = Image.new('RGB', (2*size[0], 4*size[1]))
    
    #paste all images and save as png
    new_im.paste(images[0], (0,0))
    new_im.paste(images[1], (size[0],0))
    new_im.paste(images[2], (0,size[1]))
    new_im.paste(images[3], (size[0],size[1]))
    new_im.paste(images[4], (0,2*size[1]))
    new_im.paste(images[5], (size[0],2*size[1]))
    new_im.paste(images[6], (0,3*size[1]))
    new_im.paste(images[7], (size[0],3*size[1]))
    new_im.save("merged_images.png", "PNG")
    
    
    
def deleteimages():
    i = 0
    while i < 1440:
        time = datetime.time(int(i/60), int(i % 60))
        os.remove(str(time.hour) + ('%02d' % time.minute) + '.png')
        i = i + 15
    
def pnggenerate(date):


    data = np.load('data.npz')
    tec = data['tec']
            
    i = 0
    while i < 1440:
        time = datetime.time(int(i/60), int(i % 60))
        timeint = timeconvert(time)
        Z = tec[:,:,timeint]
        Z = interpolate(Z, tec, timeint, 4)
                
        plot(Z, time, date)
        i = i + 15
    threehourplot(date)
    gifmake(date)
    deleteimages()

    
    
    
main()