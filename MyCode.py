
import numpy as np
#import pygrib
import matplotlib.pyplot as pl
from netCDF4 import Dataset
#import iris
from mpl_toolkits import basemap
import os
FileNames=[]
MRMSFileList=[]
for file in os.listdir("/home/ho0man/Desktop/Win - Desktop/Jason PPT/Data/MRMS/NCs/"):
    if file.endswith(".nc"):
        MRMSFileList.append("/home/ho0man/Desktop/Win - Desktop/Jason PPT/Data/MRMS/NCs/"+ file)
        FileNames.append(file[:-3])
        
IMERGFileList=[]
for file in os.listdir("/home/ho0man/Desktop/Win - Desktop/Jason PPT/Data/IMERG/IMERG-Early/"):
    if file.endswith(".nc"):
        IMERGFileList.append("/home/ho0man/Desktop/Win - Desktop/Jason PPT/Data/IMERG/IMERG-Early/"+ file)    
for i in range(1):        
    #file = '/home/z5194283/MRMS.grib2' #example filename
    #dataset = pygrib.open(file)
    dataset1 = Dataset(MRMSFileList[i])
    dataset2 = Dataset(IMERGFileList[i])
       
    
#    filename = "C:/Users/Ho0maN/Desktop/MRMS/NCs/MRMSGC6H2018091400.nc"
    with Dataset(MRMSFileList[i], mode='r') as fh:
       print(fh.variables.keys())
       lons = fh.variables['longitude'][:]
       lats = fh.variables['latitude'][:]
    #   nlats = fh.variables['lat_0'][:]
    #   nplats=np.array(nlats)
    #   lats=np.flipud(nplats)
       sst = fh.variables['GaugeCorrQPE06H_0mabovemeansealevel'][:].squeeze()/6
    
#    filename2 = "C:/Users/Ho0maN/Desktop/IMERG/IMERG-Early/140000.nc"
    with Dataset(IMERGFileList[i], mode='r') as fh2:
       print(fh2.variables.keys())
       print(fh2.variables['lon'])
       lons2 = fh2.variables['lon'][:]
       lats2 = fh2.variables['lat'][:]
    #   nlats = fh.variables['lat_0'][:]
    #   nplats=np.array(nlats)
    #   lats=np.flipud(nplats)
       sst2 = fh2.variables['precipitationCal'][:]#.squeeze()
    
    lons_sub, lats_sub = np.meshgrid(lons2[::1], lats2[::1])
    
    sst_coarse = basemap.interp(sst, lons, lats, lons_sub, lats_sub, order=1)
#    sst_coarse[sst_coarse == -3] = 0
#    y,x=sst2.shape
    sss=np.transpose(sst2)
#    delta_sst=np.zeros((x,y))
    delta_sst=sst_coarse-sss   
#    image = pl.imshow(sst_coarse)
#    pl.show()
#    image2=pl.imshow(sss)
#    pl.show()
    
    #with Dataset(filename2) as fh3:
    #   fh3.variables['precipitationCal']=delta_sst
    #   fh3.close()
    
    
    #image3=pl.show(sst_coarse-sss)
    #pl.show()
    
#    f = Dataset('C:/Users/Ho0maN/Desktop/MRMS/NC-Regridded/'+FileNames[i]+'regrrided.nc','w', format='NETCDF4')
#    tempgrp = f.createGroup('Diff_data')
#    tempgrp.createDimension('lon', len(lats2))
#    tempgrp.createDimension('lat', len(lons2))
#    tempgrp.createDimension('time', None)
#    longitude = tempgrp.createVariable('Longitude', 'f4', 'lon')
#    latitude = tempgrp.createVariable('Latitude', 'f4', 'lat') 
#    time = tempgrp.createVariable('Time', 'i4', 'time')
#    temp = tempgrp.createVariable('Diff_data', 'f4', ('time','lon', 'lat'))
#    longitude[:] = lats2 #The "[:]" at the end of the variable instance is necessary
#    latitude[:] = lons2
#    temp[0,:,:] = sst_coarse
#    f.close()
#    
#    f1 = Dataset('C:/Users/Ho0maN/Desktop/Differences/'+FileNames[i]+'Delta.nc','w', format='NETCDF4')
#    tempgrp = f1.createGroup('Diff_data')
#    tempgrp.createDimension('lon', len(lats2))
#    tempgrp.createDimension('lat', len(lons2))
#    tempgrp.createDimension('time', None)
#    longitude = tempgrp.createVariable('Longitude', 'f4', 'lon')
#    latitude = tempgrp.createVariable('Latitude', 'f4', 'lat') 
#    time = tempgrp.createVariable('Time', 'i4', 'time')
#    temp = tempgrp.createVariable('Diff_data', 'f4', ('time','lon', 'lat'))
#    longitude[:] = lats2 #The "[:]" at the end of the variable instance is necessary
#    latitude[:] = lons2
#    temp[0,:,:] = delta_sst
#    f1.close()    

    
toexclude = ['ExcludeVar1', 'ExcludeVar2']

with Dataset(IMERGFileList[0], mode='r') as src, Dataset('/home/ho0man/Desktop/Win - Desktop/Jason PPT/Data/MRMS/NC-Regridded/'+FileNames[0]+'regrrided.nc','w', format='NETCDF4', 'w') as dst:
    # copy global attributes all at once via dictionary
    dst.setncatts(src.__dict__)
    # copy dimensions
    for name, dimension in src.dimensions.items():
        dst.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))
    # copy all file data except for the excluded
    for name, variable in src.variables.items():
        if name not in toexclude:
            x = dst.createVariable(name, variable.datatype, variable.dimensions)
            dst[name][:] = src[name][:]
            # copy variable attributes all at once via dictionary
#            dst[name].setncatts(src[name].__dict__)

#f = Dataset('/home/ho0man/Desktop/Win - Desktop/Jason PPT/Data/MRMS/NC-Regridded/'+FileNames[0]+'regrrided.nc','w', format='NETCDF4')
#tempgrp = f.createGroup('MRMS-Regrided')
#tempgrp.createDimension('lon', len(lons2))
#tempgrp.createDimension('lat', len(lats2))
#tempgrp.createDimension('time', None)
#longitude = tempgrp.createVariable('Longitude',  np.float32, 'lon')
#latitude = tempgrp.createVariable('Latitude',  np.float32, 'lat') 
#time = tempgrp.createVariable('Time', 'i4', 'time')
#temp = tempgrp.createVariable('Precip. Regrid', np.float32, ('time','lon', 'lat'))
##temp.coordinates="lon lat"
#longitude[:] = lons2#The "[:]" at the end of the variable instance is necessary
#latitude[:] = lats2
##
##longitude.units = 'degrees east'
##latitude.units = 'degrees north'
#temp[0,:,:] = np.transpose(sst_coarse)
##temp[0,:,:] = sst_coarse
#f.close()

#IRISLoadedMRMS = iris.load_cube('/home/z5194283/MRMS.nc')


#IRISLoadedIMERG = iris.load_cube('/home/z5194283/IMERG.nc.nc')

#for i in dataset:
#    for j in i.keys():
#        print(j)
#print(dataset[1].latitudes)
#np.savetxt('test.txt', dataset[1].values, delimiter=',')
#print(dataset[1].latitudeOfFirstGridPointInDegrees, '<--LATF')
#print(dataset[1].longitudeOfFirstGridPointInDegrees, '<--LONF')
#print(dataset[1].latitudeOfLastGridPointInDegrees, '<--LATL')
#print(dataset[1].longitudeOfLastGridPointInDegrees, '<--LONL')

#LATF=dataset[1].latitudeOfFirstGridPointInDegrees
#LONF=dataset[1].longitudeOfFirstGridPointInDegrees
#LATL=dataset[1].latitudeOfLastGridPointInDegrees
#LONL=dataset[1].longitudeOfLastGridPointInDegrees


#print(dataset2.variables['lat'][0])
#print(dataset2.variables['lon'][0])
#print(dataset2.variables['lat'][140])
#print(dataset2.variables['lon'][200])



#Dlat=(dataset2.variables['lat'][1]-dataset2.variables['lat'][0])
#Dlon=(dataset2.variables['lon'][1]-dataset2.variables['lon'][0])
#print(np.shape(dataset[1].distinctLatitudes))
#print(dataset[1].distinctLatitudes[0])
#print(dataset[1].distinctLongitudes[0])
#print(dataset[1].distinctLatitudes[3499])
#print(dataset[1].distinctLongitudes[6999])
#y,x=np.shape(dataset[1].values)
#x1,y1=np.shape(dataset2.variables['precipitationCal'])
#MRMSave=np.zeros((x1, y1))


#IRISRegridedMRMS=IRISLoadedIMERG.regrid(IRISLoadedMRMS, iris.analysis.linear())

#for ii in range(x1):
#    for jj in range(y1):
#        k=0
#        for i in range(x):
#            for j in range(y):
#                if dataset[1].distinctLongitudes[i]<dataset2.variables['lat'][$
#                    k=k+1
#                    MRMSave[i,j] = MRMSave[i,j]+dataset[1].values[i,j]
#        MRMSave[i,j]/k
#                    print(MRMSave[i,j])
#np.savetxt('MRMSave.txt', MRMSave, delimiter=',')
#np.savetxt('DiffSave.txt',MRMSave/dataset2.variables['precipitationCal'],delim$
#image = pl.imshow(IRISRegridedMRMS)
#pl.show()



#print(No0fLonPix1)


#dataset2 = Dataset("/home/hooman/Desktop/MRMS/IMERG/3B-DAY-E.MS.MRG.3IMERG.201$
#print(dataset2.variables.keys())
#print(dataset2.variables['precipitationCal'])
#print(dataset2.variables['lat'][2])
#print(dataset2.variables['lon'][2])
#y1,x1=np.shape(dataset2.variables['precipitationCal'])
#print(y1,x1)
#image = pl.imshow(dataset2.variables['precipitationCal'])
#pl.show()


#image = pl.imshow(dataset[1].latLonValues)
#pl.show()
#image = pl.imshow(dataset[1].values)
#pl.show()




#image = pl.imshow(dataset[1])
#pl.show()
#print inventory



#print(dataset.variables.keys())
#print(dataset.variables['precipitationCal'])
#print

#image = pl.imshow(dataset)#.variables['precipitationCal'])
#pl.show()




