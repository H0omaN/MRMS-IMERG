import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.colors as mcolors
from scipy import stats

import numpy as np
import netCDF4 as nc
#import pygrib

from netCDF4 import Dataset
#import iris
from mpl_toolkits import basemap
import os
FileNames=[]
MRMSFileList=[]
#motherlink = '/home/z5194283/Mygit/'
motherlink = '/home/ho0man/Mygit/'
for file in os.listdir(motherlink+"Data/MRMS/NC-Regridded/"):
    if file.endswith(".nc"):
        MRMSFileList.append(motherlink+"Data/MRMS/NC-Regridded/"+ file)
        FileNames.append(file[:-3])
MRMSFileList.sort()        
IMERGFileList=[]
for file in os.listdir(motherlink+"Data/IMERG/IMERG-Early/"):
    if file.endswith(".nc"):
        IMERGFileList.append(motherlink+"Data/IMERG/IMERG-Early/"+ file)    
IMERGFileList.sort()
for i in range(4):        
    #file = '/home/z5194283/MRMS.grib2' #example filename
    #dataset = pygrib.open(file)
    dataset1 = Dataset(MRMSFileList[i])
    dataset2 = Dataset(IMERGFileList[i])   

    
#    filename = "C:/Users/Ho0maN/Desktop/MRMS/NCs/MRMSGC6H2018091400.nc"
    with Dataset(MRMSFileList[i], mode='r') as fh:
       print(fh.variables.keys())
       lons0 = fh.variables['lon'][:]
       lats0 = fh.variables['lat'][:]
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
       sst2=np.transpose(sst2)
    
#    lons_sub, lats_sub = np.meshgrid(lons2[::1], lats2[::1])
#    
#    sst_coarse = basemap.interp(sst, lons0, lats0, lons_sub, lats_sub, order=1)
#    sst_coarse[sst_coarse == -3] = 0
#    y,x=sst2.shape    
#    delta_sst=np.zeros((x,y))
      
#    x=np.percentile(sst2,stats.percentileofscore(sst, for a in sst))   
    Qdiff=np.zeros_like(sst)
#    sstrank = sst.rank(pct=True)   
#    sst2rank= sst2.quantile(sstrank)
    sst2Raveled = sst2.ravel()
    for ii in range(len(sst2)):
        for jj in range(len(sst2[ii])):
            sst2Qvalue=stats.percentileofscore(sst2Raveled, sst2[ii,jj])/100
            sstEQ=np.quantile(sst,sst2Qvalue)
            Qdiff[ii,jj]=sstEQ-sst2[ii,jj]
    

#    delta_sst=sst-sst2  
    
#######plotting##############
    Map=[Qdiff,sst,sst2]
    Map_names=['Difference-Quantile','MRMS-Regridded','IMERG']
    for j in range(3):  
        lat = lats2
        lon = lons2
        #u = nc.variables['GaugeCorrQPE06H_0mabovemeansealevel'][:]
        u = Map[j]    
        
        minlat=np.min(lat)
        maxlat=np.max(lat)
        minlon=np.min(lon)
        maxlon=np.max(lon)
        
        map = Basemap(projection='merc',llcrnrlon=minlon,llcrnrlat=minlat,urcrnrlon=maxlon,urcrnrlat=maxlat,resolution='i') # projection, lat/lon extents and resolution of polygons to draw
        # resolutions: c - crude, l - low, i - intermediate, h - high, f - full

#        map.drawstates(antialiased=0.1)
        #map.drawlsmask(land_color='Linen', ocean_color='#CCFFFF') # can use HTML names or codes for colors
        #map.drawcounties()
        parallels = np.arange(minlat,maxlat,10.) # make latitude lines ever 5 degrees from 30N-50N
        meridians = np.arange(minlon,maxlon,10.) # make longitude lines every 5 degrees from 95W to 70W
        map.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
        map.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
        
        
        
        lons,lats= np.meshgrid(lon,lat) # for this dataset, longitude is 0 through 360, so you need to subtract 180 to properly display on map
        x,y = map(lons,lats)
        #a=np.max(u[0,:,:])
        a=np.max(u)
#        clevs = np.arange(0,a,(a-0)/20)
        cmap_data = [(1.0, 1.0, 1.0),
                     (0.3137255012989044, 0.8156862854957581, 0.8156862854957581),
                     (0.0, 1.0, 1.0),
                     (0.0, 0.8784313797950745, 0.501960813999176),
                     (0.0, 0.7529411911964417, 0.0),
                     (0.501960813999176, 0.8784313797950745, 0.0),
                     (1.0, 1.0, 0.0),
                     (1.0, 0.6274510025978088, 0.0),
                     (1.0, 0.0, 0.0),
                     (1.0, 0.125490203499794, 0.501960813999176),
                     (0.9411764740943909, 0.250980406999588, 1.0),
                     (0.501960813999176, 0.125490203499794, 1.0),
                     (0.250980406999588, 0.250980406999588, 1.0),
                     (0.125490203499794, 0.125490203499794, 0.501960813999176),
                     (0.125490203499794, 0.125490203499794, 0.125490203499794),
                     (0.501960813999176, 0.501960813999176, 0.501960813999176),
                     (0.8784313797950745, 0.8784313797950745, 0.8784313797950745),
                     (0.9333333373069763, 0.8313725590705872, 0.7372549176216125),
                     (0.8549019694328308, 0.6509804129600525, 0.47058823704719543),
                     (0.6274510025978088, 0.42352941632270813, 0.23529411852359772),
                     (0.4000000059604645, 0.20000000298023224, 0.0)]
        cmap = mcolors.ListedColormap(cmap_data, 'precipitation')
        #cs = map.contourf(x,y,u[0,:,:],clevs,cmap=cm.s3pcpn)#colors='blue')#,linewidths=1.)
        #cs = map.imshow(u[0,:,:],cmap=cmap)# cmap=cm.s3pcpn)#, alpha = 0.5)
        
        cs = map.imshow(u,cmap=cmap)# cmap=cm.s3pcpn)#, alpha = 0.5)
        #plt.clabel(cs, fontsize=9, inline=1) # contour labels
        map.drawcoastlines(linewidth=0.2)
        map.drawstates(linewidth=0.2)
#        map.drawcountries(antialiased=0.1)
        cbar = map.colorbar(cs,location='right',pad="5%")
        cbar.set_label('mm/hr')
#        plt.set_xlim(x-0.25, x+0.25)
#        plt.set_ylim(y-0.25, y+0.25)
        plt.title(Map_names[j]+str(140000+600*i)+' - Precpitaiton')     
        plt.savefig('/home/ho0man/Mygit/MRMS-IMERG/'+str(140000+600*i)+'/'+Map_names[j]+str(140000+600*i)+'.png',dpi=800)
        plt.clf()

        

