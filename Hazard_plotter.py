# -*- coding: utf-8 -*-
"""
Created on 2022-March-29 

Author: John Malone-Leigh


The purpose of this script is to provide generate hazard maps of geoelectric 
fields using modelled geoelectric field time series.


The maps are generated following the method described in Malone-Leigh et al 2023

This script:

1) Loads in 25 year geoelectric field time series
2) Bins the time series based on time of of day
3) Reads in Kp time series to compare to the bins
4) Rebins again based on Kp indices
5) Generate hazard maps of using the binned time series. Three types:
    o The hazard maps used in Malone-Leigh et al 2023
    o The same hazard maps mapped using basemap, for ease of plotting
    o Rose plots, which illustrate the hazard maps for each individual site on
    on a single map
"""

from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
#import cv2
import os
from scipy import interpolate
import scipy.ndimage as ndimage
import mapping_geo_library as mpl
os.environ['PROJ_LIB'] = "" #have to set or basemap won't work
from mpl_toolkits.basemap import Basemap
plt.close('all')

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]



def plot_ireland2():
    #plots Ireland using Basemap
    fig, ax = plt.subplots() #needed to name plots
    """
    fontsize=8
    plt.rc('font',size=fontsize)
    
    plt.rc('legend',fontsize=fontsize)
    plt.rc('axes',labelsize=fontsize)
    """
    #'Ortho is ortho
    print('Plotting Map')
    #m= Basemap(projection='ortho',lat_0=45,lon_0=0,resolution='h')
    m=Basemap(projection='merc', llcrnrlat=51.0,llcrnrlon=-10.75,urcrnrlat=56,urcrnrlon=-5.25, resolution='l')
    m=Basemap(projection='mill', llcrnrlat=51.0,llcrnrlon=-10.75,urcrnrlat=56,urcrnrlon=-5.25, resolution='l')
    #setting miller projection to only show Europe on the map
    
    #m.drawmapboundary(fill_color='aqua')
    
    #m.drawlsmask(land_color='coral',ocean_color='aqua',lakes=True)
    m.drawcoastlines(color='black',linewidth=0.5,zorder=-20)
    
    #m.drawcountries()
    
    
    # draw parallels
    #m.drawparallels(np.arange(10,90,2),labels=[1,1,0,1],labelstyle='loosely dotted')
    # draw meridians
    #m.drawmeridians(np.arange(-180,180,2),labels=[1,1,0,1],labelstyle='loosely dotted') #note change third figure in np.arange to change how often lines are plotted

    return fig,ax,m

def read_co(path):
    
    """ Read name of the sites and coordinates of the site from the input
        files. Latitude and longitude should be in degrees.
		
		Parameters
		-----------
		path = path of the site with the name of the sites and coordinates

		Returns
		-----------
		name = Name of the site
		lat = latitude of the site (in degrees)
		lon = longitude of the site (in degrees)

		-----------------------------------------------------------------
    """
   
    a = pd.read_csv(path, 
                    header = None, 
                    skiprows = None, 
                    sep='\s+'
                    )
    
    a = np.array(a)
    name = a[:,0]
    lat =  a[:,1]
    lon =  a[:,2]
    
    return(name, lat, lon)


def map_coords3(sites,site_folder):
    #co_fold=r"C:\Users\johnn\OneDrive\Documents\Geo_Electrics_realtime_houdini\in/sites_interest2.dat"
    name,lat,lon=read_co(site_folder)

    for i in sites:
        if i=='SW05':
            sites==['VAL']
            i='VAL'
        leng=0
        for j in name:
            leng=leng+1
            if j==i:
                length=leng
        #coordinates in degrees
        lat1,lon1=(lat[length-1], lon[length-1])

        #converting to coordinate system for plot
        #ud1=(25*(coord[0]-bottom)/(top-bottom))#-radius 
        #rd1=(25*(coord[1]-left)/(right-left))#-radius
    return lat1,lon1

import matplotlib.colors as mcolors
import matplotlib.colors as colors

def make_cmap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1). For example: make_cmap([c('white'), c('cyan'), 0.10, c('cyan'), c('blue'), 0.50, c('blue'),\
                                      c('lime'), 0.90, c('lime')]) {the color map used in Weigt et al. (in prep.)
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,)*3]
    cdict = {'red':[], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)
#usign a custom map, note only used for Rose plots
c = colors.ColorConverter().to_rgb
custom_map = make_cmap([c('darkred'), c('red'), 0.10, c('red'), c('orange'), 0.50, c('orange'), c('yellow'), 0.90, c('yellow')])

plt.rcParams.update({'font.size': 16})

############################################
#Part 0 setting inputs

#select the mode of the hazard map. Three modes are available
# 1) "thres": the default mapping approach from Malone-Leigh et al which
# devides based on whether a threshold (i.e. 500mV/km) is exceeded, in each
# direction
# 2) "sum": instead of determining if a threhold is exceeded, sum simply 
# calculates the sum within the 3hr window, in each direction
# 3) "max": instead of determining if a threhold is exceeded, max 
# calculates the max within the 3hr window, in each direction
 
mode='thres'
#set your values for threshold here. If you want to analyse multiple at once,
#add them all to the list
#note this doesnt matter if using sum or max
thres_list=[100] 
#setting max values for images
if mode=='thres':
    vmax=1 #normalised
if mode=='sum':
    vmax=0.5*100000000
if mode=='max':
    vmax=0.5*1000


#set to 'yes' for galvanic distorted, 'no' for galvanic corrected
gal='yes' #
#set kp value to analyse
kp=9
#set using galvanic corrected 'no', uncorrected 'yes' or 
#... the difference
######################################################################
#Generating plots for hazard maps
fig5,ax5,m5=plot_ireland2()
fig6,ax6,m6=plot_ireland2()
fig7,ax7=plt.subplots()
fig8,ax8=plt.subplots()
fig9,ax9=plt.subplots()
plt.rc('font', size=20)

plt.rc('font', size=16)
plt.rc('legend',fontsize=16)
plt.rc('axes',labelsize=16)

#fig8,ax8,m8=plot_ireland2()
#fig9,ax9,m9=plot_ireland2()

######################################################################
#set the field to analyse, the Electric field ("Ex" and "Ey") or magnetic field 
#("Bx" and "By")
field1='Ex'
field2='Ey'
#set input_folder to the folder where 25 year time series are saved
folder1='data/'
#the file where coordinate information is saved
site_file='sites_interest.dat'
#set output folder, folder to save output images to
#default will save to same folder as the python script
folder2='' 
#adding lists to append final data to
zero_test_l=[]
xcos=[]
ycos=[]
map_colors=[]
vmaxs=[]
prob_lists=[]

for threshold in thres_list:
    if mode=='sum' or mode=='max':
        kp=0
        
    else:
        kp=kp
        #can manually set for different thresholds
    print('Kp: '+str(kp))
    #kp=5
    #name labels for sites
    sites=["I001",'SW02',"I003","I004","I006","I007","I009",
           "I013","I014","I026","I027","I101","I102","I103",
           "I105","I107","I108","SW01","SW02","SW03","SW04",
           "SW06","SW07","SW08","SW09","SW10","SW11","SW12","SW13",
           "SW14","SW15","SW16","SW17","SW18","SW19","SW20","SW21",
           "SW22",'i111','i114','i203','I104','I010','SW05',
           "i201", "I025","I002", "I110"]
    #sites=['I006','SW11','SW14','I003','i114','SW03']
    #uncomment second one for example with only a few sites
    seg=3*60 #3*60#segment length for 3 hours

    #######################################################################
    #loading data for each site
    max_vals=[]
    leng=0
    for site in sites:
        print(site)
        gals=[gal]	
        for gal in gals:
            #folder1=r'C:\Users\johnn\Downloads\Geo_Electrics\out/'
            string1='_'+field1+'_25.dat'
            
            string2='_'+field2+'_25.dat'
            if gal=='yes': #not galvanic correctred
                #folder1=r'C:\Users\johnn\Downloads\Save\gal/'
                string1=field1+'.dat'
                string2=field2+'.dat'
                
            if field1=='Bx' or field2=='By':
                #folder1=r'C:\Users\johnn\Downloads\Geo_Electrics\in\data\25_years/'
                string1='_bx_c.txt'
                string2='_by_c.txt'
            print('Loading Ex')
            Ex=np.loadtxt(folder1+site+string1,usecols=[0])
            print('Loading Ey')
            Ey=np.loadtxt(folder1+site+string2,usecols=[0])
            time_folder=folder1+'Time_v2.txt'
            print('Loading Hour data')
            hour=np.loadtxt(time_folder,usecols=3)
            print('Loading Kp data')
            
            kp_data=np.loadtxt(folder1+'Kp_ap_since_1932.txt',skiprows=172401,usecols=[7])
            
            #Getting rid of padding at edge of series
            if gal=='no':
                Ex=Ex[479610:-479610]
                Ey=Ey[479610:-479610]
            
            
            l0=0
            l_loop=0
            Ex2=[]
            Ey2=[]
            #############################################################
            print('Data Loaded')
            print('Kp selecting and Thresholding')
            for x,y in zip(Ex,Ey):
                
                if mode=='thres':
                    if kp_data[l0]==kp:
                        Ex2.append(x)
                        Ey2.append(y)
                if mode=='sum' or mode=='max':
                    #include all vals
                    if kp_data[l0]>=kp:
                        Ex2.append(x)
                        Ey2.append(y)                
                
                l_loop=l_loop+1
                if l_loop>=180: #containing data within loops 3 hour loops
                    l0=l0+1
                    l_loop=0
            
            #kp selected data
            Ex=Ex2
            Ey=Ey2
            #converting to amplitude and direction
            #convert to degrees from rad
            rad=360/(2*np.pi) 
            print('Getting Total')
            amp_list=[]
            angle_list=[]
            
            
            amp_list=[]
            angle_list=[]
            print('Amplitude and angle')
        
            for x,y in zip(Ex,Ey):
                amplitude=np.sqrt(x**2+y**2)
                if x>=0 and y>=0:
                    angle=math.atan(y/x)*rad
                elif x<0 and y>=0:
                    angle=-math.atan(y/x)*rad+90
                elif x<0 and y<0:
                    angle=math.atan(y/x)*rad+180
                elif x>=0 and y<0:
                    angle=-math.atan(y/x)*rad+270
                elif x==0: #to avoid divide by 0
                    if y>0:
                        angle=90
                    if y<0:
                        angle=270
                
                
                amp_list.append(amplitude)
                angle_list.append(angle)     
        
                
            amp_180=[]
            ang_180=[]
            
            for i in range(0,len(amp_list)-180,180):
                amp_180.append(amp_list[i:i+180])
                ang_180.append(angle_list[i:i+180])
            #now need to get probabities at each angle
            #splitting into 30 degree segments
            amp_ang_list=[]
            #for threshold in thres_list:
            print('Sorting 30 degrees')
            for j,k in zip(amp_180,ang_180):
                amp_ang_list2=[]
                for i in range(0,360,30):
                    l=0
                    degrees=[]
                    for l,m in zip(j,k):
                        if m>=i and m<i+30:
                        #if m>=0: #use this for summary map instead,
                        #each direction will yield summary instead
                            degrees.append(l)
                        l=l+1
                    
                    amp_ang_list2.append(degrees)           
            
                l2=[]
                for i in amp_ang_list2:
                    if mode=='thres':
                        if len(i)>0:
                            if max(i)>threshold:
                                l2.append(1)
                            else:
                                l2.append(0)
                        else:
                            l2.append(0)
                    if mode=='sum':
                        l2.append(np.sum(i))
                    if mode=='max':
                        try:
                            l2.append(np.max(i))
                        except:
                            #method passes except if no values are present
                            l2.append(0)
                amp_ang_list.append(l2)
                        
                #probabilties of exceeded threshold in each direction
            print('Probability each direction')
            prob_list=[]
            for i in range(0,12,1):
                deg_slice=np.array(amp_ang_list)[:,i]
                if mode=='thres':
                    prob=np.mean(deg_slice)
                if mode=='sum':
                    prob=np.sum(deg_slice)
                if mode=='max':
                    prob=np.max(deg_slice)
                prob_list.append(prob)
                    
            if mode=='thres':
                vmax=1
            else:
                vmaxs.append(np.max(prob_list))
            
            ############################################################
            #Now plotting ROSE Plots
            
            cmap = custom_map#plt.get_cmap("autumn",6) #ten colour segments
            colors = cmap(np.linspace(0.0,1,6)) 
            #colors[-1]=colors[-2] #turing white to last yellow
            
            value_or=0
            print('Plotting maps')
            color_list=[]
            for i in prob_list:
                #selecting colour based on prob for ROSE PLOTS
                try:
                    if mode=='thres':
                        color2=colors[int(100/6*np.round(i,1))]
                    if mode=='sum':
                        #ideally set factor to pre-tested max
                        factor=vmax #set to 1 for normalised
                        color2=colors[int((100/(6*factor))*np.round(i,1))]
                    if mode=='max':
                        factor=vmax #set to 1 for normalised
                        color2=colors[int((100/(6*factor))*np.round(i,1))]
                except:
                    color2=colors[-1]
                    #exceeds max threshold (rounding error), use last value
                #color2=colors[int(10**i)]
                
                if i<0.001:
                    #setting color on rose plots to be negigibly small
                    color2='white'
                
                color_list.append(color2)
                
        
                
                if i>0.01 and color2!='white':
                            
                    yco1,xco1=map_coords3([site])
                    #converting to basemap coordinates
                    xco,yco=m5(xco1,yco1)
    
                    ax5.scatter([xco],[yco],s=20,color='black')
    
                    if mode=='thres':
                        scale=5000*5*(3**i)
                    if mode=='sum' or mode=='max':
                        scale=i/1000
    
                    ax5.arrow(xco,yco,scale*np.cos(30*value_or*np.pi/180),scale*np.sin(30*value_or*np.pi/180),
                                  fc=color2, ec='black',zorder=10,head_width=int(scale)*0.5,linewidth=0.5)
                    
                    #np.pi/180 convertes to radians
                    #angle, xmin?,??,length of arrow
                else:
                    yco1,xco1=map_coords3([site])
                            #converting to basemap coordinates
                    xco,yco=m5(xco1,yco1)
                    ax5.scatter([xco],[yco],s=50,color='black',zorder=20)
                    
                    pass
                value_or=value_or+1
                if mode=='thres':
                    ax5.set_title('Probability > Threhold: '+str(threshold)+'mV/km')
    
                    
                attributes = [1,1,1,1,1,1,1,1,1,1,1,1]
    
                if max(prob_list)>0.001:
                    pass
                    
                    a = ax6.pie(attributes,
                        center=m6(xco1, yco1), 
                        colors = color_list,
                        # wedgeprops={'alpha':1},
                        radius= 0.2)
                else:
                    #plots a dot if no sites have large vals
                    ax6.scatter([xco],[yco],s=50,color='black',zorder=20)
                values, base = np.histogram(amp_list, bins=400)
                
                #####################################################
                #recording cumulative data for cumulative plot
                cumulative2=[]
                for i in range(len(values)):        
                    csum=sum(values[i:len(values)])
                    cumulative2.append(csum)


            max_vals.append(np.max(base))
            
            sites=['Site A Gal','Site A No Gal','Site B Gal','Site B No Gal',
                   'Site C Gal','Site C No Gal']
            colors=['blue','blue','orange','orange','green','green']
            ################################################################
            #plotting cumulative data for all sites            
            if gal=='yes':
                try:
                    ax7.plot(base[:-1], cumulative2,label=sites[leng],color=colors[leng])
                    plt.axvline(500,color='red',linestyle=':')
                except:
                    pass
            else:
                try:
                    ax7.plot(base[:-1], cumulative2,label=sites[leng],linestyle='--',color=colors[leng])
                    plt.axvline(500,color='red',linestyle=':')
                except:
                    pass
            leng=leng+1
        
            
            ax9.plot(base[:-1], values,label='Site '+str(leng)+' '+gal)
            plt.axvline(500,color='red',linestyle='--')
            xcos.append(xco1)
            ycos.append(yco1)
            map_colors.append(color_list)
            prob_lists.append(prob_list)
            
            del_cumulative=[]
            for i in range(0,len(cumulative2)-1):
                del_cumulative.append(abs(cumulative2[i]-cumulative2[i-1]))
            ax8.plot(base[1:-1], del_cumulative,label=site+' '+gal)
            plt.axvline(500,color='red',linestyle=':')
            ax8.set_xlabel('Threshold (mV/km)')
            ax8.set_ylabel('Counts')
            
    if mode=='thres':
        vmax=1 #not one hundred percent necessary +would already be counter acted
    else:
        #normalising plots by setting the max scale on colorbar to be max value
        vmax=np.max(vmaxs)
    #ax7.set_yscale('log')       
    #plt.rc('font', size=24)

    ax7.set_xlabel('Threshold (mV/km)',fontsize=16)
    plt.axvline(500,color='red',linestyle=':')
    ax7.set_ylabel('Counts',fontsize=16)
    plt.legend(fontsize=16)     
    #Now making summary heatmaps of direction and probabilities
    l=0

    plt.figure()
    #'Ortho is ortho
    print('Plotting Map')
    #####################################################################
    #Plotting onto map
    #First plotting to a Rose map
    

    m=Basemap(projection='cyl', llcrnrlat=51.0,llcrnrlon=-10.75,urcrnrlat=56,urcrnrlon=-5.25, resolution='l')
    
    #m.drawlsmask(land_color='coral',ocean_color='aqua',lakes=True)
    m.drawcoastlines(color='black',linewidth=0.5,zorder=-11)
    
    m.drawmeridians(np.arange(-180,180,2),labels=[1,1,0,1],labelstyle='loosely dotted') #note change third figure in np.arange to change how often lines are plotted
    
    for i in range(0,len(xcos)):
    
        a = plt.pie(attributes,
            center=m(xcos[i], ycos[i]),#has to be the same m or won't work! 
            colors = map_colors[i],
            # wedgeprops={'alpha':1},
            radius= 0.2,
            startangle=-15)
        #startangle 15 ensures first value East is centred around East
        
        tval=0
        for j in map_colors[i]:
            if j!='white':
                tval=tval+1
        #chcks if all values are white => prob=0
        #plotting point if no value present
        if tval==0:
            plt.scatter([xcos[i]],[ycos[i]],s=50,color='black',zorder=20)
    
    axis=plt.gca()
    axis.set_xlim([-11, -5]) 
    axis.set_ylim([51, 56])
    plt.title('Probability > Threhold: '+str(threshold)+'mV/km')
    
    save_no=10000+threshold
    
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    
    #plt.savefig(r"C:\Users\johnn\OneDrive\Documents\MagIE\MagIE Images/thres"+str(save_no)+'_k'+str(kp)+'.png',bbox_inches='tight')
    if mode=='thres':
        plt.savefig(folder2+"thres"+str(save_no)+'_k'+str(kp)+'.png',bbox_inches='tight')
    if mode=='sum':
        plt.savefig(folder2+"sum"+str(save_no)+'_k'+str(kp)+'.png',bbox_inches='tight')
    if mode=='max':
        plt.savefig(folder2+"max"+str(save_no)+'_k'+str(kp)+'.png',bbox_inches='tight')
    
    
    #labelling directions for summary heatmap
    dirs=['W - E','NEE - SWW','NNE - SSW','N - S','NNW - SSE','NWW-SEE']
    
    
    #numbers corresponding to sin 30 , 60 ,90 ...
    dir_no1=[0,15000,25980,30000,25980,15000]
    dir_no2=[30000,25980,15000,0,-15000,-25980]
    #################################################################
    #Now plotting Basemap hazard map
    for k in range(0,6):
        plt.figure()
        #m= Basemap(projection='ortho',lat_0=45,lon_0=0,resolution='h')
        m=Basemap(projection='merc', llcrnrlat=51.0,llcrnrlon=-10.75,urcrnrlat=56,urcrnrlon=-5.25, resolution='l')
        
        #setting miller projection to only show Europe on the map
        
        #m.drawmapboundary(fill_color='aqua')
        
        #m.drawlsmask(land_color='coral',ocean_color='aqua',lakes=True)
        m.drawcoastlines(color='black',linewidth=0.5,zorder=11)
        
        xco,yco=m(-9.75,55.5)
        scale=50000
        plt.arrow(xco,yco,scale*np.cos(30*value_or*np.pi/180),scale*np.sin(30*value_or*np.pi/180),
                              fc='white', ec='white',zorder=10,head_width=int(scale)*0.5,linewidth=4)
        heat_list=[]
        
        for i in prob_lists:
            heat_list.append(i[k]+i[k+6]) #+6 includes both directons
        
        
        x,y=xcos,ycos
        
        x=np.array(x)
        y=np.array(y)
        
        
        xi = np.linspace(-12, -5, 50)
        yi = np.linspace(50, 58, 50)
        
        xi, yi = np.meshgrid(xi, yi)
        xi,yi=m(xi,yi)
        
        
        
        zi = griddata(m(x, y), heat_list, (xi, yi),method='linear')
        nans, p = nan_helper(zi)
        zi[nans]= np.interp(p(nans), p(~nans), zi[~nans])

        plt.contourf(xi, yi, zi,vmin=0,vmax=vmax)

        label=dirs[k]
        if mode=='thres':
            if field1=='Ex':
                plt.title(label+' '+str(threshold)+' mV/km')
            elif field1=='Bx':
                plt.title(label+' '+str(threshold)+' nT')
            else:
                plt.title()
        if mode=='sum':
            plt.title(label+' Normalised Sum')
        if mode=='max':
            plt.title(label+' Normalised Max')

        plt.savefig(folder2+"Basemap_"+str(kp)+str(save_no)+'_'+str(k)+'.png',bbox_inches='tight')
        #plt.close()
        
        fig=plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        
        #m= Basemap(projection='ortho',lat_0=45,lon_0=0,resolution='h')
        #m=Basemap(projection='merc', llcrnrlat=51.0,llcrnrlon=-10.75,urcrnrlat=56,urcrnrlon=-5.25, resolution='h')
        
        #setting miller projection to only show Europe on the map
        
        #m.drawmapboundary(fill_color='aqua')
        
        #m.drawlsmask(land_color='coral',ocean_color='aqua',lakes=True)
        #m.drawcoastlines(color='black',linewidth=0.5,zorder=11)
        #shp_path_IRL = r'C:\Users\johnn\OneDrive\Documents\Geo_Electrics_realtime_houdini\in\data/Ireland_N&S.shp' 
        shp_path_IRL = 'Ireland_N&S.shp' 
        (xmin, 
          xmax, 
          ymin, 
          ymax, 
          r_cell_size, 
          alpha_v,
          title1) = mpl.inputs.get_predefined_regional_parameters('Ireland')
        mpl.imaging.create_background_map(fig,
                              ax,
                              shp_path_IRL, 
                              xmin, 
                              xmax, 
                              ymin, 
                              ymax, 
                              title1, 
                              )   
        lon,lat=m(xcos,ycos)

        lon=np.array(lon)
        lat=np.array(lat)
        numcols, numrows = 1000, 1000

        xi = np.linspace(lon.min(), lon.max(), numcols)
        yi = np.linspace(lat.min(), lat.max(), numrows)
        xi, yi = np.meshgrid(xi, yi)    
        
        z=np.array(heat_list)
        f = interpolate.interp2d(lon, lat, z, kind='linear')
        
        x2=np.arange(400000,900000,100000)
        y2=np.arange(500000,1100000,120000)
        #x2=np.arange(0,900000,100000)
        #y2=np.arange(0,1100000,120000)
        z=f(lon,lat)
        
        xi = np.linspace(x2.min(), x2.max(), 50)
        yi = np.linspace(y2.min(), y2.max(), 50)
        
        x,y,z=lon,lat,z
        xi, yi = np.meshgrid(xi, yi)
        z2=[]
        for i in heat_list:
            if i<0.01:
                i=0.01
            z2.append((abs(i)))
        #log value first before, causes probs otherwise
        blank=np.arange(0,len(xcos))
        #path to file containing data
        np.savetxt('data/coordinates_haz.dat',np.c_[blank,ycos,xcos])
        
        df_c = mpl.inputs.read_coordinates('data/coordinates_haz.dat', 
                                   lon_head = 'lon',
                                   lat_head = 'lat'
                                   )

        x=df_c.iloc[:,3]
        y=df_c.iloc[:,4]

        
        zi = griddata((x, y), heat_list, (xi, yi),method='linear')
        nans, p = nan_helper(zi)
        zi[nans]= np.interp(p(nans), p(~nans), zi[~nans])    
        #nans, x = nan_helper(zi)
        #zi[nans]= np.interp(x(nans), x(~nans), zi[~nans])
        #contourf handles cubic in mpl.preprocess below
        #applying a gaussian filter to smooth large artefacts of interpolation
        
        zi= ndimage.gaussian_filter(zi, 
                                   sigma=2.0, 
                                   order=0)
        #generating a mask to hide data in sea, where hazard map doesnt work
        mask_rain = mpl.pre_processing.generate_mask(zi, 
                                                     xi, 
                                                     yi, 
                                                     shp_path_IRL
                                                     )
        cm = plt.cm.get_cmap('viridis')
        
        #plotting data
        #note do not set vmin to 0
        #creates division by 0 issues when plotting
        var_1 = mpl.imaging.plot_background_data(fig,
                                             ax, 
                                             xi, 
                                             yi, 
                                             zi, 
                                             mask_rain,
                                             vmin = 0.01,
                                             vmax = vmax,
                                             cmap = cm
                                             )  
     
        #plotting map
        #fig, ax = plt.subplots() #needed to name plots
        fontsize=11
        plt.rc('font', size=20)
        
        if gal=='yes':
            label2='uncorrected'
        if gal=='no':
            label2='corrected'
        
        if mode=='thres':
            plt.title(label+' '+str(threshold)+' mV/km')
        if mode=='sum':
            plt.title(label+' Normalised Sum '+label2)
        if mode=='max':
            plt.title(label+' Normalised Max '+label2)
        #ax.set_facecolor('#8A959B')
        plt.tight_layout()
        #plotting directional arrows
        ax.arrow(470000,900000,dir_no2[k],dir_no1[k],color='black',zorder=20,head_width=20000)
        ax.arrow(470000,900000,-dir_no2[k],-dir_no1[k],color='black',zorder=20,head_width=20000)
        plt.title(str(threshold)+'mV/km, Kp'+str(kp))
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.scatter([x],[y],s=25,color='white',zorder=10)
        plt.grid(False)   
        plt.savefig('hmap'+str(k)+'.png')
