# -*- coding: utf-8 -*-
"""
Created on Thu May 21 22:59:01 2020

@author: joanc
"""
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
#import cv2
#import os
import numpy as np
#import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
import matplotlib
#import time


#from os.path import isfile
#from os.path import join
from scipy.interpolate import griddata#, RectBivariateSpline
#from sklearn.gaussian_process import GaussianProcessRegressor
from osgeo import  osr # gdal, ogr,
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#from matplotlib.colors import ListedColormap#, LinearSegmentedColormap
from datetime import timedelta # datetime, 
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#from multiprocessing import Pool
#from functools import partial
#from matplotlib import gridspec



#######################################################################
#######################################################################
class inputs():    
    ###################################################################    
    def read_data(path, 
                  s_d, 
                  e_d, 
                  time_format = '%Y/%m/%d',
                  time_param_format = '%Y-%m-%d',
                  ind_col = 0
                  ):
    

        df = pd.read_csv(path,
                         index_col = ind_col
                         )
       
        df = df.loc[pd.notnull(df.index)]
        df["index"] = df.index
           
        min_date = pd.to_datetime(df.index).min()
        max_date = pd.to_datetime(df.index).max()
       
        #if date range is within the max/min range, use it. Otherwise ignore it (max/min from hydrographs will be used as default)
        if pd.to_datetime(s_d) > min_date or pd.to_datetime(e_d) < max_date:
            # Trim the dataframe by the specified dates
            df_2 = df[(pd.to_datetime(df['index'],
                                      format = time_format
                                      )
                       >=
                       pd.to_datetime(s_d,
                                      format = time_param_format
                                      ))
                       &
                       (pd.to_datetime(df['index'],
                                       format = time_format
                                       )
                       <
                       pd.to_datetime(e_d,
                                      format = time_param_format
                                      )
                       )]
                       
        else:
            df_2 = df
    
             
        sites = [x for x in df.columns]
        sites = sites[0:-1]
        
        return(df_2, 
               sites
               )
    
    #######################################################################
    def read_coordinates(path, 
                         lon_head = 'X_WGS84',
                         lat_head = 'Y_WGS84',
                         espg = 4326,
                         to_ITM = True
                         ):
        
        # 1.3 Read Coordinates
        df = pd.read_csv(path,
                         header=None,
                         names=['lat', 'lon'],
                         usecols=[1,2],
                         delimiter=' ')
        lon_head='lon'
        lat_head='lat'
        # Convert to ITM
        if to_ITM == True:
            geometry_sar = [Point(xy) for xy in zip(df[lon_head], 
                                                    df[lat_head]
                                                    )]
            crs = {'init': 'epsg:' + str(espg)}
            
            gdf = gpd.GeoDataFrame(df, 
                                   crs=crs, 
                                   geometry=geometry_sar
                                   )
        
            gdf = gdf.to_crs(epsg=2157)
            
            df['ITM_E'] = gdf['geometry'].x 
            df['ITM_N'] = gdf['geometry'].y  
    
        return(df)    
            
    #######################################################################
    def get_predefined_regional_parameters(area_int): 
    
        if area_int == 'Ireland':
            xmin = 410000
            xmax = 770000
            ymin = 510000
            ymax = 970000
            
            title1 = 'Ireland'
            
            r_cell_size = 3000
            alpha_v = 0.75   
            
        return(xmin, 
               xmax, 
               ymin, 
               ymax, 
               r_cell_size, 
               alpha_v,
               title1) 

#######################################################################
#######################################################################      
class pre_processing():
    #######################################################################
    '''
    def create_interp_rainfall_values(sites, 
                                      rain, 
                                      rain_cc, 
                                      xmin, 
                                      xmax, 
                                      ymin, 
                                      ymax,  
                                      step,
                                      date,
                                      delta_t = 0
                                      ):
    
        # Define times of interest
        d_min = pd.to_datetime(date) - timedelta(days = delta_t)
        d_max = pd.to_datetime(date) + timedelta(days = 0)
        
        # Create new grid 
        xmin = xmin - step
        xmax = xmax + step
        ymin = ymin - step
        ymax = ymax + step
        
        x_new = []
        y_new = []
        for i in range(xmin, xmax, step ):
            x_new.append(i)
        for i in range(ymin, ymax, step ):
            y_new.append(i)
        
        new_grid_x, new_grid_y = np.meshgrid(x_new, 
                                             y_new
                                             )
    
        # Get data of interest from met eirean
        z_t1 = []
        x_t1 = []
        y_t1 = []
        
        av_rain = rain[(pd.to_datetime(rain.index,
                              format = '%Y-%m-%d'
                              )
                        >=
                        d_min)
                      &
                       (pd.to_datetime(rain.index,
                               format = '%Y-%m-%d'
                               )
                         <=
                         d_max
                         )]
        
        av_rain = av_rain.replace([' '], 
                                  np.nan
                                  )
        
        av_rain = av_rain.astype(float)
        for i in sites[1::]:    
            for j in rain_cc['# site_number']:
                a = int(j)
                try:
                    b = int(i)
                except:
                    b = int(i[1:-1])
                if a == b:
                    if av_rain[i].isnull().values.all() == False: 
                        z_t1.append(np.nanmean((av_rain[i])))
                        x_t1.append(rain_cc['ITM_E'][rain_cc['# site_number'] 
                                                     == a].iloc[0]
                                    )
        
                        y_t1.append(rain_cc['ITM_N'][rain_cc['# site_number'] 
                                                     == a].iloc[0]
                                    )
        
        x_t = []
        y_t = []
        z_t = []
        for i in range(0,len(z_t1)):
            if np.isnan(z_t1[i]) != True:
                y_t.append(y_t1[i])
                x_t.append(x_t1[i])
                z_t.append(z_t1[i])
        
        points = np.column_stack([x_t, y_t])
            
        zg_l = griddata(np.array(points), 
                        np.array(z_t), 
                        (new_grid_x.astype('float'), 
                         new_grid_y.astype('float')), 
                        #method = 'cubic'
                        #method = 'linear'
                        method = 'nearest'
                        )
                        
        zg_l = ndimage.gaussian_filter(zg_l, 
                                       sigma=5.0, 
                                       order=0)
        
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(2157)
        y, x = np.mgrid[slice(ymin, 
                              ymax, 
                              step
                              ),
                        slice(xmin, 
                              xmax, 
                              step
                              )
                        ]
        
        return(x, 
               y, 
               zg_l
               )    
    
    ###################################################################
    def select_data_in_region_of_interest(df1, 
                                          df2, 
                                          xmin, 
                                          xmax, 
                                          ymin, 
                                          ymax
                                          ):
        """ df1: dataframe with time series
            df2: dataframe with coordinates of the sites in ITM
            xmin,xmax,ymin,ymax coodinates of the box in ITM
        """
        remove_list = []
        for i in range(0,len(df2)):
            for j in range(0,len(df1.columns)):
                if df2.index[i] == df1.columns[j]:
                    if df2['ITM_E'][i] < xmin:
                        remove_list.append(df2.index[i])
                        break
                    if df2['ITM_E'][i] > xmax:
                        remove_list.append(df2.index[i])
                        break
                    if df2['ITM_N'][i] < ymin:
                        remove_list.append(df2.index[i])
                        break
                    if df2['ITM_N'][i] > ymax:
                        remove_list.append(df2.index[i]) 
                        break
        
        for i in remove_list:
            df2 = df2.drop(index = [i])
            df1 = df1.drop(columns = [i])   
    
        return(df1, 
               df2)'''
    def create_interp_rainfall_values(sites, 
                                      rain, 
                                      rain_cc, 
                                      xmin, 
                                      xmax, 
                                      ymin, 
                                      ymax,  
                                      step,
                                      delta_t = 0
                                      ):
    
        
        # Create new grid 
        xmin = xmin - step
        xmax = xmax + step
        ymin = ymin - step
        ymax = ymax + step
        
        x_new = []
        y_new = []
        for i in range(xmin, xmax, step ):
            x_new.append(i)
        for i in range(ymin, ymax, step ):
            y_new.append(i)
        
        new_grid_x, new_grid_y = np.meshgrid(x_new, 
                                             y_new
                                             )
    
        # Get data of interest from met eirean
        z_t1 = []
        x_t1 = []
        y_t1 = []
        
        

        for i in sites[1::]:    
            for j in rain_cc['# site_number']:
                a = int(j)
                try:
                    b = int(i)
                except:
                    b = int(i[1:-1])
                if a == b:
                    if av_rain[i].isnull().values.all() == False: 
                        z_t1.append(np.nanmean((av_rain[i])))
                        x_t1.append(rain_cc['ITM_E'][rain_cc['# site_number'] 
                                                     == a].iloc[0]
                                    )
        
                        y_t1.append(rain_cc['ITM_N'][rain_cc['# site_number'] 
                                                     == a].iloc[0]
                                    )
        
        x_t = []
        y_t = []
        z_t = []
        for i in range(0,len(z_t1)):
            if np.isnan(z_t1[i]) != True:
                y_t.append(y_t1[i])
                x_t.append(x_t1[i])
                z_t.append(z_t1[i])
        
        points = np.column_stack([x_t, y_t])
            
        zg_l = griddata(np.array(points), 
                        np.array(z_t), 
                        (new_grid_x.astype('float'), 
                         new_grid_y.astype('float')), 
                        #method = 'cubic'
                        #method = 'linear'
                        method = 'nearest'
                        )
                        
        zg_l = ndimage.gaussian_filter(zg_l, 
                                       sigma=5.0, 
                                       order=0)
        
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(2157)
        y, x = np.mgrid[slice(ymin, 
                              ymax, 
                              step
                              ),
                        slice(xmin, 
                              xmax, 
                              step
                              )
                        ]
        
        return(x, 
               y, 
               zg_l
               )          

    #######################################################################
    def interpolation_2d(sites, 
                         rain, 
                         rain_cc, 
                         xmin, 
                         xmax, 
                         ymin, 
                         ymax,  
                         step,
                         date,
                         delta_t = 0
                         ):
        
        # Define times of interest
        d_min = pd.to_datetime(date) - timedelta(days = delta_t)
        d_max = pd.to_datetime(date) + timedelta(days = 1)
        
        # Create new grid 
        xmin = xmin - step
        xmax = xmax + step
        ymin = ymin - step
        ymax = ymax + step
        
        x_new = []
        y_new = []
        for i in range(xmin, xmax, step ):
            x_new.append(i)
        for i in range(ymin, ymax, step ):
            y_new.append(i)
        
        new_grid_x, new_grid_y = np.meshgrid(x_new, 
                                              y_new
                                              )
        
        # Get data of interest from met eirean
        z_t1 = []
        x_t1 = []
        y_t1 = []
        
        av_rain = rain[(pd.to_datetime(rain.index,
                              format = '%Y-%m-%d'
                              )
                        >=
                        d_min)
                      &
                        (pd.to_datetime(rain.index,
                                format = '%Y-%m-%d'
                                )
                          <=
                          d_max
                          )]
        
        av_rain = av_rain.replace([' '], 
                                  np.nan
                                  )
        
        av_rain = av_rain.astype(float)
        for i in sites[1::]:    
            for j in rain_cc['# site_number']:
                a = int(j)
                try:
                    b = int(i)
                except:
                    b = int(i[1:-1])
                if a == b:
                    if av_rain[i].isnull().values.all() == False: 
                        z_t1.append(np.nanmean((av_rain[i])))
                        x_t1.append(rain_cc['ITM_E'][rain_cc['# site_number'] 
                                                      == a].iloc[0]
                                    )
        
                        y_t1.append(rain_cc['ITM_N'][rain_cc['# site_number'] 
                                                      == a].iloc[0]
                                    )
        
        x_t = []
        y_t = []
        z_t = []
        for i in range(0,len(z_t1)):
            if np.isnan(z_t1[i]) != True:
                y_t.append(y_t1[i])
                x_t.append(x_t1[i])
                z_t.append(z_t1[i])
        
        points = np.column_stack([x_t, y_t])
            
        zg_l = griddata(np.array(points), 
                        np.array(z_t), 
                        (new_grid_x.astype('float'), 
                          new_grid_y.astype('float')), 
                        #method = 'cubic'
                        #method = 'linear'
                        method = 'nearest'
                        )
        
        zg_l = ndimage.gaussian_filter(zg_l, 
                                       sigma=5.0, 
                                       order=0)
        
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(2157)
        y, x = np.mgrid[slice(ymin, 
                              ymax, 
                              step
                              ),
                        slice(xmin, 
                              xmax, 
                              step
                              )
                        ]
        
        return(x, 
                y, 
                zg_l
                )
        
    #######################################################################
    def generate_mask(z_rain, 
                      x_rain, 
                      y_rain, 
                      shp_path_IRL_u,
                      atribute = 'FID',
                      #at_val = 6538,
                      at_val = 1,
                      simplify = 10000,
                      buffer = 10000
                      ):
    
        RIug = gpd.read_file(shp_path_IRL_u) 
        #RIug = RIu.unary_union  
        #print(RIu)                             
        #RIug = RIu[RIu[atribute] == at_val]
        RIug.crs =  {'init': 'epsg:2157'}
        gdf = RIug['geometry'].simplify(simplify, 
                                        preserve_topology=True
                                        ).buffer(buffer)
    
        
        mask = np.zeros(z_rain.shape) #change to zeros
        
        # Get values of rain to be considered
        for i in range(0, z_rain.shape[0]):
            for j in range(0, z_rain.shape[1]):
                pg = Point(x_rain[i,j], y_rain[i,j])
                g1 = gpd.GeoSeries([pg])
                if gdf.iloc[0].contains(g1[0]) == True:
                    mask[i,j] = 1 # northern Ireland
                if gdf.iloc[1].contains(g1[0]) == True:
                    mask[i,j] = 1#Ireland
    
        return(mask)    
    
#######################################################################
#######################################################################
class imaging():
    def def_name_based_reference_date(t, t_ref = '2000-01-01'):
        name = (t - pd.to_datetime(t_ref) ).days
        return(name)

    ###################################################################
    def add_title(fig, tit, i):
        font = matplotlib.font_manager.FontProperties(family = 'times new roman', 
                                                      size = 20
                                                      )   
        title = fig.suptitle(str(tit)
                             + ' '
                             + str(i)[0:10], 
                             x = 0.4,
                             y = 0.98
                             )
    
        title.set_font_properties(font)
        return()
    ###################################################################
    def plot_background_data(fig,
                             ax,
                             x_rain, 
                             y_rain, 
                             z_rain, 
                             mask_rain,
                             vmin = 0,
                             vmax = 15,
                             cmap = 'OrRd'
                             ): 
        
        
        z_rain_t = np.where(mask_rain == 0., 
                            np.nan, 
                            z_rain
                            )
        
    
        # new_cmap1 = plt.cm.get_cmap('Blues', 512)
        # new_cmap = ListedColormap(new_cmap1(np.linspace(0, 
        #                                                 0.6, 
        #                                                 256
        #                                                 )))

        '''
        var_1 = ax.pcolormesh(x_rain,
                               y_rain,
                               z_rain,
                               vmin = .1,
                               vmax = 6.9,
                               cmap = cmap,
                               zorder = 5,
                               shading = 'gouraud'
                               #alpha = 0.75,

                               )
        '''
        
        
        var_1 = ax.contourf(x_rain, 
                         y_rain,
                         z_rain_t, 
                         zorder=2, 
                         alpha=0.9,
                         cmap=cmap,
                         vmin=vmin,
                         vmax=vmax,
                         interpolation='linear'
                         )
        
       
   
        #print('max_rain:' + str(np.nanmax(z_rain_t)))             
      
        return(var_1)       
    #######################################################################        
    def add_shpfile_info(fig, 
                         ax, 
                         shp_path_IRL 
                         ):
            
        # Get maps of Ireland
        RI = gpd.read_file(shp_path_IRL)    
            
        RI.plot(ax = ax, 
                facecolor="none", 
                edgecolor='black',
                lw = 0.5,
                zorder = 10)
        
        return()    

    #######################################################################
    def set_area_of_interest(fig, 
                             ax1, 
                             xmin, 
                             xmax, 
                             ymin, 
                             ymax
                             ):
        
        ax1.set(xlim=(xmin, 
                      xmax), 
                ylim=(ymin, 
                      ymax
                      )
                )
        
        return()

    #######################################################################   
    def set_axes_format_and_area_of_interest(fig, 
                                             ax1
                                             ):
        # Set font for title
        font = matplotlib.font_manager.FontProperties(family = 'times new roman', 
                                                      size = 15
                                                      )   

    
        ax1.title.set_font_properties(font)
        #ax1.set_facecolor('white')
        
        # Deffine labels
        '''
        ax1.set_xlabel('Easting', 
                       fontsize=10, 
                       color = 'gray'
                       )
        
        ax1.set_ylabel('Northing', 
                       fontsize=10, 
                       rotation = 90, 
                       color = 'gray'
                       )
        
        ax1.tick_params(axis="x", 
                        labelsize=5, 
                        labelcolor = 'gray', 
                        labelrotation = 0, 
                        pad = 5
                        )
        
        ax1.tick_params(axis="y", 
                        labelsize=5, 
                        labelcolor = 'gray', 
                        labelrotation = 45, 
                        pad = 5
                        )'''
        return()   

    
    #######################################################################
    def create_background_map(fig,
                          ax1,
                          shp_path_IRL, 
                          xmin, 
                          xmax, 
                          ymin, 
                          ymax, 
                          title1, 
                          ):
    
        # 5.1 Adding counties
        print('    ...adding counties...')
        imaging.add_shpfile_info(fig,
                                 ax1,
                                 shp_path_IRL
                                 )
        
        # 5.4 Zoom area of interest
        print('    ...zooming to area of interest...')
        imaging.set_area_of_interest(fig, 
                                     ax1, 
                                     xmin, 
                                     xmax, 
                                     ymin, 
                                     ymax
                                     )
        
        # 5.5 Set the axes of the figure 
        print('    ...setting axes...')
        imaging.set_axes_format_and_area_of_interest(fig, 
                                                      ax1 
                                                      )
        
        return()
 
    #######################################################################
    def add_site_data(fig, 
                      ax1, 
                      i, 
                      df_h, 
                      df_c,
                      cmap = 'OrRd'
                      ):
        
        if i in df_h.index:    
            sel_h = df_h[df_h.index == i]
            df_c['values'] = np.zeros(len(df_c)) 
            for i in range(0,len(df_c)):
                try:
                    df_c['values'][i] = sel_h[str(df_c['# site_number'][i])]
                except:
                    df_c['values'][i] = np.nan

            var_2 = ax1.scatter(df_c['ITM_E'], 
                                df_c['ITM_N'], 
                                #s = prop_val * plot_c_sar['variable'] * 15 / vmax_val, 
                                c = df_c['values'], 
                                vmin = 0, 
                                vmax = 15,#df_c['values'].max(), 
                                cmap = cmap,
                                edgecolors= 'k',
                                lw = 0.1,
                                alpha = 0.75,
                                zorder = 15
                                ) 

            # # Highlight largest values
            #plot_c = plot_c.sort_values(by = ['variable'], ascending=False)
            #r_alpha = [1.0]
        
            # Select largest flood
            #var_3 = ax1.scatter(plot_c['ITM_E'][0:len(r_alpha)], 
            #                   plot_c['ITM_N'][0:len(r_alpha)], 
            #                   s = plot_c['variable'][0:len(r_alpha)]* 15 / vmax_val, 
            #                   facecolors='none',
             #                  edgecolors= 'lime',
             #                  lw = 1,
             #                  zorder = 16,
               #                alpha = 1
             #                  ) 
            var_3 = 1
            
            return(var_2, var_3)


        
        
    #######################################################################
    def remove_layers_from_plot(list_var):
        for i in list_var:
            try:
                i.remove()
            except:
                pass
        return()       
#     #######################################################################
#     def create_background_bars(df, 
#                                num_sites
#                                ):
    
#         # select top 10 flood area to be considered
#         msv = pd.DataFrame(data = df.columns[0:-1],
#                            columns = ['Sites']
#                            )
        
#         msv['Total_Counts'] = np.zeros(len(msv))
        
#         # Check for maximum values
#         for i in dates:
#             select_sar = pd.DataFrame(data = df[df.index == i].iloc[0,0:-1])
#             select_sar = select_sar.sort_values(by = [i], 
#                                                 ascending=False
#                                                 )
            
#             for j, jval in enumerate (msv['Sites']):
#                 if jval == select_sar.index[0]:
#                     msv.loc[j, 'Total_Counts' ] += 1
        
#         # Organize the data and create the files
#         msv = msv.sort_values(by = ['Total_Counts'], 
#                               ascending = False
#                               )
    
#         sel_msv = msv.iloc[0:num_sites]
#         max_val_floods = sel_msv['Total_Counts'].iloc[0]
    
#         sel_msv = sel_msv.sort_values(by = ['Total_Counts'], 
#                                       ascending = True)
        
#         sel_msv = sel_msv.assign(Total_Counts=0)
#         sel_msv = sel_msv.reset_index(drop=True)
        
#         return(sel_msv, 
#                max_val_floods
#                )    
    
#     #######################################################################
#     def create_images(i):
#         (x_rain, 
#          y_rain, 
#          z_rain
#          ) = interpolation_2d(sites, 
#                               rain, 
#                               rain_cc, 
#                               xmin, 
#                               xmax, 
#                               ymin, 
#                               ymax, 
#                               r_cell_size,
#                               i,
#                               delta_t = 0
#                               )
        
#         # 7.2 Add rain values to figure
#         var_1 = plot_rainfall_data(fig,
#                                     ax1, 
#                                     x_rain, 
#                                     y_rain, 
#                                     z_rain, 
#                                     mask_rain
#                                     )
        
    
#         # 7.3 Add hydrogeological data
#         var_2, var_22, var_3 = add_hydro_data_to_plot(fig, 
#                                                       ax1, 
#                                                       i, 
#                                                       df_sar, 
#                                                       df_sar_o, 
#                                                       df_c_sar, 
#                                                       prop_val, 
#                                                       vmax_val, 
#                                                       vmax_stage, 
#                                                       alpha_v, 
#                                                       fix_stage
#                                                       )
    
#         # 7.4 Add colorbar
#         set_colorbars(fig, 
#                        ax1, 
#                        var_1, 
#                        var_2, 
#                        var_3, 
#                        convert_to, 
#                        fix_stage, 
#                        vmax_val, 
#                        prop_val
#                        )
    
#         # 7.5 Add title to figure
#         add_title(title1, i)
    
#         # 7.6 Add stats to plot
#         add_stats_bars(ax2, df_sar, df_sar_sel, i)
        
        
#         # 7.6 Save figure
#         kdate = def_name_based_reference_date(i)
#         #print(kdate)
#         plt.savefig(tmp_path
#                     + str('%04d' % kdate) 
#                     + '.png', 
#                     transparent = True, 
#                     format = 'png',
#                     dpi = 300
#                     )
    
#         # 7.6 Remove hydro data from figure
#         remove_layers_from_plot([var_1, var_2, var_3, var_22])
    
#         return(kdate)


#     #######################################################################
#     def image_larges_values(path1, 
#                             title, 
#                             data_type, 
#                             path2
#                             ):
        
#         msv = pd.read_csv(path1 
#                           + 'stats_' 
#                           + str(title) 
#                           + '_'
#                           + str(data_type) 
#                           +'.csv'
#                           )
        
#         pd_legend = msv.copy()
        
#         # Plot chart
#         plt.figure()
#         N = 12
#         ind = np.arange(N)    # the x locations for the groups
#         accum = msv.iloc[0,2::].copy() * 0
#         for i in range(0,len(msv)):
#             try:
#                 if msv.iloc[i,2::].max() > 0:
#                     #if i != 0:
#                     plt.bar(ind, 
#                             msv.iloc[i,2::], 
#                             #width,
#                             bottom=accum
#                             )
                    
#                     accum = accum + msv.iloc[i,2::]
                        
#                     plt.ylabel('frequency')
#                     plt.title('larger floods per month')
#                     plt.xticks(ind, ('Jan', 
#                                      'Feb', 
#                                      'Mar', 
#                                      'Apr', 
#                                      'May', 
#                                      'Jun', 
#                                      'Jul', 
#                                      'Aug', 
#                                      'Sep', 
#                                      'Oct', 
#                                      'Nov', 
#                                      'Dec'
#                                      ))
                    
#                 else:
#                     pd_legend = pd_legend[~pd_legend["Sites"].isin(msv['Sites'].iloc[i])]
#             except:
#                 pass
        
#         max_val = msv.iloc[:,2::].sum(axis = 0)
#         #print(max_val)
#         v_max = np.nanmax(max_val)
#         plt.yticks(np.arange(0, v_max + 15, 10))    
#         plt.legend(pd_legend['Sites'].iloc[:])
        
#         mng = plt.get_current_fig_manager()
#         mng.window.showMaximized()
#         plt.savefig(path2 
#                     + 'stats_'
#                     + str(title)
#                     + '_'
#                     + str(data_type) 
#                     +'.png', 
#                     transparent = True, 
#                     format = 'png',
#                     dpi = 300
#                     )
    
#         return()

# #######################################################################
# #######################################################################
# class animation():
#     ###################################################################
#     def convert_frames_to_video(pathIn,pathOut,fps):
#         frame_array = []
#         files = [f 
#                  for f in sorted(os.listdir(pathIn)) 
#                  if isfile(join(pathIn, 
#                                 f
#                                 ))]
                 
#         for i in range(0,len(files)):
#             print(i)
#             if '.png' in files[i]:
#                 filename=pathIn + files[i]
#                 img = cv2.imread(filename)
#                 height, width, layers = img.shape
#                 size = (width,height)
                
#                 frame_array.append(img)
        
#         for i in range(0,len(fps)):
#             out = cv2.VideoWriter(pathOut[i],
#                                   cv2.VideoWriter_fourcc(*'XVID'),  
#                                   fps[i], 
#                                   size)
        
#             for j in range(0,len(frame_array)):
#                 out.write(frame_array[j])
#             out.release()
    
#         return()
    
#     #######################################################################
#     def collect_frames(pathIn):
            
#         frame_array = []
#         files = [f 
#                  for f in sorted(os.listdir(pathIn)) 
#                  if isfile(join(pathIn, 
#                                 f
#                                 ))]
#         #print(files)
#         for i in range(0,len(files),2):
#             if '.png' in files[i]:
#                 filename=pathIn + files[i]
#                 img = cv2.imread(filename)
#                 height, width, layers = img.shape
#                 size = (width,height)
                
#                 frame_array.append(img)
        
#         return(size, 
#                frame_array
#                )
        
#     ########################################################################
#     def save_videos(pathOut, fps, f_array, size):
#         for i in range(0,len(fps)):
            
#             out = cv2.VideoWriter(pathOut[i],
#                                   cv2.VideoWriter_fourcc(*'XVID'),  
#                                   fps[i], 
#                                   size)
        
#             for j in range(len(f_array)):
#                 out.write(f_array[j])
#             out.release()
    
#         return()
      
# #######################################################################
# #######################################################################
# class data_interpretation():
#     ###################################################################
#     def compute_statistics(df_sar, 
#                            tmp_path, 
#                            title1, 
#                            convert_to, 
#                            out_path, 
#                            dates
#                            ):
        
#         print('Computing statistics...')
#         # 9.1 Define dataframe for maximum values
#         create_csv_largest_values_per_month(df_sar, 
#                                             tmp_path, 
#                                             title1, 
#                                             convert_to,
#                                             dates
#                                             )
        
#         print('Plotting statisctics')
#         # 9.2 Create figure with statistics
#         image_larges_values(tmp_path, 
#                             title1, 
#                             convert_to, 
#                             out_path
#                             )
    
#         return()