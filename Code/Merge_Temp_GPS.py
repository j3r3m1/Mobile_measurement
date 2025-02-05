#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 11:15:12 2025

@author: Jérémy BERNARD, EDYTEM, CNRS

This script:
    1. Gathers data from air temperature sensors (Hobo) and GPS data (Strava)
    into a GIS file,
    2. Averages the difference of air temperature between each cell of a 
    grid and one (or several) reference cells of the same grid.
    If no grid is given as input, a regular square grid is created. 
    To account for shelter + sensor thermal inertia, two methods can be used
    for averaging:
        - "lag": data are considered in a cell only after a given duration in the grid 
        (it can be 0 for testing the sensitivity of this method),
        - "moving_avg": a moving average (along time) is used but the value of the window
        is attributed to the beginning of the window

Before running this script, make sure:
    1 - you have created a new subdirectory in the "Data" folder with a new name
    2 - to put your air temperature data in a subdirectory of the folder created step 1:
        2.1 - call this subdir "Hobo"
        2.2 - put your data as "ID.csv" or "ID.xlsx" where ID should be replaced by the data ID
    3 - to put your GPS data in a subdirectory of the folder created step 1:
        2.1 - call this subdir "Strava"
        2.2 - put your data as "ID.gpx" where ID should be replaced by the data ID
    4. Modify the variables located in section 1 of the following script

Note that 2 sensors may have been set on a same pole: one with a shelter and
the second without shelter.
    
""" 
from pathlib import Path
import os
import pandas as pd
import geopandas as gpd
from glob import glob
import functions as fct
import datetime
import numpy as np
from shapely.geometry import Point
from datetime import timedelta, timezone

###############################################################################
############### 1. DECLARE VARIABLES ##########################################
###############################################################################
# Name of the dataset (directory) of interest
# Possible solutions
#   - Training
dataset_dir_name = "Training"

# Variable of interest (current possible solutions: "T")
var = "T"

# Select steps that need to be run
step1 = False        # Gather GPS and meteorological variables
step2 = True        # Calculates average values within grid cells

# Reference sensor(s) - it is a list since it can be the mean of several sensors
ref_sens = [23, 25]

# Start and end time of the measurement campaign in UTC
campaign_start = datetime.datetime(2023,6,12,21,14,3, tzinfo = datetime.timezone.utc)
campaign_end = datetime.datetime(2023,6,12,21,33,39, tzinfo = datetime.timezone.utc)

# Path of the file to use as grid / else size of the grid used to average the data (m)
grid_info = 10
grid_srid = 2154

# Averaging method chosen ("lag" or "moving_avg")
avg_meth = "moving_avg"

# Time spent in a grid cell before taking into account the data / 
# windows size for moving average (s)
time_lag = 30

# Hobo file format the reading function differ depending the format
meteo_fileformat = "csv"

# Percentiles to use for min and max values per grid cell
pval_max = 0.9
pval_min = 0.1

###############################################################################
############### 2. ALREADY SET VARIABLES ######################################
###############################################################################

# AIR TEMPERATURE INFORMATIONS
# Directory name
meteo_dir_name = "Hobo"
meteo_utc_offset = 0
meteo_file_settings = {"csv" :  {"header": None,
                                 "parse_dates" : True,
                                 "skiprows" : [0],
                                 "index_col" : 1,
                                 "col2read" : {"T": [2,3],
                                               "RH": [4,5],
                                               "Dew": [6,7]}},
                       "xlsx" : ""}


# GPS INFORMATIONS
# GPS directory name
gps_dir_name = "Strava"

# OUTPUT INFORMATIONS
output_dir_name = "Output"

###############################################################################
############### 3. DECLARE PATHS ##############################################
###############################################################################
data_dir = os.path.join(Path(os.path.abspath(os.path.curdir)).parent, "Data")
meteo_dir = os.path.join(data_dir, dataset_dir_name, meteo_dir_name)
gps_dir = os.path.join(data_dir, dataset_dir_name, gps_dir_name)
output_dir = os.path.join(data_dir, dataset_dir_name, output_dir_name)



###############################################################################
############### 4. SCRIPT BEGINS #################################################
###############################################################################
# Get lists of GPS and meteorological data filepaths
meteo_filepaths = glob(meteo_dir + os.sep + "*." + meteo_fileformat)
gps_filepaths = glob(gps_dir + os.sep + "*.gpx")

# Convert list to pandas series with as index the sensor ID
s_meteo_filepaths = pd.Series({Path(f).stem : f 
                               for f in meteo_filepaths})
s_gps_filepaths = pd.Series({Path(f).stem : f 
                               for f in gps_filepaths})

# Identify sensors that are both in GPS and meteo (get rid of others)
list_sensors = s_meteo_filepaths.index.intersection(s_gps_filepaths.index)


#------------------------------------------------------------------
# 1. GATHER GPS AND METEOROLOGICAL DATA
#------------------------------------------------------------------
if step1: 
    print("Step 1 started") 
    if os.path.exists(output_dir):
        print(f"Output directory '{output_dir}' already exists. Please remove it before run")
    else:
        os.mkdir(output_dir)
        dic_merge = {}
        # Read meteo data
        df_all_meteo = pd.concat({f"{var}{s}" : \
                                     fct.read_meteo(filepath = s_meteo_filepaths[s],
                                                    file_settings = meteo_file_settings,
                                                    utc_offset = meteo_utc_offset,
                                                    sensor = s,
                                                    format_type = meteo_fileformat,
                                                    var = var)
                                     for s in list_sensors}.values(), 
                                 axis = 1)
        # For each sensor in the final list
        for s in list_sensors:
            # Reference meteorological condition
            df_meteo_ref = df_all_meteo[[f"{var}{s_i}" for s_i in ref_sens]].mean(axis = 1)
            
            # Calculate relative micrometeorological variable
            df_meteo = df_all_meteo[f"{var}{s}"].subtract(df_meteo_ref).rename(var)
            
            # Read GPS data
            df_gps = gpd.read_file(s_gps_filepaths[s], 
                                   layer = "track_points")[["geometry", "time"]]\
                            .drop_duplicates(["time"]).set_index("time")
            
            # Merge GPS and meteo data
            df_merge = pd.concat([df_meteo, df_gps], sort = True, axis = 1)
            
            # Temporarily convert geometries to string to better deal with the data
            df_merge["geometry"] = df_merge["geometry"].astype(str).replace("None", np.nan)
            
            # Filters data where there is no more GPS data to avoid 
            # wrong interpolation in the following
            df_merge = df_merge[df_merge.index <= df_merge["geometry"].dropna().index[-1]]
            
            # Filter only data belonging to the measurement campaign
            df_merge = df_merge[(df_merge.index >= campaign_start) \
                                * (df_merge.index <= campaign_end)]
            
            # Fill GPS values with previous positions when there is no data
            df_merge["geometry"] =  gpd.GeoSeries.from_wkt(df_merge["geometry"].ffill(),
                                                           crs = 4326)
            
            # Save the results in an output folder
            gpd.GeoDataFrame(df_merge.reset_index(names = "time")).to_file(os.path.join(output_dir, f"{s}.fgb"))
            
#------------------------------------------------------------------
# 2. AVERAGE THE VARIABLES WITHIN GRID CELL
#------------------------------------------------------------------
if step2:
    print("Step 2 started")
    # If a grid file is given as input
    if type(grid_info) == str:
        if os.path.exists(grid_info):
            gdf_grid = gpd.read_file(grid_info)
            # Set the index as ID
            gdf_grid["ID"] = gdf_grid.index
        
    # Else create a grid file using square grid cells
    else:
        xmin, ymin, xmax, ymax = \
            pd.concat({s: gpd.read_file(os.path.join(output_dir, f"{s}.fgb"))
                             for s in list_sensors}, 
                            ignore_index = True).sort_values(axis = 0,
                                                             by = "time").total_bounds
        s_bbox_corners = gpd.GeoSeries([Point(xmin, ymin), Point(xmax, ymax)],
                                     crs = 4326)
        gdf_grid = fct.create_grid(s_bbox_corners = s_bbox_corners,
                                   cell_size = grid_info,
                                   srid = grid_srid)
        gdf_grid.to_file(os.path.join(output_dir, f"grid.fgb"))
        
    # Identify measurement points intersecting each grid cell
    # Do it for each station independantly to implement rules for the averaging
    join_fin = {}
    for s in list_sensors:
        s_obs = gpd.read_file(os.path.join(output_dir, f"{s}.fgb")).to_crs(grid_srid)
        
        if avg_meth == "lag" or avg_meth == "moving_avg":
            # Identify the intersection between station location and grid cell
            # and sort in ascending datetime
            join = gpd.sjoin(left_df=s_obs[["time","geometry", "T"]],
                             right_df=gdf_grid[["ID","geometry"]], 
                             how="left", 
                             predicate="intersects").sort_values("time")        
            if avg_meth == "lag":
                # Add a column with the time where the station get into a new grid cell
                join["in_time"] = join["time"][join["ID"].diff() != 0].reindex(join.index)
                # Calculate the time spent in the grid since the entrance
                join["in_time"] = join["time"].subtract(pd.Series(join.set_index("time")["in_time"].ffill().values,
                                                                  index = join.index).dt.tz_localize(timezone(timedelta(hours = 0))))
                # Remove points being in a cell since less than 'time_lag'
                join_fin[s] = join[join["in_time"] >= pd.Timedelta(pd.offsets.Second(time_lag))]\
                    .drop(["in_time"], axis = 1)
                
            elif avg_meth == "moving_avg":
                # First move the air temperature further in the future
                join = join.set_index("time")
                join_T = join["T"].shift(periods = 1,
                                         freq = pd.offsets.Second(-time_lag))
                # Then calculate a rolling average, else there is no more geometries for time some values
                join["T"] = join_T.reindex(join_T.index.union(join.index))\
                    .rolling(window = pd.offsets.Second(time_lag), 
                             min_periods = int(time_lag / 4)).mean()\
                        .reindex(join.index)
                join_fin[s] = join.copy(deep = True)
        else:
            print(f"The averaging method '{avg_meth}' does not exists. Please select a valid one.")
        
    # Gather all measurements in a single dataframe
    grid_fin = pd.concat(join_fin, ignore_index = True)
    
    # Aggregate values by grid cell using some statistics
    gdf_grid_fin = gdf_grid.set_index("ID")
    gdf_grid_fin[f"{var}_mean"] = grid_fin.groupby("ID")[var].mean()
    gdf_grid_fin[f"{var}_q{pval_max}"] = grid_fin.groupby("ID")[var].quantile(pval_max)
    gdf_grid_fin[f"{var}_q{pval_min}"] = grid_fin.groupby("ID")[var].quantile(pval_min)
    gdf_grid_fin[f"{var}_median"] = grid_fin.groupby("ID")[var].median()
    gdf_grid_fin[f"{var}_sqrt"] = grid_fin.groupby("ID")[var].std()
    
    # Save the resulting file
    gdf_grid_fin.to_file(os.path.join(output_dir, f"grid_fin.fgb"))