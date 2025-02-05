#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 16:49:20 2025

@author: Jérémy BERNARD, EDYTEM, CNRS
"""
import pandas as pd
from datetime import timedelta, timezone
import geopandas as gpd
from shapely.geometry import Polygon
import numpy as np


def read_meteo(filepath, file_settings, utc_offset, sensor, format_type = "csv",
               var = "T"):
    if format_type == "csv":
        # Read the file without any header first
        df = pd.read_csv(filepath,
                         header = file_settings["csv"]["header"],
                         parse_dates = file_settings["csv"]["parse_dates"],
                         skiprows = file_settings["csv"]["skiprows"],
                         index_col = file_settings["csv"]["index_col"])
        
        # For each of the observed variable
        for v in file_settings["csv"]["col2read"]:
            # Attribute the header reconstructing the values with their decimal separator
            df[v] = df[file_settings["csv"]["col2read"][v]].astype(str).apply(".".join, axis=1).astype(float)
        
        # Keep only datetime and variable values
        df = df[file_settings["csv"]["col2read"].keys()]
        
        # Set the time zone of the datetime index
        df = df.tz_localize(timezone(timedelta(hours = utc_offset)))
        
        # Keep only the variable of interest
        if var == "T":
            df = pd.DataFrame(df["T"].rename(f"T{sensor}"))
        
            return df
        else:
            print(f"Dealing with {var} variable is not yet implemented !!!")
            
            return None
        
    else:
        print(f"Reading {format_type} is not yet implemented !!!")
        
        return None
    
def create_grid(s_bbox_corners, cell_size, srid):
    # Convert to local CRS
    bbox_corner_srid = s_bbox_corners.to_crs(srid)
    
    # Create a small extent to have cells starting a bit further than at 
    # the extreme location of the measurement
    x_ll = bbox_corner_srid.x[0] - cell_size / 2
    y_ll = bbox_corner_srid.y[0] - cell_size / 2
    x_ur = bbox_corner_srid.x[1] + cell_size / 2
    y_ur = bbox_corner_srid.y[1] + cell_size / 2
    
    # Calculate the number of cells along x and y axis
    nx = int(np.trunc((x_ur - x_ll) / cell_size) + 1)
    ny = int(np.trunc((y_ur - y_ll) / cell_size) + 1)
    
    gdf = gpd.GeoDataFrame({"geometry" : [Polygon([(x_ll + i * cell_size, y_ll + j * cell_size),
                                         (x_ll + (i + 1) * cell_size, y_ll + j * cell_size), 
                                         (x_ll + (i + 1) * cell_size, y_ll + (j + 1) * cell_size),
                                         (x_ll + i * cell_size, y_ll + (j + 1) * cell_size),
                                         (x_ll + i * cell_size, y_ll + j * cell_size)])
                            for i in range(nx) for j in range(ny)]},
                           crs = srid)
    gdf["ID"] = [i + 1 for i in range(gdf.index.size)]
    
    return gdf