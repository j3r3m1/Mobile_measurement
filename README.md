# Mobile_measurement

## Before running the script
Before running this script, make sure:
    1 - you have created a new subdirectory in the "Data" folder with a new name
    2 - to put your air temperature data in a subdirectory of the folder created step 1:
        2.1 - call this subdir "Hobo"
        2.2 - put your data as "ID.csv" or "ID.xlsx" where ID should be replaced by the data ID
    3 - to put your GPS data in a subdirectory of the folder created step 1:
        2.1 - call this subdir "Strava"
        2.2 - put your data as "ID.gpx" where ID should be replaced by the data ID
    4. Modify the variables located in section 1 of the following script
    
    
## What does this script does ?
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


## Can be added in the future
1. 2 sensors may have been set on a same pole: one with a shelter and
the second without shelter: code to highlight differences should be develop
2. A file containing RSU wih urban indicators could be analyzed using the same cell of analyse than the meteorological data in order to identify relationships between urban morphology and meteorology
