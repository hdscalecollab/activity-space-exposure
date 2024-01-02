### An implementation of the Density Ranking method to create density-based surface
### Method is based on [https://github.com/yenchic/density_ranking]
### 
### {Input}: a directory of .csv files, one file for each participant 
### {Location data format}: (lat, lon) pairs
### {Output}: two raster (.tiff) surface for each participant, 'EPSG:4326' and 'EPSG:4269'
### {log fiel}: the script also creates a log file .txt, recording the time and steps run for each participant
###
### Author: Jiue-An (Jay) Yang
### Last Updated: June 21, 2022
### ------------------------------------------------------------------------------------------------------------

# setwd("D:/000_User_Documents/COH/TWSA_PA/Codes-R/")
setwd("C:/Users/jyang/Documents/TWSA_PA/Codes-R/")
source("DR.R")
library(dplyr)
library(conflicted)
library(stringr)
library(raster)
library(BBmisc)
library(geosphere)
library(progress)

### Workflow of Creating DR surface from Points
### 1) Reading Points from csv into dataframe
### 2) Setting up San Diego Boundary and filtering point with it ( or by catchment area boundary)
### 3) Run Density Ranking
### 4) Create output Rasters 

### ---------------------------------------
### Set input .csv file dir
### ---------------------------------------
# full_dir = "C:/Users/Jay-PC/Documents/UCSD/PQ/Outputs/from_PY/0918_2020/Vehicle/"  # lng, lat
# full_dir = "C:/Users/Jay-PC/Documents/UCSD/PQ/Outputs/from_PY/0105_2021/allPoints/"  # lng, lat
# full_dir = "D:/000_User_Documents/UCSD/PQ/Outputs/from_PY/0608_2022/all_points_withAcc/"  # lng, lat
# full_dir = "D:/000_User_Documents/UCSD/PQ/Outputs/from_PY/0608_2022/testing_TWSA/"  # lng, lat
# full_dir = "D:/000_User_Documents/COH/TWSA_PA/Codes-R/testing_GPS_points/"  # lng, lat
# full_dir = "D:/000_User_Documents/COH/TWSA_PA/Data/all_points_forR/part1/"  # lng, lat
# full_dir = "D:/000_User_Documents/COH/TWSA_PA/Data/all_points_forR/part2/"  # lng, lat
# full_dir = "D:/000_User_Documents/COH/TWSA_PA/Data/GPS_Acc/Part1/"  # lng, lat
# full_dir = "D:/000_User_Documents/COH/TWSA_PA/Data/GPS_Acc/Part2/"  # lng, lat
# full_dir = "C:/Users/jyang/Documents/TWSA_PA/Codes-R/Data/Part1/" # lng, lat
full_dir = "C:/Users/jyang/Documents/TWSA_PA/Codes-R/Data/GPS_Acc_rfh/" # lon,lat


### Read .csv file paths to a list
point_files <- list.files(path=full_dir, pattern="*.csv", full.names=TRUE, recursive=FALSE)  ### Full Run
# point_files <- point_files[c(1)]                                                             ### TEST Run - 1st PT

### ---------------------------------------
### Set output dir
### ---------------------------------------
# out_dir  = "C:/Users/Jay-PC/Documents/UCSD/PQ/Data_for_testing/pts_SDonly_csv_forDR/output/"
# out_dir  = "D:/PQ/DR_output/all_points_01082021/"
# out_dir  = "D:/000_User_Documents/COH/TWSA_PA/Codes-R/all_points_20220621/"  ### Testing
# out_dir  = "D:/000_User_Documents/COH/TWSA_PA/Codes-R/testing/"  ### Testing
# out_dir  = "D:/000_User_Documents/COH/TWSA_PA/Codes-R/DR_Rasters/50cell_100r/"
# out_dir  = "D:/000_User_Documents/COH/TWSA_PA/Codes-R/DR_Rasters/50cell_100r_total/"
# out_dir  = "C:/Users/jyang/Documents/TWSA_PA/Codes-R/DR_Rasters/50cell_400r_total_PQ/"
out_dir  = "C:/Users/jyang/Documents/TWSA_PA/Codes-R/DR_Rasters/50cell_400r_total_RFH/"


### Set X,Y range to San Diego County
# SD_lon_range = c(-117.596554, -116.080939)
# SD_lat_range = c(32.534154,33.505024)
### Set X,Y range to COH Catchment Area
COH_lon_range = c(-119.579046, -114.130793)
COH_lat_range = c(32.534226,35.809260)

### ---------------------------------------
### Set global vars
### ---------------------------------------
stats_type = 'total'
dataset = "rfh"
col_lon = 'lon' #'lng'
col_lat = 'lat'
# col_dateID = 'date_id'
lon_range = COH_lon_range
lat_range = COH_lat_range
kernelType = 'Gaussian'
# kernelType = 'quartic'
target_radius = 400
target_radius_txt = paste0(target_radius, 'm')
cell_size = 50
# data_seg = 'p1'
data_seg = 'all'

### Gaussian Kernal : Set smoothing bandwidth 
# h0 = 0.001
# h0 = 0.0005  ### for 100m
# h0 = 0.0007  ### for 200m
# h0 = 0.0009  ### for 400m
# h0 = 0.001  ### for 400m

if (kernelType ==  'Gaussian'){
  ### Gaussian Kernal : Set search radius
  if (target_radius == 100) {
    h0 = 0.0005  ### for 100m
  }else if (target_radius == 200) {
    h0 = 0.0007  ### for 200m
  } else if (target_radius == 400) {
    h0 = 0.0009  ### for 400m
  }
}

if (kernelType ==  'quartic'){
  ### Quartic Kernal : Set search radius
  if (target_radius == 100) {
    h0 = 0.0009  ### for 100m
  }else if (target_radius == 200) {
    h0 = 0.0018  ### for 200m
  } else if (target_radius == 400) {
    h0 = 0.0036  ### for 400m
  }
}

### Create a log file to record the process
start_time_total <- Sys.time()
sink(paste0('LOG_DR-', dataset ,'_', data_seg, '_', format(start_time_total, "%Y%m%d_%H%M%S"), ".txt"), append=TRUE, split=TRUE)
print (paste("Task started at : ", start_time_total, sep=""))
print ("--------------------------------------------------------")

### Create a DF to store error PT/Dates
df_err <- data.frame(identifier=character(),
                     date_ID=integer(),
                     stringsAsFactors=FALSE)

### Loop through the participant files and create DR surface for each
# for(i in 1:1){    #### TEST Run (first item)
for(i in 1:length(point_files)){    #### FULL Run
  start_time_loop <- Sys.time()
  
  ### Get PT_ID
  # File Path : point_files[i]
  # File Name : basename(point_files[i])
  # File Name without Extension: tools::file_path_sans_ext(basename(point_files[i]))
  pt_ID = tools::file_path_sans_ext(basename(point_files[i]))
  pt_ID = unlist(strsplit(pt_ID, split = "_"))[1]
  print(paste("Working on: [", pt_ID, "]", " (", i,"/", length(point_files), ")", sep=""))
  
  ######### ---- Part 1: Calculate DR Surface ----  #########
  
  ### Read point files
  D0 =read.csv(point_files[i])
  
  ### Select only lon, lat columns
  # D1 =  dplyr::select(D0, all_of(col_dateID), all_of(col_lon),all_of(col_lat))
  D1 =  dplyr::select(D0, all_of(col_lon),all_of(col_lat))
  
  ### Remove point with no GPS info
  D1 = dplyr::filter(D1, .data[[col_lon]] != -180)
  
  ### Keep point within COH Catchment Area
  D1 = dplyr::filter(D1,
                     (.data[[col_lon]] >= COH_lon_range[1]) & (.data[[col_lon]] <= COH_lon_range[2]) &
                       (.data[[col_lat]] >= COH_lat_range[1]) & (.data[[col_lat]] <= COH_lat_range[2]))
  
  ### Only Run DR when at least one point in data
  if(nrow(D1)> 0){

    # Get Min/Max range for Lat/Lon on the day
    lon_range = c( min(D1[[col_lon]] - 0.01 ), max(D1[[col_lon]]) + 0.01 )
    lat_range = c( min(D1[[col_lat]] - 0.009), max(D1[[col_lat]]) + 0.009)
    dist_x = distm(c(lon_range[1], lat_range[1]), c(lon_range[2], lat_range[1]), fun = distHaversine)
    dist_y = distm(c(lon_range[1], lat_range[1]), c(lon_range[1], lat_range[2]), fun = distHaversine)
    mean_xy_dist = ( (dist_x / cell_size) + (dist_y / cell_size) ) /2 
        
    # Set the grid resolution
    # The resolution of grids. 201 means a 201 by 201 matrix over the range of xlim and ylim.
    # n_res0 = 4000 ### 242m X 180m
    n_res0 = as.integer(mean_xy_dist)
  
        
    ################  Adding some error handling tricks #################
        
    ##### Step 1:  Run Density Ranking within SD County
    D_DR = tryCatch({
      DR(data=D1, kernel= kernelType, h=h0, n_res=n_res0, xlim=lon_range, ylim=lat_range)
    }, warning=function(w){
      print(paste0("!!!!! Warning for :", pt_ID))
      print(w)
      return(NA)
    }, error = function(e){
      print(paste0("xxxxx Error for :", pt_ID))
      print(e)
      return(NA)
    })
        
    #### Step 2: Create Raster File if DR run is good.
    if (!is.list(D_DR)){
      
      ## Add error to output file
      df_err[nrow(df_err) + 1,] <- c(pt_ID)
          
    } else{
      
      ## Set row and col of the raster
      len_x <- length(D_DR$x_grid)
      len_y <- length(D_DR$y_grid)
      
      ## Set x,y limits
      x_min = min(D_DR$x_grid)
      x_max = max(D_DR$x_grid)
      y_min = min(D_DR$y_grid)
      y_max = max(D_DR$y_grid)
          
      ## Initiate the raster
      # print ("3: Creating Raster output")
      r <- raster(nrow=len_y, ncol=len_x, xmn=x_min, xmx=x_max, ymn=y_min, ymx=y_max)
      set.seed(1)
          
      ##### ----- NORMALIZATION (start) ----- #####
      
      ## Option 1: No Normalization
      # values(r) <- D_DR$gr_alpha
      # r[r == 0] <- NA
      
      ## Option 2: Normalized output value to 0-1
      # raster_value <- normalize(D_DR$gr_alpha, method = "range", range = c(0, 1), margin = 1L, on.constant = "quiet")
      # values(r) <- raster_value
      # r[r == 0] <- NA
      
      ## Option 3: No Normalization and Don't fill NA at 0
      values(r) <- D_DR$gr_alpha
      
      ##### ----- NORMALIZATION (end) ----- #####
      
      ## Note!! Need to flip the Raster on y-axis
      r <- flip(r, direction='y')
      
      ## Convert to Raster and write output to GeoTiff
      xyz <- rasterToPoints(r)
      rst <- rasterFromXYZ(xyz)
      
      ##### ----- Write to output raster (start) ----- #####
      
      ## Output 1: with 4326 projection
      # crs(rst) <- CRS('+init=EPSG:4326')
      # # raster_type_name = '_norm_1000'
      # raster_type_name = paste("_norm", target_radius, n_res0, h0, DATE_ID,  sep="_")
      # rater_output_filename <- file.path(out_dir, paste(pt_ID, raster_type_name, ".tif", sep=""))
      # writeRaster(rst, filename=rater_output_filename, overwrite=TRUE, options=c('TFW=YES'))
      # print ("Step 4: 4325 Raster created")
      
      ## Output 2: with 4269 projection
      crs(rst) <- CRS('+init=EPSG:4269')
      # raster_type_name = '_conic_1000'
      # raster_type_name = paste("_conic", target_radius, n_res0, h0, DATE_ID, sep="_")
      raster_type_name = paste('_DR', target_radius_txt, cell_size, sep="_")
      rater_output_filename <- file.path(out_dir, paste(pt_ID, raster_type_name, ".tif", sep=""))
      writeRaster(rst, filename=rater_output_filename, overwrite=TRUE, options=c('TFW=YES'))
      # print ("Step 5: 4269 Raster created")
      
      ##### ----- Write to output raster (end) ----- #####
    }
      
  }else{
    print (paste("No data for: ", pt_ID))
    print ("--------------------------------------------------------")
  }R
  
}
##### PT-Processing Ends Here
### Record execution time for each round
# end_time_loop <- Sys.time()
# time_taken_loop <- end_time_loop - start_time_loop
# print (paste("Time spent in round: ",time_taken_loop, sep=""))
# print ("--------------------------------------------------------")  


### Export the DF recording error PT/Dates
write.csv(df_err, paste0(paste("DR_errors_PtDates", dataset ,kernelType, target_radius_txt, cell_size, data_seg, sep="_"), '.csv'),  row.names = FALSE)

### Record total time spent in log file
print ("--------------------------------------------------------")
end_time_total <- Sys.time()
time_taken_total <- end_time_total - start_time_total
print (paste0("All completed at: ", end_time_total))
print (paste("Total Time Spent: ",time_taken_total, sep=""))
sink()
