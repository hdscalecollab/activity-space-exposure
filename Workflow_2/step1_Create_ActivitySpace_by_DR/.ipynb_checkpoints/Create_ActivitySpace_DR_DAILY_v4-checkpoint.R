#### -----------  ### An implementation of the Density Ranking method to create density-based surface
### Method is based on [https://github.com/yenchic/density_ranking]
### 
### {Input}: a directory of .csv files, one file for each participant 
### {Location data format}: (lat, lon) pairs
### {Output}: two raster (.tiff) surface for each participant, 'EPSG:4326' and 'EPSG:4269'
### {log fiel}: the script also creates a log file .txt, recording the time and steps run for each participant
###
### Author: Jiue-An (Jay) Yang, @JiueAnYang
### Last Updated: Dec 20, 2023
### ------------------------------------------------------------------------------------------------------------

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


#### -----------  Part A: Environment and Parameter Setup  -----------------------------------

### Set input .csv file dir, this is the directory with the GPS files
full_dir = "C:/Users/jyang/Documents/TWSA_PA/Codes-R/Data/GPS_PA760/rfh_all_points_for_activitySpace_PA760_10272021/"  # .csv files with GPS file (lon, lat)

### Read all .csv file paths to a list
point_files <- list.files(path=full_dir, pattern="*.csv", full.names=TRUE, recursive=FALSE)

### Set output dir, where the activity space rasters are saved
out_dir  = "C:/Users/jyang/Documents/TWSA_PA/Codes-R/DR_Rasters/50cell_200r_daily_RFH_PA760/"

### Set X,Y range of the study area (e.g. San Diego County)
SD_lon_range = c(-117.596554, -116.080939)
SD_lat_range = c(32.534154,33.505024)

### Set Global Variables
cohort = 'RfH'        # name of the study cohort, for reference
col_lon = 'lon'       # column name of the GPS longitude 
col_lat = 'lat'       # column name of the GPS latitude 
col_dateID = 'date_id_int'    # column name of the dateID 
lon_range = COH_lon_range     # use the X,Y range set above
lat_range = COH_lat_range     # use the X,Y range set above
kernelType = 'quartic'   # define the type of kernel to used for KDE/DR calculation  ('Gaussian', 'quartic'...)
target_radius = 200      # define the radius of kernel to search for points 
target_radius_txt = paste0(target_radius, 'm')   # text for file outputs
cell_size = 50        # define the output cell sie of the activity space rasters  
data_seg = 'all'      # use this for reference if working on full dataset point ('all') or subset points (e.g. 'walking','inVehicle') 


#### -----------  Part B: Define the Smoothing parameter ('h')  -----------------------------------

## The selection of smoothing parameter is corresponding to the intended targeted kernel radius, 
## for example, when using Gaussian Kernal, the empirical smoothing bandwidth ('h') and radius are:
# h0 = 0.001   ### this is the value set in the author's example 
# h0 = 0.0005  ### for 100m
# h0 = 0.0007  ### for 200m
# h0 = 0.0009  ### for 400m
# h0 = 0.001  ### for 400m

### Gaussian Kernal : Set search radius
if (kernelType ==  'Gaussian'){
  if (target_radius == 100) {
    h0 = 0.0005  ### for 100m
  }else if (target_radius == 200) {
    h0 = 0.0007  ### for 100m
  } else if (target_radius == 400) {
    h0 = 0.0009  ### for 100m
  }
}

### Quartic Kernal : Set search radius
if (kernelType ==  'quartic'){
  if (target_radius == 100) {
    h0 = 0.0009  ### for 100m
  }else if (target_radius == 200) {
    h0 = 0.0018  ### for 100m
  } else if (target_radius == 400) {
    h0 = 0.0036  ### for 100m
  }
}


#### -----------  Part C: Main Code Starts Here -----------------------------------

### Create a log file to record the process
start_time_total <- Sys.time()
sink(paste0('LOG_DR-', cohort, '_', data_seg, '_', format(start_time_total, "%Y%m%d_%H%M%S"), ".txt"), append=TRUE, split=TRUE)
print (paste("Task started at : ", start_time_total, sep=""))
print ("--------------------------------------------------------")

### Create a dataframe to store error PT/Dates
df_err <- data.frame(identifier=character(),
                     date_ID=integer(),
                     stringsAsFactors=FALSE)


#### -----------  Part D: Calculate DR Surface  -----------------------------------
### Loop through the participant files and create DR surface for each
for(i in 1:length(point_files)){ 
  
  ### Uncomment this line if you wish to record execution time for each PT   
  start_time_loop <- Sys.time()
  
  ### Get participant ID (PT_ID)
  # File Path : point_files[i]
  # File Name : basename(point_files[i])
  # File Name without Extension: tools::file_path_sans_ext(basename(point_files[i]))
  pt_ID = tools::file_path_sans_ext(basename(point_files[i]))
  pt_ID = unlist(strsplit(pt_ID, split = "_"))[1]
  print(paste("Working on: [", pt_ID, "]", " (", i,"/", length(point_files), ")", sep=""))
  
  ### Read GPS point file
  D0 =read.csv(point_files[i])
  
  ### Select only lon, lat columns
  D1 =  dplyr::select(D0, all_of(col_dateID), all_of(col_lon),all_of(col_lat))
  
  ### Remove point with no GPS info
  D1 = dplyr::filter(D1, .data[[col_lon]] != -180)
  
  ### Keep point within COH Catchment Area
  D1 = dplyr::filter(D1,
                     (.data[[col_lon]] >= COH_lon_range[1]) & (.data[[col_lon]] <= COH_lon_range[2]) &
                       (.data[[col_lat]] >= COH_lat_range[1]) & (.data[[col_lat]] <= COH_lat_range[2]))
  
  #### ------------ DAILY-Processing Starts Here -------------
  ### get a list of unique dates to be processed
  unique_dates_forPT <- as.list(unique(D1[[col_dateID]]))
  
  ### loop through each date and extract subset for that day
  for(i in 1:length(unique_dates_forPT)){  
      
    DATE_ID <- unique_dates_forPT[c(i)]
    
    print(paste("Working on: [", pt_ID, '-', DATE_ID ,"]", " (", i,"/", length(unique_dates_forPT), ")", sep=""))
    D_Day <- dplyr::filter(D1, .data[[col_dateID]] == DATE_ID)

    ### Select only LAT/LON columns
    D_Day = D_Day %>% dplyr::select(all_of(col_lon),all_of(col_lat))
    
    ### Only Run DR when at least one point in subset data for the day
    if(nrow(D_Day)> 0){

      ## Get Min/Max range for Lat/Lon on the day
      lon_range = c( min(D_Day[[col_lon]] - 0.01 ), max(D_Day[[col_lon]]) + 0.01 )
      lat_range = c( min(D_Day[[col_lat]] - 0.009), max(D_Day[[col_lat]]) + 0.009)
      dist_x = distm(c(lon_range[1], lat_range[1]), c(lon_range[2], lat_range[1]), fun = distHaversine)
      dist_y = distm(c(lon_range[1], lat_range[1]), c(lon_range[1], lat_range[2]), fun = distHaversine)
      mean_xy_dist = ( (dist_x / cell_size) + (dist_y / cell_size) ) /2 
      
      ## Set the grid resolution
      ## The resolution of grids. 201 means a 201 by 201 matrix over the range of xlim and ylim.
      ## n_res0 = 4000 ### 242m X 180m
      n_res0 = as.integer(mean_xy_dist)

      
      #### ------------ Adding some error handling steps with 'tryCatch' ------------- 
      
      ### Step 1:  Run Density Ranking within SD County
      D_DR = tryCatch({
        DR(data=D_Day, kernel= kernelType, h=h0, n_res=n_res0, xlim=lon_range, ylim=lat_range)
      }, warning=function(w){
        print(paste0("!!!!! Warning for :", pt_ID, "-" ,DATE_ID))
        print(w)
        return(NA)
      }, error = function(e){
        print(paste0("xxxxx Error for :", pt_ID, "-" ,DATE_ID))
        print(e)
        return(NA)
      })
      
      ### Step 2: Create Raster File if DR run is godd.
      if (!is.list(D_DR)){
        
        ## Add error to output file
        df_err[nrow(df_err) + 1,] <- c(pt_ID, DATE_ID)
        
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
        r <- raster(nrow=len_y, ncol=len_x, xmn=x_min, xmx=x_max, ymn=y_min, ymx=y_max)
        set.seed(1)
        
          
        #### ----- Three Options for Raster Normalization ----- ####
        
        ### Option 1: No Normalization
        # values(r) <- D_DR$gr_alpha
        # r[r == 0] <- NA
        
        ### Option 2: Normalized output value to 0-1
        # raster_value <- normalize(D_DR$gr_alpha, method = "range", range = c(0, 1), margin = 1L, on.constant = "quiet")
        # values(r) <- raster_value
        # r[r == 0] <- NA
        
        ### Option 3: No Normalization and Don't fill NA at 0
        values(r) <- D_DR$gr_alpha
        
          
        #### ----- Additional Steps for Raster post-processing ----- ####
          
        ### Important Step!! Need to flip the Raster on y-axis
        r <- flip(r, direction='y')
        
        ### Convert to Raster and write output to GeoTiff
        xyz <- rasterToPoints(r)
        rst <- rasterFromXYZ(xyz)
        
        #### ----- Write activity space raster to output raster file ----- ####
        
        ### Output option 1: with 4326 projection
        # crs(rst) <- CRS('+init=EPSG:4326')
        # raster_type_name = paste("_norm", target_radius, n_res0, h0, DATE_ID,  sep="_")
        # rater_output_filename <- file.path(out_dir, paste(pt_ID, raster_type_name, ".tif", sep=""))
        # writeRaster(rst, filename=rater_output_filename, overwrite=TRUE, options=c('TFW=YES'))
        
        ### Output option 2: with 4269 projection
        crs(rst) <- CRS('+init=EPSG:4269')
        raster_type_name = paste('_DR', target_radius_txt, cell_size, DATE_ID, sep="_")
        rater_output_filename <- file.path(out_dir, paste(pt_ID, raster_type_name, ".tif", sep=""))
        writeRaster(rst, filename=rater_output_filename, overwrite=TRUE, options=c('TFW=YES'))

      }
      
    }else{
      print (paste("No data for: ", pt_ID, "-", DATE_ID))
      print ("--------------------------------------------------------")
    }
  
  }
  #### ------------ DAILY-Processing Ends Here -------------
  
  ### Uncomment the next few lines here to record execution time for each PT
  # end_time_loop <- Sys.time()
  # time_taken_loop <- end_time_loop - start_time_loop
  # print (paste("Time spent in round: ",time_taken_loop, sep=""))
  # print ("--------------------------------------------------------")  
}


### Export the DF recording error PT/Dates to file, so you can check and revisit if needed
write.csv(df_err, paste0(paste("DR_errors_PtDates", kernelType, target_radius_txt, cell_size, data_seg, sep="_"), '.csv'),  row.names = FALSE)


### Record total time spent for processing into the log file
print ("--------------------------------------------------------")
end_time_total <- Sys.time()
time_taken_total <- end_time_total - start_time_total
print (paste0("All completed at: ", end_time_total))
print (paste("Total Time Spent: ",time_taken_total, sep=""))
sink()
