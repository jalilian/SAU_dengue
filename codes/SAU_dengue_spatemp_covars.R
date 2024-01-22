
library("tidyverse")
library("sf")
library("ecmwfr")
library("ncdf4")
library("ncdf4.helpers")

setwd(tempdir())
# =========================================================
# read map data

download.file(
  url="https://geodata.ucdavis.edu/gadm/gadm4.1/shp/gadm41_SAU_shp.zip",
  destfile="SAU.zip"
  )

unzip("SAU.zip", exdir="SAU/")

SAU_map <- read_sf("SAU/",
                  layer="gadm41_SAU_2")

# =========================================================
# Naif's selected regions

regions <- 
  c("Makkah Al Mukarramah", 
    "Al Madinah", 
    "Jiddah", 
    "Jazan", 
    "Al Taif", 
    "Najran")

ggplot(data=SAU_map) + 
  geom_sf() +
  geom_sf(data=SAU_map %>% 
            filter(NAME_2 %in% regions),
          aes(fill=NAME_2)) +
  scale_fill_discrete(name="Regions") +
  theme_light()

dev.copy2pdf(file="~/Desktop/SAU_dengue.pdf", width=8.5, height=5.5)

# =========================================================
# spatiotemporal covariates
# =========================================================
# European Centre for Medium-Range Weather Forecasts (ECMWF) 
# Copernicus's Climate Data Store (CDS)
# https://cds.climate.copernicus.eu
# DOI: 10.24381/cds.e2161bac

# user credentials for ECMWF data access
user <- "****************"
cds.key <- "********************************"

# set secret ECMWF token
wf_set_key(user=user, key=cds.key, service="cds")

land_data_fun <- function(year, map, map_var, 
                          area=NULL, temp_dir=NULL)
{
  # create a temporary directory to extract the downloaded file
  if (is.null(temp_dir))
  {
    temp_dir <- tempdir()
  }
  # set the working directory to the temporary directory
  setwd(temp_dir)
  
  # geographical region
  if (is.null(area))
  {
    area <- st_bbox(map)
    area <- c(
      # North
      ceiling(10 * area["ymax"]) / 10,
      # West
      floor(10 * area["xmin"]) / 10,
      # South
      floor(10 * area["ymin"]) / 10,
      # East
      ceiling(10 * area["xmax"]) / 10
    )
    area <- as.numeric(area)
  }
  
  # request for getting land data
  request <- list(
    # dataset name
    dataset_short_name = "reanalysis-era5-land",
    # climate variables 
    variable = c("10m_u_component_of_wind", 
                 "10m_v_component_of_wind", 
                 "leaf_area_index_high_vegetation", 
                 "leaf_area_index_low_vegetation", 
                 "skin_temperature", 
                 "surface_pressure",
                 "total_precipitation", 
                 "volumetric_soil_water_layer_1"),
    # temporal framework: year, month, day, hour
    year = as.character(year),
    month = sprintf("%02d", 1:12),
    day = sprintf("%02d", 1:31),
    time = "12:00",
    # geographical region
    #      North, West, South, East
    area = area,
    # output file format
    format = "netcdf.zip",
    # output file name
    target = paste0("landvars_hourly_", year, ".zip")
  )
  
  # check the validity of a data request and login credentials
  wf_check_request(user=user, request=request)
  
  # download the data request
  dpath <- wf_request(user=user, 
                      request=request,
                      transfer=TRUE, 
                      path=getwd(),
                      # waiting time for download to start
                      time_out=3 * 60 * 60,
                      verbose=TRUE)
  print("download completed!")
  
  # extract downloaded Zip file
  unzip(zipfile=dpath, 
        exdir=paste0(temp_dir, "/", year, "/"), 
        overwrite=TRUE)
  
  # open netCDF file containing the data
  nc_data <- nc_open(paste0(temp_dir, "/", year, "/data.nc"))
  
  # extract longitude
  lon <- ncvar_get(nc_data, "longitude")
  # extract latitude
  lat <- ncvar_get(nc_data, "latitude")
  # extract date and time
  dt <- nc.get.time.series(nc_data)
  # list of names of data variables
  vars <- nc.get.variable.list(nc_data)
  
  # create a spatial grid using longitude and latitude
  grid <- 
    expand_grid(lon_idx=1:length(lon),
                lat_idx=1:length(lat)) %>%
    rowwise() %>% 
    mutate(points=list(st_point(c(lon[lon_idx], 
                                  lat[lat_idx])))) %>%
    st_as_sf()
  
  # set the coordinate reference system for the grid
  st_crs(grid) <- st_crs(map)
  
  # conduct a spatial join to determine points inside specific map regions
  grid <- 
    st_join(grid, map, join=st_within) %>%
    filter(!is.na(!!as.name(map_var)))
  
  # convert nc data to an R data.frame
  dat <- vector("list", length=length(vars))
  for (i in 1:length(vars))
  {
    vals <- ncvar_get(nc_data, vars[i])
    if (length(dim(vals)) > 3)
    {
      idx_3 <- 
        apply(vals, MARGIN=3, 
              function(x){ mean(is.na(x)) } 
        )
      idx_3 <- which.min(idx_3)
      vals <- vals[, , idx_3, ]
    }
    
    # define dimension names and indices
    dimnames(vals) <- 
      list(lon_idx=1:length(lon), 
           lat_idx=1:length(lat),
           time_idx=1:length(dt))
    
    # convert data array to a data frame
    dat[[i]] <- as.data.frame.table(vals) %>%
      mutate(lon_idx=as.integer(lon_idx),
             lat_idx=as.integer(lat_idx),
             time_idx=as.integer(time_idx))
    
    # aggregate by map regions
    # determine grid points inside map regions
    dat[[i]] <- dat[[i]] %>%
      left_join(grid, by=join_by(lon_idx, lat_idx)) %>%
      filter(!is.na(!!as.name(map_var)))
    # retrieve data from lon, lat, and time index
    dat[[i]] <- dat[[i]] %>%
      mutate(longitude=lon[lon_idx],
             latitude=lat[lat_idx],
             time=dt[time_idx]) 
    
    # aggregate by the year and week
    #   week is the number of complete 7 day periods that have 
    #   occurred between the date and January 1st, plus one.
    dat[[i]] <- dat[[i]] %>%
      mutate(year=year(time), 
             week=week(time)) %>%
      group_by(!!as.name(map_var), 
               year, week) %>%
      summarise(mean=mean(Freq, na.rm=TRUE),
                min=min(Freq, na.rm=TRUE),
                max=max(Freq, na.rm=TRUE),
                sd=sd(Freq, na.rm=TRUE)) %>%
      ungroup()
    # rename variables for clarity
    dat[[i]] <- dat[[i]] %>% 
      setNames(c(
        colnames(dat[[i]])[1:3],
        paste(vars[i], 
              c("mean", "min", "max", "sd"), sep="_")
      ))
  }
  
  # close the open netCDF file
  nc_close(nc_data)
  
  # combine data frames for different variables using a full join
  dat <- dat %>% 
    reduce(full_join, 
           by = join_by(!!as.name(map_var), 
                        year, week))
  return(dat)
}

# retrieve land covariate data for the years 2013 to 2022
land_covars <- 
  lapply(2019:2023, land_data_fun, 
         map=SAU_map %>% filter(NAME_2 %in% regions), 
         map_var="NAME_2")

# garbage collection to reduce memory usage
gc(verbose=TRUE, full=TRUE)

# combine data frames for different years
land_covars <- land_covars %>% 
  reduce(full_join)

# save the extracted land covariates
saveRDS(land_covars, file="spat_temp_covars.rds")

land_covars %>% write_csv(file="~/SAU_land_covars.csv")
