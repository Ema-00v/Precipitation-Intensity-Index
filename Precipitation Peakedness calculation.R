pacman::p_load(pacman, stars, dplyr, tidyverse, ggplot2, cubelyr, data.table)

#Data input
{
#Load the .nc file
prec <- read_stars("/home/admin/Desktop/Research/ATO4Water/R Index Analysis/NWIOIprecDAY.nc")

#Change the name of the NWIOIprecDAY attribute
attributes(prec)$names <- "Prec"

#Define reference system
crs <- st_crs('EPSG:4326')
prec <- st_set_crs(prec, crs)

#Obtain maximum daily precipitation for each month
prec <- aggregate(prec, by="months", FUN=max)
prec <- aperm(prec, c(2,3,1)) #Change dimensions back to x,y,time

#Time vector
t <- time(prec)

#Example plot
# ggplot()+
#   geom_stars(data=prec[,,,1])
}

#Define function to calculate the index and choose if monthly or yearly index
{
p_load(fitdistrplus, FAdist)
# peak_ind_calc <- function(x){ #Index comparing to monthly maximum precipitations
#   p <- unlist(x)
# 
#   #Check if all values are missing and return NA
#   if(all(is.na(p))){
#     return(rep(NA, length=length(p)))
#   }
# 
#   peak_ind <- rep(NA, length = length(p)) #Prepare output
# 
#   for(i in 1:12){
#     p_month <- p[month(t)==i] #Extract monthly values
# 
#     #Check if enough values are present
#     if(sum(is.na(p_month))>(length(p_month)-2)){
#       peak_ind[month(t)==i] <- NA
#     }else{
#       #Estimate Gumbel parameters
#       scale_est <- (sd(p_month, na.rm = TRUE)*sqrt(6))/pi
#       loc_est <- mean(p_month, na.rm = TRUE) + 0.5772157*scale_est
# 
#       #Adapt a gumbel distribution
#       gumbel_fit <- fitdist(p_month[!is.na(p_month)], dgumbel,
#                             start = list(scale = scale_est, location=loc_est))
#       #Calculate quantiles
#       p_month <- pgumbel(p_month, scale=gumbel_fit$estimate[1],
#                          location=gumbel_fit$estimate[2])
#       #Normalization
#       peak_ind[month(t)==i] <- qnorm(p_month, 0, 1)
#     }
#   }
#   return(peak_ind)
# }

peak_ind_calc <- function(x){ #Index comparing to all maximum precipitations
  p <- unlist(x)

  #Check less than 2 values are present and in case return NA
  if(sum(is.na(p))>length(p)-2){
    return(rep(NA, length=length(p)))
  }

  peak_ind <- rep(NA, length = length(p)) #Prepare output

  #Estimate Gumbel parameters
  scale_est <- (sd(p, na.rm = TRUE)*sqrt(6))/pi
  loc_est <- mean(p, na.rm = TRUE) + 0.5772157*scale_est

  gumbel_fit <- fitdist(p[!is.na(p)], dgumbel,
                        start = list(scale = scale_est, location=loc_est)) #Adapt a gumbel distribution
  p <- pgumbel(p, scale=gumbel_fit$estimate[1],
                     location=gumbel_fit$estimate[2]) #Calculate quantiles
  peak_ind <- qnorm(p, 0, 1) #Normalization
  return(peak_ind)
}
}

#Peakedness index calculation for cells
{
  # peak_ind <- st_apply(prec, c("x","y"),
  #                 function(x) peak_ind_calc(x))
  # 
  # attributes(peak_ind)$names <- "Peak_Ind"
  # peak_ind <- aperm(peak_ind, c(2,3,1))
  # 
  # ggplot()+
  #   geom_stars(data=peak_ind[,,,400])
  # peak_ind_thr <- peak_ind
  # peak_ind_thr[peak_ind_thr <=1.5] <- NA
}

#Load a basins .shp file and aggregate the precipitation data
{
  p_load(terra, tidyterra)
  basins <- vect("")
  basins <- project(basins, "EPSG:4326")
  
  #Sort in decreasing area order
  basins <- basins[order(basins$area, decreasing = TRUE)]
  # basins <- filter(basins, codice=="BANSA")
  
  #Plot to see basins
  # plot(basins)
  
  #Calculate mean maximum precipitation for basins
  prec_basins <- extract(as(prec, "SpatRaster"), basins, mean, exact=TRUE)
  
  #Reformat results as data frame with time series for each basin
  peak_ind_bas <- transpose(prec_basins[,2:dim(prec_basins)[2]])
names(peak_ind_bas) <- basins$codice
}

#Calculate index for basins
peak_ind_bas <- lapply(peak_ind_bas, FUN = peak_ind_calc)

#Plot a map of basins and indices
peak_ind_display <- transpose(peak_ind_bas)
date <- 500
basins[["Peak_ind"]] <- peak_ind_display[[date]]

ggplot(basins)+
  geom_spatvector(aes(fill=Peak_ind))+
  scale_fill_gradient2(
    name = waiver(),
    low = "red",
    mid = "white",
    high = "blue",
    na.value = "grey50",
    aesthetics = "fill",
    limits = c(-3,3))+
    labs(fill = "Peakedness Index",
         title = t[date])
