rm(list=ls()) 

### -- Package Required -- ###
suppressMessages({
  library("ggplot2")
  library("parallel")
  library("fields")
  library("mgcv")
  library("maps")
  library("dplyr")
  library("plyr")
  library("ggpubr")
  library("viridis") # viridis colour scale
  library("reshape2")
  
  library("colorspace")
  library("grDevices")
  
  ## Load ismev library for marginal GPD fits:
  library("ismev")
  library("evd")
  
  
  library("tictoc")
  library("extRemes")
  
  library("scales") # plotting against time
  
  # Removed these dependencies as they were causing problems. They are just
  # needed to make the plots pretty (e.g., showing degrees on the axis).
  #library("ggspatial")
  #library("ggOceanMapsData")
  #library("ggOceanMaps")
})
########################

# -- Exploratory Analysis -- :
# load("~/Downloads/DATA.Rdata") 
load("~/Downloads/redseatemperature.rdata")
print(dim(loc))
# dimension: 16703     2

loc_full <- loc

## Extract observations from the Summer (July, August, September):
summer_months <- month %in% c(7,8,9)
data  <- data[summer_months, ]
month <- month[summer_months]
year  <- year[summer_months]
time  <- time[summer_months]

## Extract observations from the southern part of interest:
idx  <- loc[, "lat"] < 20.2 & loc[, "lat"] > 15.75
data <- data[, idx]
loc  <- loc[idx, ]

## Transform locations:
loc[, "lon"] <- loc[, "lon"] * 1.04
loc[, "lat"] <- loc[, "lat"] * 1.11

loc_full[, "lon"] <- loc_full[, "lon"] * 1.04
loc_full[, "lat"] <- loc_full[, "lat"] * 1.11


##  -- Plot the full data set -----------------------------: 
# Create a ggplot with separate colors for loc_full and loc
loc_full_dataframe <- data.frame( lon = loc_full[, "lon"] ,
                                  lat = loc_full[, "lat"] ,
                                  group = rep("loc_full", length(loc_full[, "lat"]))
)

loc_dataframe <- data.frame( lon = loc[, "lon"] ,
                             lat = loc[, "lat"] ,
                             group = rep("loc", length(loc[, "lat"]))
)

figure <- ggplot() +
  geom_point(data= loc_full_dataframe, aes(x = lon, y = lat, color = group), size = 0.3, ) +
  geom_point(data = loc_dataframe, aes(x = lon, y = lat, color = group), size = 0.3) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() + coord_fixed()
print(figure)

# Save the ggplot as an image file (e.g., PNG)
ggsave("spatialDomain.png",
       plot = figure,
       width = 6,
       height = 6,
       units = "in",
       path= "~/Downloads/img/RedSea",
       dpi = 300)

## Southern Region: 
figure <- ggplot() +
          geom_point(aes(x = loc[, "lon"], y = loc[, "lat"]), size = 0.3) +
          labs(x = "Longitude", y = "Latitude") +
          theme_bw() + coord_fixed()
ggsave("SouthernRegion.png",
       plot = figure,
       width = 6,
       height = 4,
       units = "in",
       path = "~/Downloads/img/RedSea",
       dpi = 300)

## Full spatial domain: 
figure <- ggplot() +
  geom_point(aes(x = loc_full[, "lon"], y = loc_full[, "lat"]), size = 0.3) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() + coord_fixed()
ggsave("FullDomain.png",
       plot = figure,
       width = 6,
       height = 4,
       units = "in",
       path = "~/Downloads/img/RedSea",
       dpi = 300)
### -------------------------------------

# Rename the full data so that we can access it in each of the following sections:
full_data <- data # 2852 6239
full_loc  <- loc
print(dim(full_loc)) # 6239    2

## Function to transform to GPD margins:
# transforming daily maximum temperature into [-1,0] uniform by apply below 
# X^U = F_i(X_t) - 1 = - max{0, (1 + (xi_hat/sigma_hat) * X_t)}^(-1/xi_hat)
#
GPDMargins <- function(univariateData, thrQ = 0.95)
{
  size <- length(univariateData)
  empiricalMargins <- rank(univariateData)/(size+1)
  
  empiricalMargins[univariateData>quantile(univariateData, thrQ)] <- NA
  
  invisible(capture.output(fit <- gpd.fit(univariateData, threshold=quantile(univariateData,thrQ))$mle))
  scale <- fit[1]
  shape <- fit[2]
  
  empiricalMargins[is.na(empiricalMargins)] <- 1 - pgpd(univariateData[univariateData>quantile(univariateData, thrQ)],
                                                        quantile(univariateData, thrQ),
                                                        scale,
                                                        shape,
                                                        lower.tail=FALSE)*(1-thrQ) 
  
  # empiricalMargins[is.na(empiricalMargins)] <- 1 - pgpd(univariateData[univariateData>quantile(univariateData, thrQ)], 
  #                                                       quantile(univariateData, thrQ),
  #                                                       scale,
  #                                                       shape,
  #                                                       lower.tail=FALSE)*(1-thrQ)
  # 
  # laplaceMargins <- rep(NA, length(empiricalMargins))
  # laplaceMargins[empiricalMargins <= 0.5] <- log(2*empiricalMargins[empiricalMargins<=0.5])
  # laplaceMargins[empiricalMargins > 0.5]  <- -log(2*(1-empiricalMargins[empiricalMargins>0.5]))
  
  gpdMargins <- rep(NA, length(empiricalMargins))
  gpdMargins[empiricalMargins <= thrQ] <- empiricalMargins[empiricalMargins<=thrQ] -1
  gpdMargins[empiricalMargins > thrQ]  <- empiricalMargins[empiricalMargins > thrQ] -1
  
  return(gpdMargins)
}


### Plot:
bubblePlot(full_loc, GPDMargins(full_data[,1]), highlight=FALSE)
bubblePlot(full_loc, full_data[,1], highlight=FALSE)


## Transform to GPD margins:
tic()
data_GPD <- suppressWarnings(apply(full_data, 2, GPDMargins))
toc()
# 47.114 sec elapsed

# Save the distance matrix
D <- fields::rdist(full_loc)
save(file = file.path("~/Downloads/img/RedSea", "D.rda"), D)

## -- Geometric mean --: 
geometricMean <- function(space_obs)
{
 gm <- exp(mean(log(space_obs)))
 return(gm)
}

##  Work with the geometric mean risk function:
risk_t <- apply((data_GPD+1), 1, geometricMean) # risk varies across time

## atleast 5 observation every year is captured
# check for spatial de- clustering and threshold selection
store_quantiles <- rep(NA)
t=0
for(i in 1: 31){
  print(i)
  
  store_quantiles[i] <- quantile(risk_t[(1:92)+t], 0.90)
  u_temp <- store_quantiles[i]
  print(which(risk_t[(1:92)+t]>u_temp))
  t=t+ 92}

u_threshold <- quantile(store_quantiles, 0.40)

t=0
for(i in 1: 31)
{
  print(i)
  print(which(risk_t[(1:92)+t]> u_threshold))
  t=t+ 92
}

## threshold exceedance plot
plot(risk_t,
     pch = 19,
     col = alpha('black', alpha = 0.3)) 
abline(h = u_threshold, col = 'magenta', lty = 2)




## -- r-Pareto processes  Y_r --
Y_r <- (data_GPD[risk_t>u_threshold,]+1)/u_threshold
dim(Y_r)


## To be noted r-Pareto could be written as
# Y_r = R*W
# log Y_r = log R + log W
# log W = log Y_r - log R

R <- risk_t[risk_t>u_threshold]
length(R)

log_R <-  log(R)
log_Y_r <- log(Y_r)

log_W = log_Y_r - log_R

## Use LatticeKrig here to model the log_W processes:
library("LatticeKrig")
fit_spatial <-  LKrig(full_loc,
                      t(log_W),
                      a.wght=5,
                      nlevel=1,
                      # NC=10,
                      nu=0.5, 
                      lambda=0.1,
                      verbose = TRUE) 
# need to define the variogram from this:

summary(fit_spatial)

look_org <- bubblePlot(fit_spatial$x, log_W[1,], highlight = FALSE, main= 'log W')
look<- bubblePlot(fit_spatial$x, fit_spatial$y[,1], highlight = FALSE, main= 'Predict Surface log W')

predictSE <- predictSurfaceSE(fit_spatial, verbose=TRUE)


# tic()                             
# for(i in 1:ncol(data_GPD))
# {
#   print(i)
#   for(t in 1:nrow(data_GPD))
#   {
#     i_th_spatial_t_day <- data_GPD[t,i]
#     
#     threshold_exceed_risk[i, t] <- (i_th_spatial_t_day > risk_t[t])
#     extract_index[[t]] <- 
#   }
# }
# toc()
# 72.603 sec elapsed
# 
# # Checking exceedance per year 
# exceed_per_year <- rep(NA)
# for(t in 1:nrow(data_GPD))
# {
#   exceed_per_year[t] <- sum(threshold_exceed_risk[,t])
# }
# 
# 
# # Subset the data based on the threshold
# subset_data <- your_data[your_data$risk >= threshold_per_year, ]
# u <- risk_t[]
# 
# ## Plot some of the fields
# set.seed(1)
# plot_data(extreme_data_L, S, img_path = img_path, type = "random", geom = "tile")
# plot_data(extreme_data_L, S, img_path = img_path, type = "most_extreme", geom = "tile")
# plot_data(extreme_data_L, S, img_path = img_path, type = "least_extreme", geom = "tile")
# 
# ## Find the extreme data
# extreme_data_L <- find_extreme_fields(data_L, s0_idx, path)
# 
# 
# # -- Unique Lon and Lat -- #
# # of location combination = 16703
# x <- unique(loc[,1]) # Longitudes 
# y <- unique(loc[,2]) # Latitudes 
# 
# nx <- length(x) # number of longitude points = 233
# print(nx) # 233
# 
# ny <- length(y) # number of latitude points = 359
# print(ny) # 359
# 
# # Observe on daily basis
# n.times <- nrow(data) ## number of time points = 11315 
# print(n.times)
# 
# n.sites <- ncol(data) ## number of spatial sites = 16703
# print(n.sites)
# 
# n.years <- length(unique(year)) ## number of years = 31
# print(n.years)
# 
# dim(data) # (11315,16703)
# dim(loc) # (16703, 2)
# 
# # -- Annual Maximum --:
# # row: spatial location
# # columns: daily observation
# data_matrix <- t(data) 
# dim(data_matrix) # (16703, 11315)
# 
# # -- Focus on HOT EXTREMES during the SUMMER SEASIN (JUNE to SEPTEMBER) --:
# # RED SEA data has 31 years of daily terms: 
# # i.e., extract the summer observation for 31 years across the 16703 locations
# summer_index <- 152:273 # June to September
# length(summer_index)
# reduce_data_matrix <- array(NA, dim=c(nrow(data_matrix), 31, length(summer_index)))
# 
# tic()                             
# for(i in 1:nrow(data_matrix))
# {
#   print(i)
#   t=0
#   for(j in 1:31)
#   {
#     if(j%%5==0)
#     {
#       print(j)
#     }
#     reduce_data_matrix[i,j,] <- data_matrix[i,c(t+summer_index)]
#     t=t+365
#   }
# }
# toc()
# # 35.254 sec elapsed
# 
# # Plot TIME SERIES:
# plot(c(reduce_data_matrix[1,31,]), 
#      type= 'l')
# plot(c(reduce_data_matrix[1,1,]), 
#      type= 'l')
# 
# # -- FIT GPD to the marginals --: 
# library(ismev)
# GPD_parameter <- matrix(NA,
#                         nrow= nrow(reduce_data_matrix),
#                         ncol= 2)
# threshold_across_space <- rep(NA)
# 
# tic()
# for(i in 1:nrow(reduce_data_matrix))
# {
#   u <- quantile(c(reduce_data_matrix[i, , ]),
#                 probs=0.90)
#   
#   threshold_across_space[i] <- u
#   fit_marginal <- gpd.fit(c(reduce_data_matrix[i, , ]),
#                           npy= 122,
#                           threshold= u)
#   
#   GPD_parameter[i, ] <- fit_marginal$mle
# }
# toc()
# # 220.524 sec elapsed
# # 3.68 mins
# 
# # BubblePlot of the GPD parameters:
# bubblePlot(loc,
#            GPD_parameter[,1],
#            highlight= FALSE,
#            xlab= 'Lon',
#            ylab= 'Lat')
# title(main= expression(sigma),
#       line= 1)
# 
# bubblePlot(loc,
#            GPD_parameter[,2],
#            highlight= FALSE,
#            xlab= 'Lon',
#            ylab= 'Lat')
# title(main= expression(xi),
#       line= 1)
# 
# 
# XU <- matrix(NA,
#              nrow= nrow(data_matrix),
#              ncol= length(summer_index))
# for(i in 1:nrow(data_matrix))
# {
#   # print idiot numbers
#   if(i%%656==0)
#   {
#     print(i)
#   }
#   
#   # Estimated parameter from the marginal fit:
#   sigma <- GPD_parameter[i, 1]
#   xi <- GPD_parameter[i, 2]
#   
#   # F_i(X_t) = -(1+(xi/sigma)*X_t)_+
#   
#   XP[i, ] <- (abs(1+ (GPD_parameter[i,2]/GPD_parameter[i,1])*z))^(1/GPD_parameter[i,2])
# }
# 
# 
# 
# # # Y_r process:
# # risk_across_time <- matrix(NA,
# #                            nrow=years,
# #                            ncol=days)
# # 
# # for(t in 1:years)
# # { 
# #  for(d in 1:days)
# #  {
# #    risk_across_time[t, d] <- as.numeric(exp(mean(log(matrix_3d[,t,d]))))
# #  }
# # }
# # 
# # 
# # # Pick the observation for which the risk exceeds the thresholds:
# # index_threshold_exceedance <- matrix(NA,
# #                                      nrow= years,
# #                                      ncol=days)
# # for(t in 1:years)
# # { 
# #   for(d in 1:days)
# #   {
# #     index_threshold_exceedance[t, d] <- (risk_across_time[t, d]>u)
# #     
# #   }
# # }
# # # Check: 
# # # apply(index_threshold_exceedance, 1, sum)
# # 
# # 
# # # Looking at the threshold exceedance index for each year
# # threshold_exceedance_index <- list()
# # tic()
# # for(t in 1:n.years)
# # {
# #   threshold_exceedance_index[[t]] <- (which(index_threshold_exceedance[1,]==TRUE))
# # }
# # toc()
# # # 1.784 sec elapsed
# # 
# 
# 
# 
# 
# # We are focusing on the only on the summer data to across the year and space
# # to avoid getting into the seasonal variability.
# #
# # COMMENT: reduce_data_matrix is 3D-array which store observation 
# # only for the summer months across the locations and each years
# 
# 
# # -- EXTRACTING RED SEA SST DATA FOCUSING ONLY ON A SMALL REGION -- #
# #
# # NOTE:
# # define the Lat-Lon ranges for the Southern part only
# subset <- list(lon=c(37.5, 44.5),
#                lat=c(17, 22))
# 
# # extract observations within the specified lat and lon ranges
# lon_ind <- (loc[,1]> subset$lon[1] & loc[,1] <= subset$lon[2])
# lat_ind <- (loc[,2]> subset$lat[1] & loc[,2] <= subset$lat[2]) 
# 
# # - subset index -: 
# subset_index <- intersect(which(lat_ind==TRUE),
#                         which(lat_ind==TRUE))
# lon_subset <-  loc[subset_index, 1]
# lat_subset <-  loc[subset_index, 2]
# 
# # -- Lon-Lat subset --:
# loc_subset_ <- loc[subset_index, ]
# print(dim(loc_subset_)) 
# # dimension of subset location: 6290    2
# 
# # -- Further subsetting --:
# every_third <- seq(1, length(lon_subset_), by=9)
# loc_subset <- loc_subset_[every_third,]
# print(dim(loc_subset)) # 699   2
# 
# # file_path_space <- "~/Desktop/Project2+RedSea/Paper 2-dataset/subset_spatial_loc.csv"
# # write.csv(subset_loc_nr, file = file_path_space, row.names=FALSE)
# 
# 
# # ## -- Image Plot --
# # # png("~/Desktop/Project2+RedSea/grid_plot.png",
# # #     width=640,
# # #     height=440,
# # #     res=300)
# # # # Set plot margins
# # # # Set smaller margins and layout parameters
# # # par(mai = c(1, 1, 1, 1),
# # #     mar = c(2, 2, 2, 2),
# # #     oma = c(0, 0, 0, 0))
# # plot(loc[,1],
# #      loc[,2],
# #      cex= 0.1,
# #      pch=19,
# #      asp=1,
# #      col=alpha('black', alpha=0.2),
# #      xlab='Lon',
# #      ylab='Lat')
# # points(subset_loc[,1],
# #        subset_loc[,2],
# #        cex= 0.1,
# #        col=alpha('gray', alpha=0.2),
# #        xlab='Lon',
# #        ylab='Lat')
# ## Save the plot as a PNG file
# # dev.off()
# 
# 
# # Grid Plot of the subset locations: 
# plot(loc_subset[,1],
#      loc_subset[,2],
#      cex= 0.5,
#      pch= 19,
#      col= alpha('gray', alpha=0.6),
#      xlab= 'Lon',
#      ylab= 'Lat',
#      main= 'Subset Grid',
#      line= 1)
# 
# # -- Southern region data --:
# observations_summer <- reduce_data_matrix[loc_subset_, , ] # 6280   31  122
# dim(observations_summer)
# 
# # Subset: 
# observations_summer <- observations_summer[every_third, , ]
# dim(observations_summer) # 698  31 122
# 
# 
# # Histogram: 
# random_locs <- floor(runif(5, min=1, max=nrow(subset_loc)))
# print(random_locs)
# 
# # Plot
# # threshold for ith spatial location and jth year
# hist(observations_summer[random_locs[1],31,])
# temp_obs <- observations_summer[random_locs[1],31,]
# u <- quantile(temp_obs, probs= 0.90)
# plot(temp_obs, pch=19)
# points(which(temp_obs>u),
#        temp_obs[which(temp_obs>u)],
#        col='red',
#        pch=19)
# 
# hist(observations_summer[random_locs[2],31,])
# temp_obs <- observations_summer[random_locs[2],31,]
# u <- quantile(temp_obs, probs=0.90)
# plot(temp_obs, pch =19)
# points(which(temp_obs>u), temp_obs[which(temp_obs>u)], col='red', pch=19)
# 
# hist(observations_summer[random_locs[3],31,])
# temp_obs <- observations_summer[random_locs[3],31,]
# u <- quantile(temp_obs, probs=0.90)
# plot(temp_obs, pch =19)
# points(which(temp_obs>u), temp_obs[which(temp_obs>u)], col='red', pch=19)
# 
# 
# hist(observations_summer[random_locs[4],31,])
# temp_obs <- observations_summer[random_locs[4],31,]
# u <- quantile(temp_obs, probs=0.90)
# plot(temp_obs, pch =19)
# points(which(temp_obs>u), temp_obs[which(temp_obs>u)], col='red', pch=19)
# 
# 
# hist(observations_summer[random_locs[5],31,])
# temp_obs <- observations_summer[random_locs[5],31,]
# u <- quantile(temp_obs, probs=0.90)
# plot(temp_obs, pch =19)
# points(which(temp_obs>u), temp_obs[which(temp_obs>u)], col='red', pch=19)
# 
# 
# # defining threshold and extracting information from each years
# # along with index and deriving common days
# nr <- nrow(observations_summer)
# common_index <- vector("list", length=nr)
# tic()
# for(i in 1:nr)
# {
#   print(i)
#   for(j in 1:31)
#   {
#     if(j%%5==0)
#     {
#       print(j)
#     }
#     # observation from ith spatial location and jth year
#     temp <- observations_summer[i,j, ]
#     # threshold for ith spatial location and jth year
#     u <- quantile(temp, probs= 0.90)
#     common_index[[i]] <- which(temp>u)
#   }
# }
# toc()
# # 8.59 sec elapsed
# 
# # find the union of indexes across all matrices
# exceed_indices <- unique(unlist(common_index))
# print(exceed_indices)
# 
# 
# extreme_obs <- array(NA, dim=c(nr, 31, length(exceed_indices)))
# # extracting extreme obs:
# tic()
# for(i in 1:nr)
# {
#   print(i)
#   for(j in 1:31)
#   {
#     if(j%%5==0)
#     {
#       print(j)
#     }
#     # observation from ith spatial location and jth year
#     temp <- observations_summer[i,j, ]
#     extreme_obs[i,j,] <- temp[exceed_indices]
#   }
# }
# toc()
# # 2.586 sec elapsed
# 
# 
# # -- Building Spatial Model--:
# # work with year: 31 (other option, 14, 26)
# # spatial extreme observation:
# y <- extreme_obs[,31,]
# dim(y) # 698  25
# 
# # Store the CSV file
# file_path <- "~/Desktop/Project2+RedSea/Paper 2-dataset/Threshold_exceed_year31.csv"
# write.csv(y, file = file_path, row.names=FALSE)
# 
# 
# 
# 
# quilt.plot(subset_loc_nr, extreme_obs[,31,1])
# plot(subset_loc_nr, cex=0.2)
# 
# # nr = #spatial grids on the subset
# set.seed(123)
# nos <- floor(runif(n=5, min=1, max=nr))
# print(nos)
# 
# 
# 
# library(LatticeKrig)
# 
# # -- Fit LatticeKrig Package  --:
# x <- subset_loc_nr
# fit <- LatticeKrig(x, y)
# plot(y, predict(fit, x), col=alpha('black',alpha=0.2))
# 
# # -- Histogram --: 
# # randomly selected spatial locations
# hist(y[nos[1],], main='')
# hist(y[nos[2],], main='')
# hist(y[nos[3],], main='')
# hist(y[nos[4],], main='')
# hist(y[nos[5],], main='')
# 
# # Likelihood function: 
# # SAR model: 
# # install.packages("~/Downloads/LatticeKrig_9.0.tar.gz",
# #                  repos=NULL,
# #                  type="source")
# 
# library(LatticeKrig)
# 
# # TRYING out LK on the extreme spatial dataset
# #look<- LKrigLatticeCenters(LKinfo, Level=1) 
# # look <- make.surface.grid( look)
# # plot( look,  cex=.5) #alpha=c(.2) )
# # image.plot(B)
# 
# # -- Log likelihood function --: 
# 
# # maximum-likelihood structure of a, alpha
# log_like <- function(a, theta, ln_y)
# {
#   d <- nrow(ln_y) # spatial dimension 
#   # cat('dimension of the spatial process', d)
#   
#   n <- ncol(ln_y) # no of observation
#   # cat('# of Obs', n)
#   
#   # initialize alpha
#   alph <- as.numeric(theta)
#   a <- as.numeric(a)
#   NC <- (d/4)^( 1/2)
#   
#   LKinfo <- LKrigSetup(subset_loc_nr,
#                        #LKGeometry="LKInterval",
#                        LKGeometry="LKRectangle",
#                        nlevel=1,
#                        a.wght=a,
#                        NC= NC,
#                        NC.buffer=0,
#                        noCheck=TRUE)
#   
#   B_a <- LKrigSAR( LKinfo, Level=1)
#   B_a <-spind2full(B_a)
#   
#   # dim(B_a)
#   # print(d)
#   
#   llike <-  0
#   for(j in 1:n)
#   {
#     llike <- llike + (d*log(alph) - alph*sum(B_a%*%log(y[1:nrow(B_a),j])))
#   }
#   # log-likelihood structure:)
#   llike <- -llike # negative likelihood structure
#   return(llike)
# }
# 
# 
# alpha_est <- rep(NA)
# neglik_vec <- rep(NA)
# for(a in 4:20)
# {
#   print(a)
#   para_start <- 1
#   optim_result <- optim(par=para_start,
#                         fn=log_like,
#                         a=a,
#                         ln_y=log(y), 
#                         # Pass y as ln_y to match the function definition
#                         method="BFGS")
#   
#   alpha_est[a] <- optim_result$par
#   neglik_vec[a] <- optim_result$value
# }
# 
# #control=list(parscale=c(0.01, 0.1, 0.01)))
# ## alpha hat <- 4.182372e+14 
# ## a = -8.960444e+14
# plot(neglik_vec,
#      pch=19)
# title(main='neg log-likelihood value', 
#       line=2)
# plot(alpha_est,
#      pch=19)
# title(main='alpha hat', 
#       line=2)
# 
# cat('alpha for min neg-like',
#     alpha_est[4],
#     'neg lik val',
#     neglik_vec[4],
#     'a.weight', a=5)
# 
# dist_matrix <- dist(subset_loc_nr)
# dim(dist_matrix)
# 
# compute_spatial_lags <- function(y)
# {
#   d <- nrow(y)
#   qu <- 
#     chi_hat <- rep(0, d)
#   for (s in 1:d)
#   {
#     chi_hat[s] <- ifelse(y[s + hs] > qu & y[s] > qu, 1, 0)
#   }
#   return(sum(chi_hat))
# }
# 
# 
# # Compute spatial lags
# spatial_lags <- compute_spatial_lags(x, hs, ht, qu)
# print(spatial_lags)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # 
# # ## -- Model Approach -- ##
# # # r = GM(.)
# # # u = 95th percentile
# # # log(Y_r) = log(R) + log(W)
# # # where log(W) = G - G_bar + const(gamma, semi-variogram(exponential covariance for the moment))
# # 
# # # (1)
# # # # # Fit GPD to the marginal, on each grid on the temporal:
# # # # # which is 31 yrs of avg temperature observation from avg_seasonality: 92 observations per year
# # # # # Fitting extremes to first 100 locations
# # # parameter_gpd <- matrix(NA,
# # #                         nrow=nrow(data_matrix),
# # #                         ncol=2) # save the scale, shape parameter
# # # tic()
# # # for(i in 1:nrow(data_matrix))
# # # {
# # #   print(i)
# # #   data_i_loc <- as.numeric(data_matrix[i,])
# # #   
# # #   # threshold
# # #   # Determine an appropriate threshold 'u' based on your data
# # #   u <- as.numeric(quantile(data_i_loc, prbs=0.90))
# # # 
# # #   # Fit the GPD model using the threshold 'u'
# # #   fit_gpd <- gpd.fit(x=data_i_loc,
# # #                      threshold=u,
# # #                      npy=92)
# # #   parameter_gpd[i,] <- fit_gpd$mle
# # # }
# # # toc()
# # # # 1460.925 sec elapsed
# # # # 24.34875 mins
# # # setwd("~/Desktop/Paper 2-dataset/Red Sea ")
# # # write.csv(parameter_gpd, "parameterGPD.csv", row.names=FALSE)
# # 
# # setwd("~/Desktop/Paper 2-dataset/Red Sea ")
# # parameter_gpd <- read.csv("parameterGPD.csv",
# #                           header=TRUE)
# # 
# # # Converting the marginals to "Standard Pareto":
# # # by dividing by the its scale parameter
# # data_matrix_scaled <- data_matrix/parameter_gpd[,1]
# # dim(data_matrix_scaled) # 16703x11315
# # 
# # # Convert the 2D matrix to a 3D array
# # days <- 365  # The number of 2D matrices in the 3D array
# # space <- nrow(data_matrix_scaled)
# # years <- n.years
# # 
# # # Create the 3D array
# # matrix_3d <- array(NA, dim=c(space, years, days))
# # risk <- matrix(NA, nrow=years, ncol=days)
# # tic()
# # for(i in 1:space)
# # {
# #   a <- 1 # lower limit
# #   b <- 365 # upper limit
# #   if(i%%1225==0)
# #   {
# #     print(i)
# #   }
# #   
# #   for(t in 1:years)
# #   {
# #     matrix_3d[i, t, ] <- data_matrix_scaled[i, a:b]
# #     #risk[t,] <- as.numeric(exp(apply(log(data_matrix_scaled[, a:b]), 2, FUN=mean)))
# #     # calculate the risk at across the space for each realization
# #     if(t<n.years)
# #     {
# #       a=b+1
# #       b=(b+365)
# #     }
# #   }
# # }
# # toc()
# # 
# # # Fix the threshold first: u
# # thresholds <- apply(data_matrix_scaled,
# #                     MARGIN=2,
# #                     FUN=quantile,
# #                     probs = 0.95)
# # u <- min(thresholds)
# # 
# # # Converting Standard Pareto Distribution to Y_r: r-Pareto Process
# # # where r() is the risk function for our case is the Geometric Mean 
# # # we work with spatial threshold, such that each spacial threshold
# # # carries atleast 10 observations per year and union of that set 
# # # Y_r process:
# # risk_across_time <- matrix(NA,
# #                            nrow=years,
# #                            ncol=days)
# # 
# # for(t in 1:years)
# # { 
# #  for(d in 1:days)
# #  {
# #    risk_across_time[t, d] <- as.numeric(exp(mean(log(matrix_3d[,t,d]))))
# #  }
# # }
# # 
# # 
# # # Pick the observation for which the risk exceeds the thresholds:
# # index_threshold_exceedance <- matrix(NA,
# #                                      nrow= years,
# #                                      ncol=days)
# # for(t in 1:years)
# # { 
# #   for(d in 1:days)
# #   {
# #     index_threshold_exceedance[t, d] <- (risk_across_time[t, d]>u)
# #     
# #   }
# # }
# # # Check: 
# # # apply(index_threshold_exceedance, 1, sum)
# # 
# # 
# # # Looking at the threshold exceedance index for each year
# # threshold_exceedance_index <- list()
# # tic()
# # for(t in 1:n.years)
# # {
# #   threshold_exceedance_index[[t]] <- (which(index_threshold_exceedance[1,]==TRUE))
# # }
# # toc()
# # # 1.784 sec elapsed
# # 
# # #  Note: for year each, I have set of days for which the risk exceeds the threshold.
# # 
# # # Combine the values of rows into a single vector
# # combined_vector <- do.call(c, threshold_exceedance_index)
# # working_index <- unique(combined_vector) # common working index
# # print(working_index)
# # length(working_index)
# # 
# # ## Y_r: r-Pareto Processes
# # # XP|r(XP)>u=Y_r
# # # Generate r-Pareto process:
# # Y_r <- array(NA, dim=c(space, years, length(working_index)))
# # for(i in 1:space)
# # {
# #   if(i%%1225==0)
# #   {
# #     print(i)
# #   }
# #   for(t in 1:years)
# #   {
# #    Y_r[i,t,] <- matrix_3d[i,t,working_index]/u
# #   }
# # }
# # # 3D array= Y_r
# # dim_3d <- dim(Y_r)
# # print(dim_3d) # 16703    31   136
# # 
# # dim_2d <- c(dim_3d[1],prod(dim_3d[-1]))
# # print(dim_2d) # 16703  4216
# # 
# # y_r <- matrix(Y_r, nrow= dim_2d[1], byrow = TRUE)
# # dim(y_r) # 16703  4216
# # 
# # # Y_r= RW, where 
# # R <- as.numeric(exp(apply(log(y_r), 1, FUN=mean)))
# # length(R)
# # 
# # log_R <- log(R)
# # length(log_R)
# # 
# # hist(R)
# # hist(log_R)
# # 
# # log_y <- log(y_r)
# # dim(log_y)
# # hist(log_y)
# # 
# # # Assume to fit a Gaussian Processes to log W
# # log_W <- (log_y-log_R)
# # hist(log_W)
# # 
# # # FIT SAR model to log_W with exponential covariance term:
# # 
# 
# #
# # ts.plot(data_matrix[500,1:(2*365)], 
# #         ylab='Loc: 500')
# # 
# # ts.plot(data_matrix[5000,],
# #         ylab='Loc: 5000')
# # 
# # ts.plot(data_matrix[10000,],
# #         ylab='Loc: 10000')
# 
# # Clear we see some seasonal patterns here.
# 
# # Deseasonality:
# # Assuming seasonal_period is the number of periods in a season (e.g., 12 for monthly data)
# # deseasonalized_data <- time_series_data - c(rep(NA, seasonal_period), time_series_data[1:(length(time_series_data) - seasonal_period)])
# 
