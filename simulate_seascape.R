####################################################################################################
# The aims of this code are:
#1) simulating a landscape/seascape of autocorrelated noise to represent the variable of trawlability
#2) simulating populations correlated to that underlying environment at various strengths
#3) simulating sampling of these populations according to stratified random design to mimic the NMFS survey
#4) saving these simulated datasets for generating model based abundance in VAST 
#
# Author: Madison B Hall / Jan 2021
####################################################################################################
library(raster)
library(gstat)
library(faux)
library(tidyr)
library(plyr)
library(dplyr)
library(survey)
library(rgeos)
library(rgdal)

#####################################################################################################

#importing most recent GOA sampling grid
grid<-shapefile("GOA_GRID_2019/grid2019.shp")

#creating a goa raster of depth data
load('Extrapolation_Depths.RData')

goa = SpatialPointsDataFrame(
  coords = Extrapolation_depths[,c('E_km', 'N_km')],
  data = data.frame(depth= Extrapolation_depths$depth) )
goa_ras = raster(goa, resolution = 5) #each grid cell is 5km X 5km
goa_ras =rasterize(x = goa, y = goa_ras, field = 'depth')
plot(goa_ras, col = rev(terrain.colors(1000)), axes = F)

########################
# below from http://santiago.begueria.es/2010/10/generating-spatially-correlated-random-fields-with-r/
# also see https://www.r-bloggers.com/2019/07/exploring-spatial-autocorrelation-in-r/

#create imaginary 10K cell grid seascape
xy <- expand.grid(1:100, 1:100)
names(xy) <- c('x','y')
g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, model=vgm(psill=0.025, range=5, model='Exp'), nmax=20)
gyy <- predict(g.dummy, newdata=xy, nsim=4)
gridded(yy) = ~x+y
spplot(obj=gyy[1]
       )
#coarser autocorrelation pattern
g.dummy1 <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, model=vgm(psill=0.025, range=20, model='Exp'), nmax=20)
gyy1 <- predict(g.dummy1, newdata=xy, nsim=4)
gridded(gyy1) = ~x+y
spplot(obj=gyy1[1])


#ID descibes the grid cell ID for every cell in the GOA sampling grid
# formatted as column - row. The separate fxn splits the ID at "-". I am trying to use the column
# and row numbers in the same way that we used x and y for the spatial pixels above
idlist<-grid@data[["ID"]]
idlist<-as.data.frame(idlist)
idlists<-separate(idlist, idlist[1], sep = "-", into = c("x", "y"))
idlists$x<-as.numeric(idlists$x)
idlists$y<-as.numeric(idlists$y)
#create df of grid id with no NAs or fxn will fail
idnona_x<-idlists$x[!is.na(idlists$x)]
idnona_y<-idlists$y[!is.na(idlists$y)]
id_nona_df<-as.data.frame(cbind(idnona_x, idnona_y))
#rename columns to match gstat formula above
names(id_nona_df)<- c("x","y")
yy <- predict(g.dummy, newdata=id_nona_df, nsim=4)
gridded(yy) = ~x+y
spplot(obj=yy)
#SUCCESS 

#each layer of yy is a different simulated seascape
#spr1<-raster(yy, layer=1, values=T)
#spr2<-raster(yy, layer=2, values=T)
#spr3<-raster(yy, layer=3, values=T)
#spr4<-raster(yy, layer=4, values=T)


df1<-as.data.frame(yy[1])
df2<-as.data.frame(yy[2])
df3<-as.data.frame(yy[3])
df4<-as.data.frame(yy[4])

#### GRAB GRID CELL CENTROIDS
grid_centers <- SpatialPointsDataFrame(gCentroid(grid, byid=TRUE), 
                                       grid@data, match.ID=FALSE)

griddf<-as.data.frame(grid_centers)
names(griddf)[10:11]<-c("centroid_x", "centroid_y")
count(griddf, "trawlable")

grid_nona<-griddf[!is.na(griddf$ID),]
count(grid_nona, "trawlable")
head(griddf)


# Pulling in data from Zack's abundance models, will use 
# the mean and stdev in order to simulate the imaginary pop

ZO_estimates<-read.csv("ZO_estimates.csv")
mu2013<-mean(ZO_estimates$polyspinis_2013)
sd2013<-sd(ZO_estimates$polyspinis_2013)

#created ID column and join simulated data df with grid 
df1$ID<-paste(df1$x, df1$y, sep="-")
df1<-as.data.frame(df1)
joined_df1<-right_join(df1, grid_nona, by = "ID")
joined_df1<-unique(joined_df1)
# approx 22% of resolved GOA grid cells are untrawlable
# setting the top 22% of simulated cells as untrawlable (4750/21588)
df1_desc<-joined_df1[order(-joined_df1$sim1),]
df1_desc$row_number <- seq.int(nrow(df1_desc))
df1_desc$trawlable<- 1
df1_desc$trawlable[df1_desc$row_number < 4750] <- -1
w=as.data.frame(table(df1_desc$STRATUM))
names(w)[1]<- "STRATUM"
w$STRATUM<-as.character(w$STRATUM)
w$sample<-w$Freq*0.025
w$sample<-round(w$sample, 0)
df1_final<-right_join(w, df1_desc, by = "STRATUM")

mu1<-mean(df1_final$sim1)
std1<-sd(df1_final$sim1)
# working with first sim enviro for now...
data1<-df1_final$sim1

# the data is not normal so simulating population generates negative
# abundance numbers in some cells. This problem might be fixable by using a diff distribution fxn, e.g. lognormal
# and there is the rlnorm fxn for simulating lognormal data, but it doesn't work here bcs this fxn
# doesn't generate (pop) data correlated to the 
# simulated environmental data, so sticking with rnrom_pre and
# just replacing negative abundance cells with 0 for now.. will have 
# to deal with this later

#below code generates a vector with defined mean and stdev, correlated at value r to dataset data1

v1 <- rnorm_pre(data1, mu = mu2013, sd = sd2013, r = 0.10)
v1[v1<0]<- 0
v2 <- rnorm_pre(data1, mu = mu2013, sd = sd2013, r = 0.25)
v2[v2<0]<- 0
v3 <- rnorm_pre(data1, mu = mu2013, sd = sd2013, r = 0.50)
v3[v3<0]<- 0
v4<- rnorm_pre(data1, mu = mu2013, sd = sd2013, r = 0.75)
v4[v4<0]<-0
v5<- rnorm_pre(data1, mu = mu2013, sd =sd2013, r = 0.90)
v5[v5<0]<-0

#pull generated population data into df with sim enviro, label columns with corr value
df1_final$pop10<-v1
df1_final$pop25<-v2
df1_final$pop50<-v3
df1_final$pop75<-v4
df1_final$pop90<-v5





#the sample size numbers below are equal to roughly 2.5% of the totaly number
#of grid cells in each stratum. A 2 boat survey typically selects 550 which
# is approx 2.5% of the total number of grid cells. This method selects 
# a total of 539
mydesigna<- sampling::strata(df1_final, stratanames=c("STRATUM"),
                             size= c(14, 16,  6, 14,  3, 13, 14, 14,  8, 18, 12,  9, 13, 15,  9, 11,  8, 10, 10,  7, 10, 15,  9,  8, 14,
                                     12, 4, 10,  5,  4,  6,  7,  3,  7, 12, 19,  8, 12,  4, 12,  7,  2,  6,  6,  4,  4, 15, 10,  8,  5,
                                     6, 6, 5, 16, 8, 6, 8,  7, 5), 
                             method=c("srswor"), description = T)

df_trawlonly<-df1_final[(df1_final$trawlable == 1), ]

mydesignb<- sampling::strata(df_trawlonly, stratanames=c("STRATUM"),
                             size= c(14, 16,  6, 14,  3, 13, 14, 14,  8, 18, 12,  9, 13, 15,  9, 11,  8, 10, 10,  7, 10, 15,  9,  8, 14,
                                     12, 4, 10,  5,  4,  6,  7,  3,  7, 12, 19,  8, 12,  4, 12,  7,  2,  6,  6,  4,  4, 15, 10,  8,  5,
                                     6, 6, 5, 16, 8, 6, 8,  7, 5), 
                             method=c("srswor"), description = T)

surveydata_both<-sampling::getdata(df1_final, mydesigna)
names(surveydata_both)[25]<-"STRAT_COUNT"
surveydata_trawlonly<-sampling::getdata(df_trawlonly, mydesignb)
names(surveydata_trawlonly)[25]<-"STRAT_COUNT"

count(surveydata_both$STRATUM)
count(surveydata_trawlonly$STRATUM)

############################################################################
#### simulating area swept #################################################
############################################################################
realdata<-read.csv("survey_area_swept.csv")
realdata$AREA_SWEPT_KM2<-realdata$DISTANCE_FISHED*realdata$NET_WIDTH.KM.
realdata<-realdata[(realdata$YEAR > 2010), ]
shapiro.test(realdata$AREA_SWEPT_KM2)
hist(realdata$AREA_SWEPT_KM2, breaks = 20)
areamu<-mean(realdata$AREA_SWEPT_KM2)
areasd<-sd(realdata$AREA_SWEPT_KM2)
sim_area_swept<-rnorm(539, areamu, areasd)

surveydata_both$sim_area_swept_km2<-sim_area_swept
# remove rows where area swept exceeds the actual area of the grid cell...
#(these are very small slivers of grid cells bisected by Statum lines)
surveydata_both<-surveydata_both[(surveydata_both$AREA_KM2 > surveydata_both$sim_area_swept_km2), ]
surveydata_both$area_conversion<-(surveydata_both$sim_area_swept_km2/surveydata_both$AREA_KM2)
surveydata_both$sim_survey_pop10_perfect<-(surveydata_both$pop10*surveydata_both$area_conversion)
surveydata_both$sim_survey_pop25_perfect<-(surveydata_both$pop25*surveydata_both$area_conversion)
surveydata_both$sim_survey_pop50_perfect<-(surveydata_both$pop50*surveydata_both$area_conversion)
surveydata_both$sim_survey_pop75_perfect<-(surveydata_both$pop75*surveydata_both$area_conversion)
surveydata_both$sim_survey_pop90_perfect<-(surveydata_both$pop90*surveydata_both$area_conversion)


#repeat for trawl only
surveydata_trawlonly$sim_area_swept_km2<-sim_area_swept
surveydata_trawlonly<-surveydata_trawlonly[(surveydata_trawlonly$AREA_KM2 >surveydata_trawlonly$sim_area_swept_km2), ]
surveydata_trawlonly$area_conversion<-(surveydata_trawlonly$sim_area_swept_km2/surveydata_trawlonly$AREA_KM2)
surveydata_trawlonly$sim_survey_pop10_perfect<-(surveydata_trawlonly$pop10*surveydata_trawlonly$area_conversion)
surveydata_trawlonly$sim_survey_pop25_perfect<-(surveydata_trawlonly$pop25*surveydata_trawlonly$area_conversion)
surveydata_trawlonly$sim_survey_pop50_perfect<-(surveydata_trawlonly$pop50*surveydata_trawlonly$area_conversion)
surveydata_trawlonly$sim_survey_pop75_perfect<-(surveydata_trawlonly$pop75*surveydata_trawlonly$area_conversion)
surveydata_trawlonly$sim_survey_pop90_perfect<-(surveydata_trawlonly$pop90*surveydata_trawlonly$area_conversion)
#########################################################################
###### plotting  ########################################################
#########################################################################


#separate df for plotting
plotdf<-as.data.frame(df1_final)
plotdf1<-plotdf[, (c(4:6, 14:20))]

gridded(plotdf1) = ~x+y
spplot(plotdf1[1]) #plot simulated environment
spplot(plotdf1[4:8]) #plot 5 simulated populations corr to that environment

#########################################################################
### TO BE DONE ##########################################################
#########################################################################


##############################################################################
#### gradb centroids of each grid cell #######################################
##############################################################################
# get the centroids and then convert them to a SpatialPointsDataFrame
grid_centers <- SpatialPointsDataFrame(gCentroid(grid, byid=TRUE), 
                                      grid@data, match.ID=FALSE)


head(grid_centers)
grid_centers<-as.data.frame(grid_centers)

#######################################################################################################
### TRASH CODE ########################################################################################
#######################################################################################################

library(spatialEco)
library(raster)

rr<-raster(nrow= 100, ncol= 230) #same number as GOA sampling grid
grid<-raster()

#random noise raster, but not autocorrelated
ras<-random.raster(r = goa_ras,
                   n.layers = 1,
                   min = 1,
                   max = 100,
                   distribution = "random")
plot(ras)
summary(ras)
ras

yras<-raster(ncol=100, nrow=100)
yras<-rasterize(yy[1], yras, data= "sim1")
yras<-yras$sim1
#getting cell values from raster into df format
rvec<-as.data.frame(yras$sim1)

ddply(df,.(ID),function(x) x[sample(nrow(x),500),])
#OR
new_df <- df %>% group_by(ID) %>% sample_n(500)


strata<-grid$STRATUM
strata<-strata[!is.na(strata)]








