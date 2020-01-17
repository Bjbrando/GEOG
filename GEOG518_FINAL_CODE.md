# GEOG
Autocorrelation_Code
library("lubridate")
library("maptools")
library("raster")
library("tidyverse")
library("adehabitatLT")
library("adehabitatHR")
library("chron")
library("spatstat")
library("tkrplot")
library("amt")
library("rgdal")
library("dplyr")
library("ggplot2")
library("tmap")
library("sf") 
library("move")
library("sp")
library("ecodist")
library("mapview")
library("maps")
library("ggmap")
library("ggsn")
library("grid")
library("gridExtra")
library("gtable")
library("REdaS")
library("plyr")
library("spdep")
library("GISTools")
library("tmap")
library("BAMMtools")
library("shinyjs")
library("spatstat") #d. tesselation
library("tmaptools")
library("ggspatial")
#LIKELY MORE LIBRARIES THAN NECESSARY

setwd("/Users/briannabrandon/GEOG518_Code/Final")
getwd()
#Read in muskox GPS collar data
dat<-read.csv("./YMu2015-19.csv")

#Read in Predictive Ecosystem Model ( 6m x 6m resolution)
PEM<-raster("./NS_PEM_final.tif")

#set Coordinate Reference System from PEM
#crs(PEM)
crs=CRS("+proj=utm +zone=7 +ellps=GRS80 +units=m +no_defs") #projection of PEM

#Remove letters from device name
dat$Device.Name<-as.integer(gsub('[a-zA-Z]','', dat$Device.Name))

#Remove muskox 1612 from dataset due to 10 repeated date-time stamps
dat_m<-dat %>%
  filter(Device.Name %in% c(1501, 1502, 1503, 1504, 1505, 1506, 1507, 1508,1609, 1610, 1613, 1614, 1615, 1816, 1817, 1818, 1827, 1819,1820, 1821, 1822, 1823, 1824, 1825))

#Convert dates to POSIXct
t_local<-as.character(dat_m$Date...Time..Local.)
t_local<-as.POSIXct(strptime(as.character(dat_m$t_local),"%m/%d/%Y %H:%M", tz="America/Dawson"))

t_GMT<-as.character(dat_m$Date...Time..GMT.)
t_GMT<-as.POSIXct(strptime(as.character(t_GMT), "%m/%d/%Y %H:%M", tz="GMT"))


dat_m$tGMT<-t_GMT
#====
coords<-dat_m[c("Longitude", "Latitude")]
xysp<-SpatialPoints(coords, proj4string = CRS("+init=epsg:4326"))
xysp<-spTransform(xysp,CRS("+proj=utm +zone=7 +ellps=GRS80 +units=m +no_defs"))
dat_m<-data.frame(xysp, dat_m)
#====

cols<-c(2, 7:12)

#Create ltraj 
musk<- as.ltraj(xy = dat_m[,c("Longitude","Latitude")], date =t_GMT, id = dat_m$Device.Name,  infolocs= dat_m[cols], typeII = TRUE)

foo <- function(dt) {
  return(dt> (3600*5.2))
  }

musk_5hr <- cutltraj(musk, "foo(dt)", nextr = TRUE)

#set reference date-time to begin counting
refda <- strptime("2016-05-01 00:00", "%Y-%m-%d %H:%M", tz="GMT")

#SetNA and when time interval >5hrs
musk_5hr_reg<- setNA(musk_5hr, refda, 5, units = "hour")
musk_5hr_reg <- sett0(musk_5hr_reg, refda, 5, units = "hour")


#Subset muskox id:1610 
musk1610_5hr<-musk_5hr_reg[30]

#Convert to dataframes from ltraj

df_musk1610_5hr<-ld(musk1610_5hr)

#======Calculate additional infolocs=========#
#Subtract 8 hours from GMT to get local time, excluding daylight savings

df_musk1610_5hr$date_local<-as.POSIXct(as.character(df_musk1610_5hr$date - 8*60*60), "%d/%m/%Y %H:%M")


#Classify Seasons
#Note, seasons currently classified based on GMT (Local America/ Dawson= -8 hr GMT)

getSeason <- function(DATES) {
  LW <- as.Date("2012-01-01", format = "%Y-%m-%d") # Late winter start
  CS <- as.Date("2012-04-01",  format = "%Y-%m-%d") # Calving spring start
  SS <- as.Date("2012-06-01",  format = "%Y-%m-%d") # Summer start
  FR <- as.Date("2012-08-01",  format = "%Y-%m-%d") # Fall rut start
  EW <- as.Date("2012-10-01", format = "%Y-%m-%d") # Early winter start
  
  # Convert dates from any year to 2012 dates
  d <- as.Date(strftime(DATES, format="2012-%m-%d"))
  
  ifelse (d >= LW & d < CS, "Late.Winter",
          ifelse (d >= CS & d < SS, "Calving",
                  ifelse (d >= SS & d < FR, "Summer", 
                          ifelse (d>=FR & d<EW, "Rut",
                                  "Early.Winter"))))
}


test<-as.Date(as.POSIXct(df_musk1610_5hr$date, "Y%/%m-%d", tz="GMT"), tz="GMT")
df_musk1610_5hr$Season<-getSeason(test)


#Calculate tod
df_musk1610_5hr$doy<-strftime(df_musk1610_5hr$date, format="%j")

#Extract infolocks before interpolating values removed. 
info<-df_musk1610_5hr[,c(14:22)]
#===== Screen outlyers=====#
#Find row of location with distance >1000m 
summary(df_musk1610_5hr$dist)
which(df_musk1610_5hr$dist>= 10000, arr.ind=TRUE)

#Replove row with distance >1000m
df_musk1610_5hr<-df_musk1610_5hr[-4599,]


#===Convert back into ltraj and recalculate step metrics====#
musk1610_df<-data.frame(df_musk1610_5hr) #Convert spatial points dataframe to dataframe
musk1610_traj1<-dl(df_musk1610_5hr, crs) #Convert from df to ltraj
musk1610_traj1<-rec(musk1610_traj1) #Recalculate step metrics 

#Extract infolocs with number of rows =5589
info<-infolocs(musk1610_traj1)

#Set missing fix at 5hr interval to na (rows =5590)
musk1610.traj <- setNA(musk1610_traj1, refda, 5, units = "hour")
#One na value

#Interpolate single na location (rows = 5589.. why??)
musk1610.traj <- redisltraj(na.omit(musk1610.traj[1]), (3600*5), type="time")

#Rejoin infolocs (rows= 5589)
infolocs(musk1610.traj)<-info


#plot(musk1610.traj)
#plotltr(musk1610.traj, "dt")
is.regular(musk1610.traj)
#TRUE

#===Extract habitat covariates and subset based on year===#
#convert back to dataframe
musk1610.traj_df<-ld(musk1610.traj)

#Create column for year
musk1610.traj_df$year<-as.character(format(musk1610.traj_df$date, "%Y"))

#Create spatical points dataframe
coords<-musk1610.traj_df[,c("x", "y")]

musk1610.traj_df.sp<-SpatialPointsDataFrame(coords=coords, data=musk1610.traj_df, proj4string = crs)

#Extract PEM classfication for each point
musk1610.traj_df.sp$PEM<-raster::extract(PEM, musk1610.traj_df.sp)

#Subset based on year

musk1610.traj_df2016<-musk1610.traj_df %>% filter(date >= as.Date('2016-01-01') & date <= as.Date('2016-12-31') )
musk1610.traj_df2017<-musk1610.traj_df %>% filter(date >= as.Date('2017-01-01') & date <= as.Date('2017-12-31') )
musk1610.traj_df2018<-musk1610.traj_df %>% filter(date >= as.Date('2018-01-01') & date <= as.Date('2018-12-31') )
musk1610.traj_df2019<-musk1610.traj_df %>% filter(date >= as.Date('2019-01-01') & date <= as.Date('2019-12-31') )

#Convert year subsets to itraj

musk1610.traj_2016<-dl(musk1610.traj_df2016, crs)
musk1610.traj_2017<-dl(musk1610.traj_df2017, crs)
musk1610.traj_2018<- dl(musk1610.traj_df2018, crs)
musk1610.traj_2019<-dl(musk1610.traj_df2019, crs)


#Subset 2017 into Seasons

musk1610.traj_df2017.calv<-musk1610.traj_df2017[musk1610.traj_df2017$Season=="Calving",]
musk1610.traj_df2017.sum<-musk1610.traj_df2017[musk1610.traj_df2017$Season=="Summer",]
musk1610.traj_df2017.rut<-musk1610.traj_df2017[musk1610.traj_df2017$Season=="Rut",]
musk1610.traj_df2017.ewin<-musk1610.traj_df2017[musk1610.traj_df2017$Season=="Early.Winter",]
musk1610.traj_df2017.lwin<-musk1610.traj_df2017[musk1610.traj_df2017$Season=="Late.Winter",]



#Convert 2017 seasons to ltraj
musk1610.traj_2017.calv<-dl(musk1610.traj_df2017.calv, crs)
musk1610.traj_2017.sum<-dl(musk1610.traj_df2017.sum, crs)
musk1610.traj_2017.rut<-dl(musk1610.traj_df2017.rut, crs)
musk1610.traj_2017.ewin<-dl(musk1610.traj_df2017.ewin, crs)
musk1610.traj_2017.lwin<-dl(musk1610.traj_df2017.lwin, crs)
#======VISUALIZE DATA=============#
#======SITE MAP===================#

#GGmaps
#Get maps
yukon2<-get_map(location="yukon", zoom=5, maptype="terrain", source="stamen")
terrain<-get_map(c(left = -138, bottom = 68.6, right = -136.9, top = 69), maptype="terrain", source = "stamen")

#Site map
#Calculate mean coordinates
mean(geo_plot$Latitude.1)
mean(geo_plot$Longitude.1)


jpeg("./Site_Map.jpeg")
site_map<-ggmap(yukon2)+
  geom_point(data=site, aes(x=x, y=y), size=5, shape=22, fill="red")+
  scalebar(x.min=-149, x.max= -121, y.min=57.5, y.max=69.7, 
           dist=200, dist_unit="km", border.size=0.3, st.bottom=FALSE, transform=TRUE, model="WGS84", location="bottomleft", anchor=c(x=x_scale_loc, y=y_scale_loc))
north2(site_map, 0.85, 0.9, scale=0.1, symbol= 3)
dev.off()

#=====DATA FOR MOVEMENT MAPS======#
geo_plot<-musk1610.traj_df2017[,c("Latitude.1", "Longitude.1", "Season", "date")]
geo_plot.lwin<-musk1610.traj_df2017.lwin[,c("Latitude.1", "Longitude.1", "Season", "date")]
geo_plot.calv<-musk1610.traj_df2017.calv[,c("Latitude.1", "Longitude.1", "Season", "date")]
geo_plot.sum<-musk1610.traj_df2017.sum[,c("Latitude.1", "Longitude.1", "Season", "date")]
geo_plot.rut<-musk1610.traj_df2017.rut[,c("Latitude.1", "Longitude.1", "Season", "date")]
geo_plot.ewin<-musk1610.traj_df2017.ewin[,c("Latitude.1", "Longitude.1", "Season", "date")]

#======MAP 2017=======#
#season<- c("#FF6EB4CC", "#0000EECC", "#63B8FFCC", "#00CD66CC", "#FF7F00CC")
season<-alpha(c("hotpink1", "steelblue1", "blue2","darkorange1", "springgreen3"), 0.5)


move_map2017<-ggmap(terrain)+
  geom_path(data=geo_plot, aes(x=Longitude.1, y=Latitude.1, colour=Season, show.legend=NULL), size= 0.2)+
  geom_point(data=geo_plot, aes(x=Longitude.1, y=Latitude.1, colour=Season, show.legend=NULL), size= 0.5, shape=22)+
  geom_point(data=median_centre, aes(x=long, y=lat), colour="black", size=5, shape=13)+
  scale_colour_manual(values=c("hotpink1", "steelblue1", "blue2","darkorange1", "springgreen3"))+
  theme(legend.position = "none")+
  scalebar(x.min=-138, x.max= -137, y.min=68.6, y.max=69, 
           dist=5, dist_unit="km", border.size=0.3, st.bottom=FALSE, transform=TRUE, model="WGS84", location="bottomright", anchor=c(x=x_scale_loc1, y=y_scale_loc1))

jpeg("./Move_Map2017.jpg")
move_map2017
north2(move_map2017, 0.20, 0.9, scale=0.1, symbol= 3)
dev.off()

calving_map2017<-ggmap(terrain)+
  geom_path(data=geo_plot.calv, aes(x=Longitude.1, y=Latitude.1, show.legend=NULL), size= 0.2, colour="hotpink1")+
  geom_point(data=geo_plot.calv, aes(x=Longitude.1, y=Latitude.1, show.legend=NULL), size= 0.5, shape=22, colour="hotpink1")+
  geom_point(data=median_centre.calv, aes(x=long, y=lat), colour="black", size=5, shape=13)+
  theme(legend.position="none")+
  scalebar(x.min=-138, x.max= -137, y.min=68.6, y.max=69, 
           dist=5, dist_unit="km", border.size=0.3, st.bottom=FALSE, transform=TRUE, model="WGS84", location="bottomright", anchor=c(x=x_scale_loc1, y=y_scale_loc1))

jpeg("./Calving_Map2017.jpg")
calving_map2017
north2(calving_map2017, 0.20, 0.9, scale=0.1, symbol= 3)
dev.off()

summer_map2017<-ggmap(terrain)+
  geom_path(data=geo_plot.sum, aes(x=Longitude.1, y=Latitude.1, show.legend=NULL), size= 0.2, colour="springgreen3")+
  geom_point(data=geo_plot.sum, aes(x=Longitude.1, y=Latitude.1, show.legend=NULL), size= 0.5, shape=22, colour="springgreen3")+
  geom_point(data=median_centre.sum, aes(x=long, y=lat), colour="black", size=5, shape=13)+
  theme(legend.position="none")+
  scalebar(x.min=-138, x.max= -137, y.min=68.6, y.max=69, 
           dist=5, dist_unit="km", border.size=0.3, st.bottom=FALSE, transform=TRUE, model="WGS84", location="bottomright", anchor=c(x=x_scale_loc1, y=y_scale_loc1))

jpeg("./Summer_Map2017.jpg")
summer_map2017
north2(summer_map2017, 0.20, 0.9, scale=0.1, symbol= 3)
dev.off()

rut_map2017<-ggmap(terrain)+
  geom_path(data=geo_plot.rut, aes(x=Longitude.1, y=Latitude.1, show.legend=NULL), size= 0.2, colour="darkorange1")+
  geom_point(data=geo_plot.rut, aes(x=Longitude.1, y=Latitude.1, show.legend=NULL), size= 0.5, shape=22, colour="darkorange1")+
  geom_point(data=median_centre.rut, aes(x=long, y=lat), colour="black", size=5, shape=13)+
  theme(legend.position="none")+
  scalebar(x.min=-138, x.max= -137, y.min=68.6, y.max=69, 
           dist=5, dist_unit="km", border.size=0.3, st.bottom=FALSE, transform=TRUE, model="WGS84", location="bottomright", anchor=c(x=x_scale_loc1, y=y_scale_loc1))

jpeg("./Rut_Map2017.jpg")
rut_map2017
north2(rut_map2017, 0.20, 0.9, scale=0.1, symbol= 3)
dev.off()

e.wint_map2017<-ggmap(terrain)+
  geom_path(data=geo_plot.ewin, aes(x=Longitude.1, y=Latitude.1, show.legend=NULL), size= 0.2, colour="steelblue1")+
  geom_point(data=geo_plot.ewin, aes(x=Longitude.1, y=Latitude.1, show.legend=NULL), size= 0.5, shape=22, colour="steelblue1")+
  geom_point(data=median_centre.ew, aes(x=long, y=lat), colour="black", size=5, shape=13)+
  theme(legend.position="none")+
  scalebar(x.min=-138, x.max= -137, y.min=68.6, y.max=69, 
           dist=5, dist_unit="km", border.size=0.3, st.bottom=FALSE, transform=TRUE, model="WGS84", location="bottomright", anchor=c(x=x_scale_loc1, y=y_scale_loc1))

jpeg("./eWint_Map2017.jpg")
e.wint_map2017
north2(e.wint_map2017, 0.20, 0.9, scale=0.1, symbol= 3)
dev.off()


lwint_map2017<-ggmap(terrain)+
  geom_path(data=geo_plot.lwin, aes(x=Longitude.1, y=Latitude.1, show.legend=NULL), size= 0.2, colour="dodgerblue4")+
  geom_point(data=geo_plot.lwin, aes(x=Longitude.1, y=Latitude.1, show.legend=NULL), size= 0.5, shape=22, colour="dodgerblue4")+
  geom_point(data=median_centre.lw, aes(x=long, y=lat), colour="red", size=5, shape=13)+
  theme(legend.position="none")+
  scalebar(x.min=-138, x.max= -137, y.min=68.6, y.max=69, 
           dist=5, dist_unit="km", border.size=0.3, st.bottom=FALSE, transform=TRUE, model="WGS84", location="bottomright", anchor=c(x=x_scale_loc1, y=y_scale_loc1))

jpeg("./lWint_Map2017.jpg")
lwint_map2017
north2(lwint_map2017, 0.20, 0.9, scale=0.1, symbol= 3)
dev.off()




#====PLOT TRAJECTORIES=====#
#Calculate median center for each season 
med_centre_LAT=median(geo_plot$Latitude.1)
med_centre_LONG=median(geo_plot$Longitude.1)

med_centre_LAT.lw=median(geo_plot.lwin$Latitude.1)
med_centre_LONG.lw=median(geo_plot.lwin$Longitude.1)

med_centre_LAT.calv=median(geo_plot.calv$Latitude.1)
med_centre_LONG.calv=median(geo_plot.calv$Longitude.1)

med_centre_LAT.sum=median(geo_plot.sum$Latitude.1)
med_centre_LONG.sum=median(geo_plot.sum$Longitude.1)

med_centre_LAT.rut=median(geo_plot.rut$Latitude.1)
med_centre_LONG.rut=median(geo_plot.rut$Longitude.1)

med_centre_LAT.ew=median(geo_plot.ewin$Latitude.1)
med_centre_LONG.ew=median(geo_plot.ewin$Longitude.1)

median_centre<-data.frame(name="median_centre", long=med_centre_LONG, lat=med_centre_LAT)
coords_median<-median_centre[,c("long", "lat")] #seperate lat and long

median_centre.lw<-data.frame(name="median_centre", long=med_centre_LONG.lw, lat=med_centre_LAT.lw)
coords_median.lw<-median_centre.lw[,c("long", "lat")] #seperate lat and long
#medianPoint<-SpatialPointsDataFrame(coords=coords_median, data=median_centre, proj4string=crs)
median_centre.calv<-data.frame(name="median_centre", long=med_centre_LONG.calv, lat=med_centre_LAT.calv)
coords_median.calv<-median_centre.calv[,c("long", "lat")] #seperate lat and long

median_centre.sum<-data.frame(name="median_centre", long=med_centre_LONG.sum, lat=med_centre_LAT.sum)
coords_median.sum<-median_centre.sum[,c("long", "lat")] #seperate lat and long

median_centre.rut<-data.frame(name="median_centre", long=med_centre_LONG.rut, lat=med_centre_LAT.rut)
coords_median.rut<-median_centre.rut[,c("long", "lat")] #seperate lat and long

median_centre.ew<-data.frame(name="median_centre", long=med_centre_LONG.ew, lat=med_centre_LAT.ew)
coords_median.ew<-median_centre.ew[,c("long", "lat")] #seperate lat and long








#======PLOT TRAJECTORIES INDEPENDENTY===========#
jpeg("./traj_2017.jpg")
plot(musk1610.traj_2017, main="2017", xlab="Easting (mE)", ylab="Northing (mN)")
dev.off()

jpeg("./traj_2017calv.jpg")
plot(musk1610.traj_2017.calv, main="Calving", xlab="Easting (mE)", ylab="Northing (mN)")
dev.off()

jpeg("./traj_2017sum.jpg")
plot(musk1610.traj_2017.sum, main="Summer", xlab="Easting (mE)", ylab="Northing (mN)")
dev.off()

jpeg("./traj_2017rut.jpg")
plot(musk1610.traj_2017.rut, main="Fall Rut", xlab="Easting (mE)", ylab="Northing (mN)")
dev.off()

jpeg("./traj_2017ewin.jpg")
plot(musk1610.traj_2017.ewin, main="Early Winter", xlab="Easting (mE)", ylab="Northing (mN)")
dev.off()

jpeg("./traj_2017lwin.jpg")
plot(musk1610.traj_2017.lwin, main="Late Winter", xlab="Easting (mE)", ylab="Northing (mN)")
dev.off()

#====CALCULATE RESIDENCE TIMES=====#
#(Barraquand & Benhamou 2008). Two parameters are required to determine RT, namely the circle’s radius and the cut-off time. TtoR is the uninterrupted time an animal spends before its first return to the circle centred on each location. Returns within a time shorter than the cut-off time are not classified as returns; instead, such a ‘return’ is considered part of the uninterrupted residency. We


plot(musk1610.traj_2017.calv)

#calculate RT over timeseries
#Patch size(radius=1000m), maximum time to return (maxt=10hrs)

res.1000.10 <- residenceTime(musk1610.traj_2017, radius = 1000, maxt=10, units="hour", addinfo = TRUE)

res.1000.10.calv <- residenceTime(musk1610.traj_2017.calv, radius = 1000, maxt=10, units="hour", addinfo = TRUE)

res.1000.10.sum <- residenceTime(musk1610.traj_2017.sum, radius = 1000, maxt=10, units="hour", addinfo = TRUE)

res.1000.10.rut <- residenceTime(musk1610.traj_2017.rut, radius = 1000, maxt=10, units="hour", addinfo = TRUE)

res.1000.10.ewin <- residenceTime(musk1610.traj_2017.ewin, radius = 1000, maxt=10, units="hour", addinfo = TRUE)

res.1000.10.lwin <- residenceTime(musk1610.traj_2017.lwin, radius = 1000, maxt=10, units="hour", addinfo = TRUE)

#Convert to dataframes

res.1000.10_df<-ld(res.1000.10)

res.1000.10.calv_df<-ld(res.1000.10.calv)

res.1000.10.sum_df<-ld(res.1000.10.sum)

res.1000.10.rut_df<-ld(res.1000.10.rut)

res.1000.10.ewin_df<-ld(res.1000.10.ewin)

res.1000.10.lwin_df<-ld(res.1000.10.lwin)


#=====histograms (dist and residence time)

hist.dist_2017<-ggplot(res.1000.10_df, aes(x=dist))+geom_histogram()

jpeg("./musk_explore/hist.dist_2017.jpg")
hist_dist
dev.off()


hist.RT_2017<-ggplot(res.1000.10_df, aes(x=RT.1000))+geom_histogram()

jpeg("./musk_explore/hist_RT_2017.jpg")
hist.RT_2017
dev.off()


#COLOUR SCHEMES

season<-alpha(c("hotpink1", "steelblue1", "blue2","darkorange1", "springgreen3"), 0.5)

#Boxplots of years + dist, rel.angle, abs.angle

geom_boxplot(outlier.colour="black", outlier.shape=16,
             outlier.size=2, notch=FALSE, show.legend = FALSE)

box_dist_year<-ggplot(res.1000.10_df, aes(x=Season, y=dist, fill=Season))+scale_fill_manual(values=alpha(season, 0.8))+geom_boxplot()+theme(legend.position="none")

jpeg("./musk_explore/box_dist2017.jpg")
box_dist_year
dev.off()

box_RT.1000_year<-ggplot(res.1000.10_df, aes(x=Season, y=RT.1000, fill=Season))+scale_fill_manual(values=alpha(season, 0.8))+geom_boxplot()+theme(legend.position="none")

jpeg("./musk_explore/box_RT.1000_year.jpg")
box_RT.1000_year
dev.off()


scatter_doy_dist<-ggplot(res.1000.10_df, aes(x=doy, y=dist, colour=Season))+geom_point()+ scale_x_discrete(breaks=c(180,360))+ scale_colour_manual(values=season)

jpeg("./musk_explore/scatter_doy_dist.jpg")
scatter_doy_dist
dev.off()


scatter_doy_RT.1000<-ggplot(res.1000.10_df, aes(x=doy, y=RT.1000, colour=Season))+geom_point()+ scale_x_discrete(breaks=c(180,360))+ scale_colour_manual(values=season)

jpeg("./musk_explore/scatter_doy_RT.1000.jpg")
scatter_doy_RT.1000
dev.off()
#======SUMMARY STATS=====#
#Distance
table_summary.all<-res.1000.10_df %>%
  dplyr::summarise(
    count = n(),
    mean = mean(dist, na.rm = TRUE),
    sd = sd(dist, na.rm = TRUE),
    median = median(dist, na.rm = TRUE)
  )

table_summary.all<-data.frame(table_summary.all)
table_summary.all$Season<-"Full Year"
table_summary.all<-table_summary.all[,c(5, 1, 2, 3, 4)]


table_summary.season<-dplyr::group_by(res.1000.10_df, Season) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(dist, na.rm = TRUE),
    sd = sd(dist, na.rm = TRUE),
    median = median(dist, na.rm = TRUE)
  )
table_summary<-rbind(table_summary.all, table_summary.season)

table_summary<-data.frame(rbind(table_summary.all, table_summary.season))

table_summary$CoV<-(table_summary$sd/table_summary$mean)*100
table_summary<-table_summary[c("Season", "count", "mean", "median", "sd", "CoV")]

Summary_Stats<-tableGrob(table_summary, rows= c("Full Year","Calving", "Early Winter", "Late Winter", "Rut", "Summer"), cols=c("Season", "Count", "Mean", "Median", "sd", "CoV"))
Caption<-textGrob("Table 1: Summary Statistics for Step Distance \nrecorded at regular 5hr intervals in 2017", gp= gpar(fontsize=10))
padding<-unit(5,"mm")
Summary_Stats<-gtable_add_rows(Summary_Stats,
                               heights=grobHeight(Caption)+padding,
                               pos=0)
Summary_Stats_dist<- gtable_add_grob(Summary_Stats,
                                     Caption, t = 1, l = 2, r = ncol(table_summary)+1)

dev.off()
grid.arrange(Summary_Stats_dist,newpage=TRUE)


png("./Summary_Stats_dist.png")
Summary_Stats_dist
grid.arrange(Summary_Stats_dist,newpage=TRUE)
dev.off()

#Residence Time
table_summary.all<-res.1000.10_df %>%
  dplyr::summarise(
    count = n(),
    mean = mean(RT.1000, na.rm = TRUE),
    sd = sd(RT.1000, na.rm = TRUE),
    median = median(RT.1000, na.rm = TRUE)
  )

table_summary.all<-data.frame(table_summary.all)
table_summary.all$Season<-"Full Year"
table_summary.all<-table_summary.all[,c(5, 1, 2, 3, 4)]


table_summary.season<-dplyr::group_by(res.1000.10_df, Season) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(RT.1000, na.rm = TRUE),
    sd = sd(RT.1000, na.rm = TRUE),
    median = median(RT.1000, na.rm = TRUE)
  )
table_summary<-rbind(table_summary.all, table_summary.season)
table_summary<-data.frame(rbind(table_summary.all, table_summary.season))
table_summary$CoV<-(table_summary$sd/table_summary$mean)*100
table_summary<-table_summary[c("Season", "count", "mean", "median", "sd", "CoV")]

Summary_Stats<-tableGrob(table_summary, rows= c("Full Year","Calving", "Early Winter", "Late Winter", "Rut", "Summer"), cols=c("Season", "Count", "Mean", "Median", "sd", "CoV"))
Caption<-textGrob("Table 2: Summary statistics for \nResidence Time (Area radius=1000m, ToR=10hr) in 2017", gp= gpar(fontsize=10))
padding<-unit(5,"mm")
Summary_Stats<-gtable_add_rows(Summary_Stats,
                               heights=grobHeight(Caption)+padding,
                               pos=0)
Summary_Stats_RT<- gtable_add_grob(Summary_Stats,
                                   Caption, t = 1, l = 2, r = ncol(table_summary)+1)

dev.off()
grid.arrange(Summary_Stats_RT,newpage=TRUE)


png("./musk_explore/Summary_Stats_RT.png")
Summary_Stats_RT
grid.arrange(Summary_Stats_RT,newpage=TRUE)
dev.off()


#====Is difference in distance travelled significantly different between seasons?===#
#Use Kruskal-Wallis test

kruskal.test(RT.1000~Season, data=res.1000.10_df)
#Kruskal-Wallis rank sum test

#data:  dist by Season
#Kruskal-Wallis chi-squared = 538.57, df = 4, p-value < 2.2e-16

#Kruskal-Wallis rank sum test

#data:  RT.1000 by Season
#Kruskal-Wallis chi-squared = 1021.3, df = 4, p-value < 2.2e-16


pairwise.wilcox.test(res.1000.10_df$RT.1000, res.1000.10_df$Season,
                     p.adjust.method = "BH")

#Pairwise comparisons using Wilcoxon rank sum test 

#data:  musk1610.traj_even$dist and musk1610.traj_even$Even.Year 

#         2016-17 2017-18
# 2017-18 0.25    -      
# 2018-19 6.1e-10 2.0e-07


#Distance in 2018-19 was significantly different than 2016-17 and 2017-18. 


#======plot temporal autocorrelation of Residence Time with 95% CI


acf_RT_2017<-acf(res.1000.10_df$RT.1000, n.sim=999, lag.max=20, na.action=na.pass, level=95, plot= FALSE)


jpeg("./ResidenceT/acf_RT.2017.jpg")
plot(acf_RT_2017, main="2017")
dev.off()


acf_RT_calv<-acf(res.1000.10.calv_df$RT.1000, n.sim=999, lag.max=20, na.action=na.pass, level=95, plot=FALSE)
jpeg("./ResidenceT/acf_RT.calv.jpg")
plot(acf_RT_calv, main="Spring Calving")
dev.off()


acf_RT_sum<-acf(res.1000.10.sum_df$RT.1000, n.sim=999, lag.max=20, na.action=na.pass, level=95, plot=FALSE)

jpeg("./ResidenceT/acf_RT.sum.jpg")
plot(acf_RT_sum, main="Summer")
dev.off()


acf_RT_rut<-acf(res.1000.10.rut_df$RT.1000, n.sim=999, lag.max=20, na.action=na.pass, level=95, plot=FALSE)

jpeg("./ResidenceT/acf_RT.rut.jpg")
plot(acf_RT_rut, main="Fall Rut")
dev.off()


acf_RT_ewin<-acf(res.1000.10.lwin_df$RT.1000, n.sim=999, lag.max=20, na.action=na.pass, level=95, plot=FALSE)
jpeg("./ResidenceT/acf_RT.ewin.jpg")
plot(acf_RT_ewin, main="Early Winter")
dev.off()

acf_RT_lwin<-acf(res.1000.10.ewin_df$RT.1000, n.sim=999, lag.max=20, na.action=na.pass, level=95, plot=FALSE)

jpeg("./ResidenceT/acf_RT.lwin.jpg")
plot(acf_RT_lwin, main="Late Winter")
dev.off()

#====acf for distance===#


acf_dist_2017<-acf(res.1000.10_df$dist, n.sim=999, lag.max=20, na.action=na.pass, level=95, plot= FALSE)


jpeg("./acf_dist.2017.jpg")
plot(acf_dist_2017, main="2017")
dev.off()


acf_dist_calv<-acf(res.1000.10.calv_df$dist, n.sim=999, lag.max=20, na.action=na.pass, level=95, plot=FALSE)
jpeg("./acf_dist.calv.jpg")
plot(acf_dist_calv, main="Spring Calving")
dev.off()


acf_dist_sum<-acf(res.1000.10.sum_df$dist, n.sim=999, lag.max=20, na.action=na.pass, level=95, plot=FALSE)

jpeg(".//acf_dist.sum.jpg")
plot(acf_dist_sum, main="Summer")
dev.off()


acf_dist_rut<-acf(res.1000.10.rut_df$dist, n.sim=999, lag.max=20, na.action=na.pass, level=95, plot=FALSE)

jpeg("./acf_dist.rut.jpg")
plot(acf_dist_rut, main="Fall Rut")
dev.off()


acf_dist_lwin<-acf(res.1000.10.lwin_df$dist, n.sim=999, lag.max=20, na.action=na.pass, level=95, plot=FALSE)

jpeg("./acf_dist.lwin.jpg")
plot(acf_dist_lwin, main="Late Winter")
dev.off()

acf_dist_ewin<-acf(res.1000.10.ewin_df$dist, n.sim=999, lag.max=20, na.action=na.pass, level=95, plot=FALSE)

jpeg("./acf_dist.ewin.jpg")
plot(acf_dist_ewin, main="Early Winter")
dev.off()

#Plot distance ACF
#does not accept non-linear argument 

jpeg("./acfdist_2019_RT.jpg")
acfdist.ltraj(res.1000.10, "RT.100", lag=20)
dev.off()

jpeg("./acfdist_calving_RT.jpg")
acfdist.ltraj(res.1000.10.calv, "RT.100", lag=20)
dev.off()

jpeg("./acfdist_summer_RT.jpg")
acfdist.ltraj(res.1000.10.sum, "RT.100", lag=20)
dev.off()

jpeg("./acfdist_rut_RT.jpg")
acfdist.ltraj(res.1000.10.rut, "RT.100", lag=20)
dev.off()

jpeg("./acfdist_e.wint_RT.jpg")
acfdist.ltraj(res.1000.10.ewin, "RT.100", lag=20)
dev.off()

jpeg("./acfdist_l.wint_RT.jpg")
acfdist.ltraj(res.1000.10.lwin, "RT.100", lag=20)
dev.off()
#NOTE: cyclical pattern of temporal autocorrelation (~5 lags= 1 day)




#===Calculate Global and local Moran's I for residence time====#

#Convert df to Spatial Points Data Frame to begin Moran's I Calculation 
coords=res.1000.10_df[c("x", "y")]
res.1000.10_sp<-SpatialPointsDataFrame(coords=coords, data=res.1000.10_df, proj4string=crs)

coords=res.1000.10.calv_df[c("x", "y")]
res.1000.10.calv_sp<-SpatialPointsDataFrame(coords=coords, data=res.1000.10.calv_df, proj4string=crs)

coords=res.1000.10.sum_df[c("x", "y")]
res.1000.10.sum_sp<-SpatialPointsDataFrame(coords=coords, data=res.1000.10.sum_df, proj4string=crs)

coords=res.1000.10.rut_df[c("x", "y")]
res.1000.10.rut_sp<-SpatialPointsDataFrame(coords=coords, data=res.1000.10.rut_df, proj4string=crs)

coords=res.1000.10.ewin_df[c("x", "y")]
res.1000.10.ewin_sp<-SpatialPointsDataFrame(coords=coords, data=res.1000.10.ewin_df, proj4string=crs)

coords=res.1000.10.lwin_df[c("x", "y")]
res.1000.10.lwin_sp<-SpatialPointsDataFrame(coords=coords, data=res.1000.10.lwin_df, proj4string=crs)

plot(res.1000.10.ewin_sp$y, res.1000.10.ewin_sp$x)

#Calculate maximum extent of movement



t=res.1000.10.ewin_sp
title="Early_Winter"

#Adjust maximum spatial extent for each time span.

mcp<-mcp(t, percent=100)
mc_b<-gBuffer(mcp, width = 100)
area(mcp)
area(mc_b)

hr<-mc_b

plot(hr)
#Establish neighbourhoods
th  <-as(dirichlet(as.ppp(t)), "SpatialPolygons")
proj4string(th) <- proj4string(t)

th.z     <- over(th, t)
th.spdf  <-  SpatialPolygonsDataFrame(th, th.z)
th.spdf<-th.spdf[!is.na(th.spdf$dist),]
proj4string(th.spdf)<-crs

#trim surface to maximum extent of movement
hr   <- raster::intersect(mcp_b,th.spdf)

plot(hr$y, hr$x)

#palette_explorer()
#Plot values
map_var<-tm_shape(hr)+
  tm_polygons(col = "dist",
              title="Distance (m)",
              style = "jenks",
              palette = "Blues", n=6, border.alpha = 0)+
  tm_layout(title=title, legend.outside=TRUE, legend.outside.position = c("right", "bottom"))+
  tm_scale_bar()

png(paste("./Dist_", title ,".png", sep=""))
map_var
dev.off()

#Assign neighbours based on Queen's Case
nb<-poly2nb(hr, queen=FALSE)

#Create neighbourhood network
net<-nb2lines(nb, coords=coordinates(hr))
proj4string(net)<-proj4string(hr)

plot(net)



#list of lagged weights
lw<-nb2listw(nb, zero.policy=TRUE, style="W")

print.listw(lw, zero.policy = TRUE)


#create lagged means
hr$Lag_Mean= lag.listw(lw, hr$dist, zero.policy=TRUE)

#map lagged means
#palette_explorer()

map_lag_Mean<-tm_shape(hr)+
  tm_polygons(col="Lag_Mean",
              main.title="Lagged Means of Step Distance",
              style="jenks",
              palette="OrRd", n=6, border.alpha=0)+
  tm_layout(title=title, legend.outside=TRUE)+
  tm_scale_bar()

jpeg(paste("./lag.m_", title, ".jpeg", sep=""))
map_lag_Mean
dev.off()


#===Moran range function to return range in values to compare===#
moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}

#================================================================#

moran.range(lw)
#-0.5595129  1.0217957

#mi<-moran.test(th.spdf$RT.1000, lw, zero.policy=TRUE)
#weights: lw    

#Note: default moran.test() requires independent data. 
#Moran I statistic standard deviate = 58.738, p-value < 2.2e-16
#alternative hypothesis: greater
#sample estimates:
#  Moran I statistic       Expectation          Variance 
#       0.8499273519     -0.0005977286      0.0002096664 

#mI<-mi$estimate[[1]]
#mIe<-mi$estimate[[2]]
#mIvar<-mi$estimate[[3]]
#mIp<-mi$p.value

#z=(mI-mIe)/sqrt(mIvar)
#z=58.73848
#p= 2*pnorm(-abs(z))
#0

#Monte Carlo permetation significance test 


mi_mc<-moran.mc(hr$dist, lw, nsim= 999, na.action(na.omit), zero.policy=TRUE)
mi_mc
#Monte-Carlo simulation of Moran I
#data:  th.spdf$RT.1000 
#weights: lw  
#number of simulations + 1: 1000 

#statistic = 0.84993, observed rank = 1000, p-value = 0.001
#alternative hypothesis: greater


#data:  th.spdf$RT.1000 
#weights: lw  
#number of simulations + 1: 100 
#
#statistic = 0.84993, observed rank = 100, p-value = 0.01
#alternative hypothesis: greater

#No difference between 100 and 1000 repetitions for statistic, but decrease in p value with incresed repetitions. 


#======Calculate local moran's i=========#

#note that a monte carlo significance permetation method is not available. Z-test and p-value likley inflated. 
lisa.test<-localmoran(hr$dist, lw)

hr$Ii<-lisa.test[,1]
hr$E.Ii<-lisa.test[,2]
hr$Var.Ii<-lisa.test[,3]
hr$Z.Ii<-lisa.test[,4]
hr$P<-lisa.test[,5]

map_LISA<-tm_shape(hr)+
  tm_polygons(col="Ii",
              title="Local Moran's I",
              style="jenks",
              palette="YlGnBu", n=6, border.alpha=0)+
  tm_layout(title=title, legend.outside=TRUE)+
  tm_scale_bar()

map_LISA

png(paste("./LocalM_dist", title ,".png", sep=""))
map_LISA
dev.off()

#Plot local moran's i scatter plot

png(paste("./LocalM_dist",title,"_scatter.png"))
moran.plot(hr$dist,lw, zero.policy=TRUE, spChk=NULL, labels=NULL,
           xlab="SD", ylab="Spatially Lagged SD", quiet=NULL)
title(main=title)
dev.off()
#=========CODE STOPP==============#


