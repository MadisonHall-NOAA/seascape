
#====================================================================================================
#================ Get VAST data
#====================================================================================================

# Put in your AFSC username/password here
username=""
password=""

# Define path
path<-getwd()

# Enter species code (example here for Pacific ocean perch)
species<-30060

# Enter survey region (example here for Gulf of Alaska)
region<-"'GOA'"

# Enter survey years for GOA (have this part hardwired for alck of better idea at the moment)
srv_yrs<-c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015,2017,2019)
srv_nms<-c("1999 GULF OF ALASKA TRIENNIAL SURVEY","2001 Gulf of Alaska Groundfish Survey","2005 Gulf of Alaska Bottom Trawl Survey","2007 Gulf of Alaska Bottom Trawl Survey","2009 Gulf of Alaska Bottom Trawl Survey","2011 Gulf of Alaska Bottom Trawl Survey","2013 Gulf of Alaska Bottom Trawl Survey","2015 Gulf of Alaska Bottom Trawl Survey","2017 Gulf of Alaska Bottom Trawl Survey","2019 Gulf of Alaska Bottom Trawl Survey","GOA TRIENNIAL SURVEY","Gulf of Alaska Biennial Survey","GULF OF ALASKA TRIENNIAL SURVEY")

# Call RODBC library and connect to AFSC database
library(RODBC)
channel=odbcConnect("afsc",uid=username,pwd=password,believeNRows=FALSE)

# Get haul/catch/cruise data
haul_data<-sqlQuery(channel,paste0("SELECT * FROM RACEBASE.HAUL WHERE ((RACEBASE.HAUL.REGION)=",region,")"),believeNRows=FALSE)
catch_data<-sqlQuery(channel,paste0("SELECT * FROM RACEBASE.CATCH WHERE ((RACEBASE.CATCH.SPECIES_CODE)=",species," AND (RACEBASE.CATCH.REGION)=",region,")"),believeNRows=FALSE)
cruise_data<-sqlQuery(channel,paste0("SELECT * FROM RACEBASE.CRUISE WHERE ((RACEBASE.CRUISE.REGION)=",region,")"),believeNRows=FALSE)

# Define which surveys were GOA surveys after 1990
cruise_data<-cruise_data[which(cruise_data$SURVEY_NAME %in% srv_nms),]
cruise_data<-cruise_data[which(as.numeric(substr(cruise_data$CRUISE,1,4)) %in% srv_yrs),]

# Now sort out haul and catch data that are associatesd with GOA surveys
haul_data<-haul_data[which(haul_data$CRUISE %in% cruise_data$CRUISE),]
catch_data<-catch_data[which(catch_data$CRUISE %in% cruise_data$CRUISE),]

# Add POP catch (and other info) to haul data (basically, make sure we have hauls with 0 catch in data along with location,etc)
YEAR<-as.numeric(substr(haul_data$CRUISE,1,4))
SPECIES_CODE<-rep(species,length(YEAR))
WEIGHT<-vector(mode="numeric",length=length(haul_data$HAULJOIN))
WEIGHT[which(haul_data$HAULJOIN %in% catch_data$HAULJOIN)]<-catch_data$WEIGHT
NUMBER<-vector(mode="numeric",length=length(haul_data$HAULJOIN))
NUMBER[which(haul_data$HAULJOIN %in% catch_data$HAULJOIN)]<-catch_data$NUMBER

# Define data for GOA POP and remove hauls not used for abundance est
POP_VAST_data<-cbind(haul_data,YEAR,SPECIES_CODE,WEIGHT,NUMBER)
POP_VAST_data_Y<-subset(POP_VAST_data,POP_VAST_data$ABUNDANCE_HAUL=="Y")

# Write out data
write.csv(POP_VAST_data,paste(path,"/POP_VAST_data.csv",sep=""))
write.csv(POP_VAST_data_Y,paste(path,"/POP_VAST_data_Y.csv",sep=""))
