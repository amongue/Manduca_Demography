library(raster)
library(sp)
#get worldclim data
all_clim <- getData("worldclim",var="bio",res=10)
#temp and precip
r <- all_clim[[c(1,12)]]
names(r) <- c("Temp","Prec")

#let's find our lats and lons

#rocky mount, nc is 35.892667,-77.677166
#madera canyon, az is 31.723079, -110.880087
#williamstown,ks is 39.092097, -95.334951

lats <- c(35.892667 , 31.723079, 39.092097)
lons <- c(-77.677166, -110.880087, -95.334951) 

coords <- data.frame(x=lons,y=lats)

points <- SpatialPoints(coords, proj4string = r@crs)

values <- extract(r,points)

df <- cbind.data.frame(coordinates(points),values)

#precip
plot(r[[2]], ylim=c(10,60),xlim=c(-130,-50), main="mean precip")
plot(points,add=T)

#####temp stuff####
temps <- all_clim[[1:9]]
names(temps) <- c("Mean.Temp.Annual","Mean.Diurnal.Rng","Isotherm","Seasonality","Max.Hot","Min.Cold","Temp.Rng.Annual","Temp.Wet", "Temp.Dry")

values <- extract(temps,points)

df_temp <- cbind.data.frame(coordinates(points),values)

par(mfrow=c(1,1))
plot(temps[[4]], ylim=c(10,60),xlim=c(-130,-50), main="",xaxt='n',yaxt='n')
plot(points,add=T,cex=2, pch=16)
mtext("Climatic variation across sampling range", side=3, line=1, cex=2)
text(lons,lats + 2,c("n = 12", "n = 8", "n = 4"),cex=1.2)
mtext("Seasonal temperature variation", side=4, line=0, cex=1.05, las=0)
text(lons,lats + 4,c("NC", "AZ", "KS"),cex=1.2)


