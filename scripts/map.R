library(maps)
library(mapdata)
map('worldHires','Australia')

library(mapproj)

map <- data.frame(
  
  "category" <- c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3),
  
  "name" <- c("Snowy Creek",
              "Gibbo River",
              "Nariel Creek",
              "Goodradigbee River",
              "Jacobs River",
              "Tuross River at Belowra",
              "Genoa River",
              "Wallagaraugh River",
              "Mann River",
              "Cataract Creek",
              "Jilliby Creek",
              "Sportsmans Creek",
              "Mammy Johnsons River",
              "Wadbilliga River",
              "Tuross River D/S Wadbilliga Junction"),

  "long" <-  c(147.413,
              147.709,
              147.826,
              148.731,
              148.427,
              149.709,            
              149.321,
              149.714,
              152.105,
              152.217,
              151.389,
              152.981,
              151.979,
              149.694,
              149.761),
  
  "lat"  <-  c(-36.569,
              -36.756,
              -36.444,            
              -35.421,            
              -36.727,
              -36.201,           
              -37.174,
              -37.371,
              -29.695,
              -28.934,
              -33.246,
              -29.467,           
              -32.244,
              -36.259,
              -36.197),
  
  stringsAsFactors = FALSE
  
)          
            
colnames(map) <- c("category","name","long","lat")            
          

map(database= "worldHires", 
    regions="Australia",
    xlim=c(135,155), 
    ylim=c(-45,-25),
    col="gray98",
    fill=TRUE,
    resolution=0)
points((subset(map, category == 1))$long, 
       (subset(map, category == 1))$lat, 
       pch=21, 
       cex=1,
       bg="gray40")  
points((subset(map, category == 2))$long, 
       (subset(map, category == 2))$lat, 
       pch=22, 
       cex=1,
       bg="gray60")  
points((subset(map, category == 3))$long, 
       (subset(map, category == 3))$lat, 
       pch=23, 
       cex=1,
       bg="gray80")  

map("rivers", add=TRUE)






ogrInfo("data/drainage divisions", "rbasin_polygon")
ogrListLayers("data/drainage divisions/rbasin_polygon.shp")
drainageMap <- readOGR("data/drainage divisions", "rbasin_polygon")

attributes <- c("A

