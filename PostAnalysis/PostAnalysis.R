#-----------#
#-Libraries-#
#-----------#

library(tidyverse)
library(reshape2)
library(abind)
library(ggthemes)
library(inlabru)
library(RColorBrewer)
library(grid)
library(extrafont)
loadfonts(quiet = TRUE)

#--------------------#
#-Simulation example-#
#--------------------#

source('~/Monarchs/PostAnalysis/simfn.R', echo=TRUE)

attach(simfn())

#----------------#
#-Compile output-#
#----------------#

filenames <- list.files(path = "Z:/Monarch/DataAnalysis/Simulation/Output", pattern = "output", full.names = TRUE)

load(filenames[1])

est <- out[[1]]

convergence <- out[[2]]

for(i in 2:length(filenames)){
  load(filenames[i])
  est <- abind(est, out[[1]], along = 3)
  convergence <- abind(convergence, out[[2]], along = 3)
}

#Remove non-converged iterations (n = 10)
nonconverged <- which(convergence[,,] != TRUE, arr.ind = TRUE)[,3]

if(length(nonconverged) != 0){
  est <- est[,,-nonconverged]
}

est <- est[,,1:1000]

#------------------------#
#-Root mean square error-#
#------------------------#

MSE <- apply((est[,c(2,4),] - est[,rep(1,2),])^2, MARGIN = c(1,2), FUN = mean, na.rm = TRUE)

RMSE <- sqrt(MSE)

#---------------#
#-Relative bias-#
#---------------#

y75 <- apply(est[,c(2,4),] - est[,rep(1,2),], MARGIN = c(1,2), FUN = quantile, probs = 0.75, na.rm = TRUE)
y50 <- apply(est[,c(2,4),] - est[,rep(1,2),], MARGIN = c(1,2), FUN = quantile, probs = 0.5, na.rm = TRUE)
y25 <- apply(est[,c(2,4),] - est[,rep(1,2),], MARGIN = c(1,2), FUN = quantile, probs = 0.25, na.rm = TRUE)
ymax <-  ((y75 - y25) * 1.5) + y75
ymin <- y25 - ((y75 - y25) * 1.5)

rb75 <- apply((est[,c(2,4),] - est[,rep(1,2),])/est[,rep(1,2),], MARGIN = c(1,2), FUN = quantile, probs = 0.75, na.rm = TRUE)
rb50 <- apply((est[,c(2,4),] - est[,rep(1,2),])/est[,rep(1,2),], MARGIN = c(1,2), FUN = quantile, probs = 0.5, na.rm = TRUE)
rb25 <- apply((est[,c(2,4),] - est[,rep(1,2),])/est[,rep(1,2),], MARGIN = c(1,2), FUN = quantile, probs = 0.25, na.rm = TRUE)
rbmax <-  ((rb75 - rb25) * 1.5) + rb75
rbmin <- rb25 - ((rb75 - rb25) * 1.5)

data <- rbind(ymin, y25, y50, y75, ymax)
params <- as.factor(rownames(data))
data <- data.frame(data)
rownames(data) <- NULL
data$params <- params
data$quantile <- rep(c("qmin","q25","q50","q75","qmax"), each = 9)
data <- melt(data, id = c("params", "quantile"))
data <- data %>% rename(type = variable) %>%
        mutate(params = recode(params, "N[1]" = "N1", "N[2]" = "N2"),
               type = recode(type, "IM.Estimate" = "Integrated", "PO.Estimate" = "Presence-only"))
data <- data %>% pivot_wider(names_from = quantile, values_from = value)


rb.data <- rbind(rbmin, rb25, rb50, rb75, rbmax)
params <- as.factor(rownames(rb.data))
rb.data <- data.frame(rb.data)
rownames(rb.data) <- NULL
rb.data$params <- params
rb.data$quantile <- rep(c("qmin","q25","q50","q75","qmax"), each = 9)
rb.data <- melt(rb.data, id = c("params", "quantile"))
rb.data <- rb.data %>% rename(type = variable) %>%
  mutate(params = recode(params, "N[1]" = "N1", "N[2]" = "N2"),
         type = recode(type, "IM.Estimate" = "Integrated", "PO.Estimate" = "Presence-only"))
rb.data <- rb.data %>% pivot_wider(names_from = quantile, values_from = value)

cbind(rb.data[,1:2], round(rb.data[,3:7] * 100, digits = 2))

#----------#
#-Figure 1-#
#----------#

Fig1A <- ggplotGrob(ggplot() +
  gg(intensity1) +
  scale_fill_gradientn(limit = c(0,2), colors = brewer.pal(n = 9, name =  "YlGnBu")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_sf(data = pop1, fill = NA, col = "Black", size = 0.25, linetype = 0) +
  # geom_sf(data = st_buffer(domain, 0.05, endCapStyle = "SQUARE"), fill = NA, col = "Black", size = 1, linetype = "dashed") +  
  theme_few() +
  theme(plot.margin = unit(c(-0.25, -0.25, 0, -0.25), "in"),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(0, "in")))

Fig1B <- ggplotGrob(ggplot() +
  gg(intensity2) +
  scale_fill_gradientn(limit = c(0,2), colors = brewer.pal(n = 9, name =  "YlGnBu")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_sf(data = pop2, fill = NA, col = "Black", size = 0.25, linetype = 0) +
  # geom_sf(data = st_buffer(domain2, 0.05, endCapStyle = "SQUARE"), fill = NA, col = "Black", size = 1, linetype = "dashed") +  
  theme_few() +
  theme(plot.margin = unit(c(-0.25, -0.25, 0, -0.25), "in"),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(0, "in")))

Fig1C <- ggplotGrob(ggplot() +
  gg(intensity1) +
  scale_fill_gradientn(limit = c(0,2), colors = brewer.pal(n = 9, name =  "YlGnBu")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_sf(data = po1, fill = NA, col = "White", size = 0.75, linetype = 0) +
  # geom_sf(data = st_buffer(domain1, 0.05, endCapStyle = "SQUARE"), fill = NA, col = "Black", size = 1, linetype = "dashed") +  
  geom_sf(data = sites, fill = NA, col = "Black", size = 0.25, linetype = 2) +
  geom_text(data = countdf, aes(Var1, Var2, label = count), color = "Black") +
  theme_few() +
  theme(plot.margin = unit(c(-0.25, -0.25, 0, -0.25), "in"),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(size = 14),
        panel.spacing = unit(0, "in")))

Fig1D <- ggplotGrob(ggplot() +
  gg(intensity2) +
  scale_fill_gradientn(limit = c(0,2), colors = brewer.pal(n = 9, name =  "YlGnBu")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_sf(data = po2, fill = NA, col = "White", size = 0.75, linetype = 0) +
  # geom_sf(data = st_buffer(domain2, 0.05, endCapStyle = "SQUARE"), fill = NA, col = "Black", size = 1, linetype = "dashed") +  
  theme_few() +
  theme(plot.margin = unit(c(-0.25, -0.25, 0, -0.25), "in"),
        legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(0, "in")))

Legend <- ggplotGrob(ggplot() +
          gg(intensity2) +
          scale_fill_gradientn(limit = c(0,2), colors = brewer.pal(n = 9, name =  "YlGnBu"), name = "Expected population density",
                               guide = guide_colorbar(ticks.colour = "Black",
                                                      frame.colour = "Black",
                                                      barwidth = 20,
                                                      barheight = 0.75,
                                                      title.position = "bottom")) +
          scale_x_continuous(expand = c(0,0)) +
          scale_y_continuous(expand = c(0,0)) +
          theme_few() +
          theme(text = element_text(size = 14),
                legend.title = element_text(hjust = 0.5),
                legend.position = "bottom",
                legend.box.margin = unit(c(0.25, 0, 0, 0), "in"),
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                panel.spacing = unit(0, "in")))$grob[[15]]


Fig1E <- ggplotGrob(rb.data %>% filter(params %in% c("beta0", "beta1", "delta1")) %>%
  ggplot(., aes(x = type)) +
  geom_boxplot(aes(ymin = qmin, lower = q25, middle = q50, upper = q75, ymax = qmax, fill = type), 
               stat = "identity", size = 0.75) +
  facet_wrap(. ~ params, scales = 'free_y',) + 
  geom_hline(yintercept = 0, col = "black", size = 1) +
  scale_fill_manual(values = c("#2C7FB8", "grey")) +
  theme_few() +
  theme(plot.margin = unit(c(-0.1, 0.05, -0.1, 0), "in"),
        axis.title.y = element_text(margin = margin(0,0,0,0)),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(color = "white"),
        legend.position = "none") +
  labs(y = "", x = ""))

param.lab <- c(expression(lambda[0]), expression(beta), expression(delta))

loc <- Fig1E$layout[grep("strip", Fig1E$layout$name),1:4]
loc <- loc %>% arrange(t,l)

for(i in 1:3){
  Fig1E <- gtable::gtable_add_grob(Fig1E, textGrob(param.lab[i], gp = gpar(fontsize = 14, fontface = 2)),
                                   t = loc$t[i], l = loc$l[i])
}


Fig1F <- ggplotGrob(rb.data %>% filter(params %in% c("N1", "N2")) %>%
  ggplot(., aes(x = params)) +
  geom_boxplot(aes(ymin = qmin, lower = q25, middle = q50, upper = q75, ymax = qmax, fill = type), 
               stat = "identity", size = 0.75) +
  facet_wrap(. ~ type, scales = 'free_y') + 
  geom_hline(yintercept = 0, col = "black", size = 1) +
  scale_fill_manual(values = c("#2C7FB8", "grey")) +
  theme_few() +
  theme(plot.margin = unit(c(-0.1, 0.05, -0.1, 0), "in"),
        axis.title.y = element_text(margin = margin(0,0,0,0)),
        # axis.ticks.x = element_blank(),
        # axis.text.x = element_blank(),
        legend.position = "none",
        strip.text = element_text(color = "white")) +
  labs(y = "", x = ""))

Fig1F$heights <- Fig1E$heights

Legend2 <- ggplotGrob(data %>% filter(params == "beta0") %>%
                        ggplot(., aes(x = type)) +
                        geom_boxplot(aes(ymin = qmin, lower = q25, middle = q50, upper = q75, ymax = qmax, fill = type), 
                                     stat = "identity", size = 0.75) +
                        scale_fill_manual(values = c("#2C7FB8", "grey")) +
                        theme_few() +
                        theme(text = element_text(size = 18),
                              plot.margin = unit(c(0, 0, 0, 0), "in"),
                              legend.position = "bottom",
                              legend.title = element_blank()) +
                        labs(y = "", x = ""))$grob[[15]]

#Add letters to figures
Figure1A <- arrangeGrob(Fig1A, top = grid::textGrob("A", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0, hjust = 0,
                                                    gp=grid::gpar(fontsize=18, fontface = 2)))
Figure1B <- arrangeGrob(Fig1B, top = grid::textGrob("B", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0, hjust = 0,
                                                    gp=grid::gpar(fontsize=18, fontface = 2)))
Figure1C <- arrangeGrob(Fig1C, top = grid::textGrob("C", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0, hjust = 0,
                                                    gp=grid::gpar(fontsize=18, fontface = 2)))
Figure1D <- arrangeGrob(Fig1D, top = grid::textGrob("D", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0, hjust = 0,
                                                    gp=grid::gpar(fontsize=18, fontface = 2)))
Figure1E <- arrangeGrob(Fig1E, top = grid::textGrob("E", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0.5, hjust = 0,
                                                    gp=grid::gpar(fontsize=18, fontface = 2)))
Figure1F <- arrangeGrob(Fig1F, top = grid::textGrob("F", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0.5, hjust = 0,
                                                    gp=grid::gpar(fontsize=18, fontface = 2)))

#Save Figure 1
tiff(file = "~/Monarchs/PostAnalysis/Figure1.tiff", res = 600, width = 6.5, height = 9, units = "in")
grid.arrange(arrangeGrob(Figure1A, Figure1B, Figure1C, Figure1D, nrow = 2, padding = unit(0, "in")),
             Legend,
             arrangeGrob(Figure1E, Figure1F, nrow = 1, padding = unit(0, "in")),
             Legend2,
             layout_matrix = matrix(rep(c(1,1,1,1,1,1,1,2,3,3,3,4), 2), nrow = 12),
             padding = unit(0, "in"))
dev.off()

ex <- est[7:8,c(1,2,4),2]
ex <- as.data.frame(melt(ex))
ex$SE <- c(NA, NA, melt(est[7:8,c(3,5),2])$value)
ex$CI <- ex$SE*1.96
ex$type <- rep(2:1,c(4,2))
             
ggplot(data = ex %>% filter(!grepl("Error", Var2)), aes(x=Var1, y = value, col = Var2)) +
  geom_point(size = 2) +
  geom_errorbar(aes(x=Var1, ymin = value - CI, ymax = value + CI), width = 0.1, size = 1) +
  facet_wrap(. ~ type, scales = 'free_y', dir = "v") +
  labs(y = "Abundance", x = "Time period") +
  theme_few() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_blank(),
    legend.title = element_blank()
  )

#-Figure 2-#

#domain_sp; domain_su; mesh[[]];

domain_sp <- st_read("~/Monarchs/DataFormatting/BaseData/Domain/Domain_Spring.shp")
domain_su <- st_read("~/Monarchs/DataFormatting/BaseData/Domain/Domain_Summer.shp")

prj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"

domain_sp <- st_transform(domain_sp, crs = st_crs(prj))
domain_su <- st_transform(domain_su, crs = st_crs(prj))

domain_sp_seg <- inla.sp2segment(as(domain_sp, "Spatial"))
domain_su_seg <- inla.sp2segment(as(domain_su, "Spatial"))

mesh <- list()

mesh[[1]] <- mesh[[2]] <- inla.mesh.2d(boundary=domain_sp_seg,
                                       max.edge=75,
                                       cutoff=30)

mesh[[3]] <- inla.mesh.2d(boundary=domain_su_seg,
                          max.edge=75,
                          cutoff=30)


sp.intensity <- raster::raster(domain_sp, resolution = c(10, 10))
sp.Grid <- coordinates(sp.intensity)
sp.GridA <- inla.spde.make.A(mesh[[1]], loc = sp.Grid)

su.intensity <- raster::raster(domain_su, resolution = c(10, 10))
su.Grid <- coordinates(su.intensity)
su.GridA <- inla.spde.make.A(mesh[[3]], loc = su.Grid)

df <- data.frame(cbind(t(rep(NA, 5))))
colnames(df) <- c("value", "x", "y", "period", "year")

#NDVI

NDVIspring <- raster::stack("~/Monarchs/DataFormatting/NDVIspring.tif")
NDVIsummer <- raster::stack("~/Monarchs/DataFormatting/NDVIsummer.tif")

NDVIspring <- raster::projectRaster(NDVIspring, sp.intensity, method = 'ngb')
NDVIsummer <- raster::projectRaster(NDVIsummer, su.intensity, method = 'ngb')

sp.NDVI <- matrix(NA, nrow = dim(sp.Grid)[1], ncol = t_n*2)
su.NDVI <- matrix(NA, nrow = dim(su.Grid)[1], ncol = t_n)

sp.NDVI[,1] <- 0.0001*raster::values(raster::overlay(NDVIspring[[1]], NDVIspring[[2]], fun = mean))
sp.NDVI[,2] <- 0.0001*raster::values(raster::overlay(NDVIspring[[3]], NDVIspring[[4]], fun = mean))
sp.NDVI[,3] <- 0.0001*raster::values(raster::overlay(NDVIspring[[5]], NDVIspring[[6]], fun = mean))
sp.NDVI[,4] <- 0.0001*raster::values(raster::overlay(NDVIspring[[7]], NDVIspring[[8]], fun = mean))
sp.NDVI[,5] <- 0.0001*raster::values(raster::overlay(NDVIspring[[9]], NDVIspring[[10]], fun = mean))
sp.NDVI[,6] <- 0.0001*raster::values(raster::overlay(NDVIspring[[11]], NDVIspring[[12]], fun = mean))

su.NDVI[,1] <- 0.0001*raster::values(raster::overlay(NDVIsummer[[1]], NDVIsummer[[2]], NDVIsummer[[3]], fun = mean))
su.NDVI[,2] <- 0.0001*raster::values(raster::overlay(NDVIsummer[[4]], NDVIsummer[[5]], NDVIsummer[[6]], fun = mean))
su.NDVI[,3] <- 0.0001*raster::values(raster::overlay(NDVIsummer[[7]], NDVIsummer[[8]], NDVIsummer[[9]], fun = mean))

sp.NDVI.stan <- (sp.NDVI - mean(c(unlist(sp.NDVI), unlist(su.NDVI))))/sd(c(unlist(sp.NDVI), unlist(su.NDVI)))
su.NDVI.stan <- (su.NDVI - mean(c(unlist(sp.NDVI), unlist(su.NDVI))))/sd(c(unlist(sp.NDVI), unlist(su.NDVI)))

#GDD

load("Z:/Monarch/DataFormatting/Daymet_sp.Rds")
load("Z:/Monarch/DataFormatting/Daymet_su.Rds")

# library(httr)
# 
# prj.daymet <- "+init=epsg:4326"
# 
# sp.nodexy <-  st_coordinates(st_transform(st_as_sf(data.frame(NA, X = sp.Grid[,1], Y = sp.Grid[,2]), 
#                                                    coords = c("X", "Y"), crs = prj), crs = st_crs(prj.daymet)))
# sp.lat <- sp.nodexy[,2]
# sp.long <- sp.nodexy[,1]
# 
# su.nodexy <-  st_coordinates(st_transform(st_as_sf(data.frame(NA, X = su.Grid[,1], Y = su.Grid[,2]), 
#                                                    coords = c("X", "Y"), crs = prj), crs = st_crs(prj.daymet)))
# su.lat <- su.nodexy[,2]
# su.long <- su.nodexy[,1]
# 
# #Spring mesh dates
# sp.start_date <- c("2016-03-01", "2016-03-22",
#                    "2017-03-01", "2017-03-22",
#                    "2018-03-01", "2018-03-22")
# sp.end_date <- c("2016-04-04", "2016-05-02",
#                  "2017-04-04", "2017-05-02",
#                  "2018-04-04", "2018-05-02")
# 
# #Summer mesh dates
# su.start_date <- c("2016-04-26",
#                    "2017-04-26",
#                    "2018-04-26")
# su.end_date <- c("2016-06-06",
#                  "2017-06-06",
#                  "2018-06-06")
# 
# source("~/Monarchs/DataFormatting/GDD/gdd_fun.R")
# 
# daymet <- function(long,lat,start_date,end_date){
#   query <- list("lat" = lat,
#                 "lon" = long,
#                 # "vars" = "tmax,tmin,prcp",
#                 "vars" = "tmax,tmin",
#                 "start" = start_date,
#                 "end" = end_date)
#   GET(url = "https://daymet.ornl.gov/single-pixel/api/data?lat={}&lon={}", 
#       query = query, write_disk(path = "~/Monarchs/DataFormatting/GDD/tmp.GDD.csv", overwrite = TRUE))
#   return(read_gdd(file = "~/Monarchs/DataFormatting/GDD/tmp.GDD.csv", date = end_date, lat = lat, long = long))
# }
# 
# #Spring mesh
# daymet.sp <- NULL
# for(i in 1:6){
#   daymet.sp <- rbind(daymet.sp, do.call(rbind, mapply(FUN = daymet, long = sp.long, lat = sp.lat, start_date = sp.start_date[i], end_date = sp.end_date[i], SIMPLIFY = FALSE)))
# }
# 
# #Summer mesh
# daymet.su <- NULL
# for(i in 1:3){
#   daymet.su <- rbind(daymet.su, do.call(rbind, mapply(FUN = daymet, long = su.long, lat = su.lat, start_date = su.start_date[i], end_date = su.end_date[i], SIMPLIFY = FALSE)))
# }

daymet <- rbind(daymet.sp, daymet.su)
daymet <- daymet %>% mutate(Year = ifelse(grepl("2016", date), "2016",
                                          ifelse(grepl("2017", date), "2017", "2018")),
                            period = ifelse(grepl("04", date), 1,
                                            ifelse(grepl("05", date), 2, 3)))

daymet <- st_as_sf(daymet, coords = c("lon", "lat"))

for(i in 1:dim(daymet)[1]){
  if(!any(is.na(daymet[i,]))) next
  yrID <- daymet[i,] %>% select(Year) %>% .$Year %>% as.factor(.)
  pID <- daymet[i,] %>% select(period) %>% .$period %>% as.factor(.)
  tmp0 <- daymet %>% filter(Year == yrID & period == pID) %>% drop_na()
  tmp1 <- as.numeric(st_distance(daymet[i,], tmp0) %>% apply(., 1, order))
  tmp2 <- as.numeric(st_distance(daymet[i,], tmp0) %>% apply(., 1, sort))
  tmp3.gdd1 <- tmp0 %>% select(gdd1) %>% .$gdd1 %>% as.numeric(.)
  tmp3.gdd2 <- tmp0 %>% select(gdd2) %>% .$gdd2 %>% as.numeric(.)
  tmp4 <- (1/tmp2[1:3])/sum(1/tmp2[1:3])
  daymet[i,]$gdd1 <- tmp3.gdd1[tmp1[1:3]] %*% tmp4
  daymet[i,]$gdd2 <- tmp3.gdd2[tmp1[1:3]] %*% tmp4
}

daymet <- daymet %>% mutate(daysaccum = ifelse(period == 1, 35, 42),
                            gdd1.avg = gdd1/daysaccum * 14,
                            gdd2.avg = gdd2/daysaccum * 14)

daymet$gdd1.stan <- (daymet$gdd1.avg - mean(daymet$gdd1.avg))/sd(daymet$gdd1.avg)
daymet$gdd2.stan <- (daymet$gdd2.avg - mean(daymet$gdd2.avg))/sd(daymet$gdd2.avg)


daymet$NDVI <- NA

daymet[daymet$Year == "2016" & daymet$period == 1,"NDVI"] <- sp.NDVI[,1]
daymet[daymet$Year == "2016" & daymet$period == 2,"NDVI"] <- sp.NDVI[,2]
daymet[daymet$Year == "2016" & daymet$period == 3,"NDVI"] <- su.NDVI[,1]

daymet[daymet$Year == "2017" & daymet$period == 1,"NDVI"] <- sp.NDVI[,3]
daymet[daymet$Year == "2017" & daymet$period == 2,"NDVI"] <- sp.NDVI[,4]
daymet[daymet$Year == "2017" & daymet$period == 3,"NDVI"] <- su.NDVI[,2]

daymet[daymet$Year == "2018" & daymet$period == 1,"NDVI"] <- sp.NDVI[,5]
daymet[daymet$Year == "2018" & daymet$period == 2,"NDVI"] <- sp.NDVI[,6]
daymet[daymet$Year == "2018" & daymet$period == 3,"NDVI"] <- su.NDVI[,3]

daymet$omega <- NA

daymet[daymet$Year == "2016" & daymet$period == 1,"omega"] <- sp.GridA %*% omg_sp1
daymet[daymet$Year == "2016" & daymet$period == 2,"omega"] <- sp.GridA %*% omg_sp2
daymet[daymet$Year == "2016" & daymet$period == 3,"omega"] <- su.GridA %*% omg_su

daymet[daymet$Year == "2017" & daymet$period == 1,"omega"] <- sp.GridA %*% omg_sp1
daymet[daymet$Year == "2017" & daymet$period == 2,"omega"] <- sp.GridA %*% omg_sp2
daymet[daymet$Year == "2017" & daymet$period == 3,"omega"] <- su.GridA %*% omg_su

daymet[daymet$Year == "2018" & daymet$period == 1,"omega"] <- sp.GridA %*% omg_sp1
daymet[daymet$Year == "2018" & daymet$period == 2,"omega"] <- sp.GridA %*% omg_sp2
daymet[daymet$Year == "2018" & daymet$period == 3,"omega"] <- su.GridA %*% omg_su

#Covariate effects
daymet$beta <- NA

daymet[daymet$Year == "2016" & daymet$period == 1,"beta"] <- beta[1]
daymet[daymet$Year == "2016" & daymet$period == 2,"beta"] <- beta[2]
daymet[daymet$Year == "2016" & daymet$period == 3,"beta"] <- beta[3]

daymet[daymet$Year == "2017" & daymet$period == 1,"beta"] <- beta[4]
daymet[daymet$Year == "2017" & daymet$period == 2,"beta"] <- beta[5]
daymet[daymet$Year == "2017" & daymet$period == 3,"beta"] <- beta[6]

daymet[daymet$Year == "2018" & daymet$period == 1,"beta"] <- beta[7]
daymet[daymet$Year == "2018" & daymet$period == 2,"beta"] <- beta[8]
daymet[daymet$Year == "2018" & daymet$period == 3,"beta"] <- beta[9]

daymet$beta1 <- NA

daymet[daymet$period == 1,"beta1"] <- beta1[1]
daymet[daymet$period == 2,"beta1"] <- beta1[2]
daymet[daymet$period == 3,"beta1"] <- beta1[3]

daymet$beta2 <- NA

daymet[daymet$period == 1,"beta2"] <- beta2[1]
daymet[daymet$period == 2,"beta2"] <- beta2[2]
daymet[daymet$period == 3,"beta2"] <- beta2[3]

daymet$beta3 <- NA

daymet[daymet$period == 1,"beta3"] <- beta3[1]
daymet[daymet$period == 2,"beta3"] <- beta3[2]
daymet[daymet$period == 3,"beta3"] <- beta3[3]

daymet$beta4 <- NA

daymet[daymet$period == 1,"beta4"] <- beta4[1]
daymet[daymet$period == 2,"beta4"] <- beta4[2]
daymet[daymet$period == 3,"beta4"] <- beta4[3]

# Prediction
daymet <- daymet %>% mutate(pred.fix = beta + beta1 * NDVI + beta2 * NDVI * NDVI + beta3 * gdd2.stan + beta4 * gdd2.stan * gdd2.stan,
                            pred.ran = pred.fix + omega)

daymet <- daymet %>% mutate(pred.fix = exp(pred.fix),
                            pred.ran = exp(pred.ran))


# sp.Pred <- cbind(exp(beta[1] + beta1[1] * sp.NDVI[,1] + beta2[1] * sp.NDVI[,1] * sp.NDVI[,1] + sp.GridA %*% omg_sp1),
#                  exp(beta[2] + beta1[2] * sp.NDVI[,2] + beta2[2] * sp.NDVI[,2] * sp.NDVI[,2] + sp.GridA %*% omg_sp2),
#                  exp(beta[4] + beta1[1] * sp.NDVI[,3] + beta2[1] * sp.NDVI[,3] * sp.NDVI[,3] + sp.GridA %*% omg_sp1),
#                  exp(beta[5] + beta1[2] * sp.NDVI[,4] + beta2[2] * sp.NDVI[,4] * sp.NDVI[,4] + sp.GridA %*% omg_sp2),
#                  exp(beta[7] + beta1[1] * sp.NDVI[,5] + beta2[1] * sp.NDVI[,5] * sp.NDVI[,5] + sp.GridA %*% omg_sp1),
#                  exp(beta[8] + beta1[2] * sp.NDVI[,6] + beta2[2] * sp.NDVI[,6] * sp.NDVI[,6] + sp.GridA %*% omg_sp2))
# 
# su.Pred <- cbind(exp(beta[3] + beta1[3] * su.NDVI[,1] + beta2[3] * su.NDVI[,1] * su.NDVI[,1] + su.GridA %*% omg_su),
#                  exp(beta[6] + beta1[3] * su.NDVI[,2] + beta2[3] * su.NDVI[,2] * su.NDVI[,2] + su.GridA %*% omg_su),
#                  exp(beta[9] + beta1[3] * su.NDVI[,3] + beta2[3] * su.NDVI[,3] * su.NDVI[,3] + su.GridA %*% omg_su))


tp <- cbind(c(1,2,1,2,1,2,3,3,3), c(1,1,2,2,3,3,1,2,3))
# 
# for(t in 1:6){
#   raster::values(sp.intensity) <- as.vector(sp.Pred[,t])
#   sp.intensity <- as(sp.intensity, "SpatialPixelsDataFrame")
#   sp.intensity <- st_as_sf(sp.intensity)
#   sp.intensity <- sp.intensity[domain_sp,]
#   sp.intensity <- as.data.frame(as(sp.intensity, "Spatial"))
#   sp.intensity$period <- tp[t,1]
#   sp.intensity$year <- tp[t,2]
#   colnames(sp.intensity) <- c("value", "x", "y", "period", "year")
#   df <- rbind(df, sp.intensity)
#   rm(sp.intensity)
#   sp.intensity <- raster::raster(domain_sp, resolution = c(10, 10))
# }
# 
# for(t in 1:3){
#   raster::values(su.intensity) <- as.vector(su.Pred[,t])
#   su.intensity <- as(su.intensity, "SpatialPixelsDataFrame")
#   su.intensity <- st_as_sf(su.intensity)
#   su.intensity <- su.intensity[domain_su,]
#   su.intensity <- as.data.frame(as(su.intensity, "Spatial"))
#   su.intensity$period <- tp[t+6,1]
#   su.intensity$year <- tp[t+6,2]
#   colnames(su.intensity) <- c("value", "x", "y", "period", "year")
#   df <- rbind(df, su.intensity)
#   rm(su.intensity)
#   su.intensity <- raster::raster(domain_su, resolution = c(10, 10))
# }
# 
# df <- df[-1,]
# 
# df1 <- st_as_sf(df, coords = c("x", "y"), crs = prj)

#Figures
prj.daymet <- "+init=epsg:4326"
st_crs(daymet) <- prj.daymet
daymet <- st_transform(daymet, crs = prj)
daymet <- daymet[st_union(domain_sp, domain_su),]

pred <- daymet %>% select(Year, period, gdd1.avg, NDVI, pred.fix, pred.ran)
pre
pred$Year <- as.factor(pred$Year)
pred$period <- as.factor(pred$period)

pred <- raster::raster(domain_sp, resolution = c(10, 10))
raster::values(pred) <- as.vector(daymet %>% filter(Year == "2016" & period == 1) %>% select(pred.ran) %>% .$pred.ran)


pal2 <- colorRampPalette(RColorBrewer::brewer.pal(9, 'YlGn'))(1000)
pal3 <- colorRampPalette(RColorBrewer::brewer.pal(11, 'Spectral'))(1000)
pal4 <- colorRampPalette(RColorBrewer::brewer.pal(9, 'Greens'))(20)
pal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 9, name =  "YlGnBu"))(100)

pal6 <- c(colorRampPalette(c("#FFFFFF", RColorBrewer::brewer.pal(9, 'Greens')[1]))(4),
          colorRampPalette(RColorBrewer::brewer.pal(9, 'Greens')[1:2])(2)[-1],
          colorRampPalette(RColorBrewer::brewer.pal(9, 'Greens')[2:3])(2)[-1],
          colorRampPalette(RColorBrewer::brewer.pal(9, 'Greens')[3:4])(2)[-1],
          colorRampPalette(RColorBrewer::brewer.pal(9, 'Greens')[4:5])(2)[-1],
          colorRampPalette(RColorBrewer::brewer.pal(9, 'Greens')[5:6])(2)[-1],
          colorRampPalette(RColorBrewer::brewer.pal(9, 'Greens')[6:7])(2)[-1],
          colorRampPalette(RColorBrewer::brewer.pal(9, 'Greens')[7:8])(4)[-1],
          colorRampPalette(RColorBrewer::brewer.pal(9, 'Greens')[8:9])(8)[-1],
          colorRampPalette(c(RColorBrewer::brewer.pal(9, 'Greens')[9], "#000000"))(2)[-1])

rounder <- function(x){round(exp(x), digits = 3)}

inverse_sp <- st_difference(st_as_sfc(st_bbox(domain_sp) + c(-50, -50, 50, 50)), domain_sp)
inverse_su <- st_difference(st_as_sfc(st_bbox(domain_su) + c(-50, -50, 50, 50)), domain_su)

Legend2 <- ggplotGrob(ggplot() +
                        geom_sf(data = daymet %>% filter(Year == "2016" & period == 1), aes(color = pred.fix), shape = 15) +
                        geom_sf(data = domain_sp, fill = NA, size = 1) +
                        scale_color_gradientn(colors = pal4, limits = c(log(0.0001), log(0.075)), 
                                              oob = scales::squish, label = rounder, name = "") +
                        theme_void() +
                        theme(plot.margin = margin(0,0.25,0,0.25, "in"),
                              panel.background = element_blank(),
                              panel.grid = element_blank(),
                              axis.text = element_blank(),
                              axis.ticks = element_blank(),
                              legend.position = "left",
                              legend.key.width = unit(0.125, "in"),
                              legend.key.height = unit(0.75, "in")) +
                        guides(color = guide_colorbar(ticks.colour = "black", frame.colour = "black")))$grob[[14]]

Fig2A <- ggplotGrob(ggplot() +
                      geom_sf(data = daymet %>% filter(Year == "2016" & period == 1), aes(color = pred.fix), shape = 15) +
                      geom_sf(data = domain_sp, fill = NA, size = 1) +
                      geom_sf(data = inverse_sp, fill = "white", color = NA) + 
                      scale_color_gradientn(colors = pal4, name = expression("Population Density / 10 m"^2), limits = c(log(0.0001), log(0.075)), oob = scales::squish, label = rounder) +
                      theme_void() +
                      theme(plot.margin = margin(0,0,0,0, "in"),
                            panel.background = element_blank(),
                            panel.grid = element_blank(),
                            axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            legend.position = "none"))

# Fig2AT <- ggplotGrob(ggplot() +
#     geom_sf(data = domain_sp, fill = NA, size = 2) +
#     theme(plot.margin = margin(0,-0.5,0,-0.5, "in"),
#           panel.grid = element_blank(),
#           axis.text = element_blank(),
#           axis.ticks = element_blank()))

Fig2B <- ggplotGrob(ggplot() +
                      geom_sf(data = daymet %>% filter(Year == "2016" & period == 2), aes(color = pred.fix), shape = 15) +
                      geom_sf(data = domain_sp, fill = NA, size = 1) +
                      geom_sf(data = inverse_sp, fill = "white", color = NA) + 
                      scale_color_gradientn(colors = pal4, name = expression("Population Density / 10 m"^2), limits = c(log(0.0001), log(0.075)), oob = scales::squish, label = rounder) +
                      theme_void() +
                      theme(plot.margin = margin(0,0,0,0, "in"),
                            panel.background = element_blank(),
                            panel.grid = element_blank(),
                            axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            legend.position = "none"))

Fig2C <- ggplotGrob(ggplot() +
                      geom_sf(data = daymet %>% filter(Year == "2016" & period == 3), aes(color = pred.fix), shape = 15) +
                      geom_sf(data = domain_su, fill = NA, size = 1) +
                      geom_sf(data = inverse_su, fill = "white", color = NA) + 
                      scale_color_gradientn(colors = pal4, name = expression("Population Density / 10 m"^2), limits = c(log(0.0001), log(0.075)), oob = scales::squish, label = rounder) +
                      theme_void() +
                      theme(plot.margin = margin(0,0,0,0, "in"),
                            panel.background = element_blank(),
                            panel.grid = element_blank(),
                            axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            legend.position = "none"))

Fig2D <- ggplotGrob(ggplot() +
                      geom_sf(data = daymet %>% filter(Year == "2017" & period == 1), aes(color = pred.fix), shape = 15) +
                      geom_sf(data = domain_sp, fill = NA, size = 1) +
                      geom_sf(data = inverse_sp, fill = "white", color = NA) + 
                      scale_color_gradientn(colors = pal4, name = expression("Population Density / 10 m"^2), limits = c(log(0.0001), log(0.075)), oob = scales::squish, label = rounder) +
                      theme_void() +
                      theme(plot.margin = margin(0,0,0,0, "in"),
                            panel.background = element_blank(),
                            panel.grid = element_blank(),
                            axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            legend.position = "none"))

Fig2E <- ggplotGrob(ggplot() +
                      geom_sf(data = daymet %>% filter(Year == "2017" & period == 2), aes(color = pred.fix), shape = 15) +
                      geom_sf(data = domain_sp, fill = NA, size = 1) +
                      geom_sf(data = inverse_sp, fill = "white", color = NA) + 
                      scale_color_gradientn(colors = pal4, name = expression("Population Density / 10 m"^2), limits = c(log(0.0001), log(0.075)), oob = scales::squish, label = rounder) +
                      theme_void() +
                      theme(plot.margin = margin(0,0,0,0, "in"),
                            panel.background = element_blank(),
                            panel.grid = element_blank(),
                            axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            legend.position = "none"))

Fig2F <- ggplotGrob(ggplot() +
                      geom_sf(data = daymet %>% filter(Year == "2017" & period == 3), aes(color = pred.fix), shape = 15) +
                      geom_sf(data = domain_su, fill = NA, size = 1) +
                      geom_sf(data = inverse_su, fill = "white", color = NA) + 
                      scale_color_gradientn(colors = pal4, name = expression("Population Density / 10 m"^2), limits = c(log(0.0001), log(0.075)), oob = scales::squish, label = rounder) +
                      theme_void() +
                      theme(plot.margin = margin(0,0,0,0, "in"),
                            panel.background = element_blank(),
                            panel.grid = element_blank(),
                            axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            legend.position = "none"))

Fig2G <- ggplotGrob(ggplot() +
                      geom_sf(data = daymet %>% filter(Year == "2018" & period == 1), aes(color = pred.fix), shape = 15) +
                      geom_sf(data = domain_sp, fill = NA, size = 1) +
                      geom_sf(data = inverse_sp, fill = "white", color = NA) + 
                      scale_color_gradientn(colors = pal4, name = expression("Population Density / 10 m"^2), limits = c(log(0.0001), log(0.075)), oob = scales::squish, label = rounder) +
                      theme_void() +
                      theme(plot.margin = margin(0,0,0,0, "in"),
                            panel.background = element_blank(),
                            panel.grid = element_blank(),
                            axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            legend.position = "none"))

Fig2H <- ggplotGrob(ggplot() +
                      geom_sf(data = daymet %>% filter(Year == "2018" & period == 2), aes(color = pred.fix), shape = 15) +
                      geom_sf(data = domain_sp, fill = NA, size = 1) +
                      geom_sf(data = inverse_sp, fill = "white", color = NA) + 
                      scale_color_gradientn(colors = pal4, name = expression("Population Density / 10 m"^2), limits = c(log(0.0001), log(0.075)), oob = scales::squish, label = rounder) +
                      theme_void() +
                      theme(plot.margin = margin(0,0,0,0, "in"),
                            panel.background = element_blank(),
                            panel.grid = element_blank(),
                            axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            legend.position = "none"))

Fig2I <- ggplotGrob(ggplot() +
                      geom_sf(data = daymet %>% filter(Year == "2018" & period == 3), aes(color = pred.fix), shape = 15) +
                      geom_sf(data = domain_su, fill = NA, size = 1) +
                      geom_sf(data = inverse_su, fill = "white", color = NA) + 
                      scale_color_gradientn(colors = pal4, name = expression("Population Density / 10 m"^2), limits = c(log(0.0001), log(0.075)), oob = scales::squish, label = rounder) +
                      theme_void() +
                      theme(plot.margin = margin(0,0,0,0, "in"),
                            panel.background = element_blank(),
                            panel.grid = element_blank(),
                            axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            legend.position = "none"))

library(grid)
library(gridExtra)

Figure2A <- arrangeGrob(Fig2A, top = grid::textGrob("A", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0, hjust = 0,
                                                    gp=grid::gpar(fontsize=14, fontface = 2)))
Figure2B <- arrangeGrob(Fig2B, top = grid::textGrob("B", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0, hjust = 0,
                                                    gp=grid::gpar(fontsize=14, fontface = 2)))
Figure2C <- arrangeGrob(Fig2C, top = grid::textGrob("C", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0, hjust = 0,
                                                    gp=grid::gpar(fontsize=14, fontface = 2)))
Figure2D <- arrangeGrob(Fig2D, top = grid::textGrob("D", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0, hjust = 0,
                                                    gp=grid::gpar(fontsize=14, fontface = 2)))
Figure2E <- arrangeGrob(Fig2E, top = grid::textGrob("E", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0, hjust = 0,
                                                    gp=grid::gpar(fontsize=14, fontface = 2)))
Figure2F <- arrangeGrob(Fig2F, top = grid::textGrob("F", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0, hjust = 0,
                                                    gp=grid::gpar(fontsize=14, fontface = 2)))
Figure2G <- arrangeGrob(Fig2G, top = grid::textGrob("G", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0, hjust = 0,
                                                    gp=grid::gpar(fontsize=14, fontface = 2)))
Figure2H <- arrangeGrob(Fig2H, top = grid::textGrob("H", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0, hjust = 0,
                                                    gp=grid::gpar(fontsize=14, fontface = 2)))
Figure2I <- arrangeGrob(Fig2I, top = grid::textGrob("I", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0, hjust = 0,
                                                    gp=grid::gpar(fontsize=14, fontface = 2)))

tiff(file = "~/Monarchs/PostAnalysis/Figure2_new.tiff", res = 600, width = 6.5, height = 5.5, units = "in")
grid.arrange(Figure2A, Figure2B, Figure2C, 
             Figure2D, Figure2E, Figure2F, 
             Figure2G, Figure2H, Figure2I, 
             Legend2,
             layout_matrix = matrix(c(1,4,7,
                                      2,5,8,
                                      3,6,9,
                                      3,6,9,
                                      10,10,10), nrow = 3))
grid.text("Early Spring", x = unit(0.75, "in"), 
          y = unit(5.325, "in"),gp = gpar(fontsize=12))
grid.text("Late Spring", x = unit(2, "in"), 
          y = unit(5.325, "in"),gp = gpar(fontsize=12))
grid.text("Early Summer", x = unit(3.5, "in"), 
          y = unit(5.325, "in"),gp = gpar(fontsize=12))
grid.text("2016", x = unit(0.125, "in"), rot = 90,
          y = unit(4.5, "in"),gp = gpar(fontsize=12))
grid.text("2017", x = unit(0.125, "in"), rot = 90,
          y = unit(2.75, "in"),gp = gpar(fontsize=12))
grid.text("2018", x = unit(0.125, "in"), rot = 90,
          y = unit(1, "in"),gp = gpar(fontsize=12))
grid.text(expression("Population Density / 100 m"^2), x = unit(5.325, "in"), rot = 90,
          y = unit(2.5, "in"),gp = gpar(fontsize=14))
dev.off()


# #Figure 3
# pred.GDD <- matrix(seq(0, 14, length.out = 100) , 100 , 100)
# pred.NDVI <- t(matrix(seq(0, 1, length.out = 100) , 100 , 100))
# 
# mean.GDD <- mean(c(Data$gdd1.avg, NodeDF$gdd1.avg))
# sd.GDD <- sd(c(Data$gdd1.avg, NodeDF$gdd1.avg))
# 
# mean.NDVI <- mean(c(Data$NDVI, NodeDF$NDVI))
# sd.NDVI <- sd(c(Data$NDVI, NodeDF$NDVI))
# 
# pred.GDD <- (pred.GDD - mean.GDD)/sd.GDD
# pred.NDVI <- (pred.NDVI - mean.NDVI)/sd.NDVI
# 
# sp1.Pred <- (exp(beta1[1] * pred.NDVI + beta2[1] * pred.NDVI * pred.NDVI + beta3[1] * pred.GDD + beta4[1] * pred.GDD * pred.GDD) - 1)
# 
# sp2.Pred <- (exp(beta1[2] * pred.NDVI + beta2[2] * pred.NDVI * pred.NDVI + beta3[2] * pred.GDD + beta4[2] * pred.GDD * pred.GDD) - 1)
# 
# su.Pred <- (exp(beta1[3] * pred.NDVI + beta2[3] * pred.NDVI * pred.NDVI + beta3[3] * pred.GDD + beta4[3] * pred.GDD * pred.GDD) - 1)
# 
# pal4 <- colorRampPalette(RColorBrewer::brewer.pal(9, 'Greens'))(100)
# pal5 <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(100)
# 
# ggplot() +
#   geom_raster(data = reshape2::melt(sp1.Pred), aes(x = Var1, y = Var2, fill = value)) +
#   scale_fill_gradientn(colors = pal5, limits = c(-2,2), na.value = RColorBrewer::brewer.pal(11, 'RdBu')[11]) + 
#   geom_point(data = as.data.frame(sp1.points), aes(x = V1, y = V2), shape = 1, size = 2) +
#   scale_y_continuous(expand = c(0,0)) +
#   scale_x_continuous(expand = c(0,0)) +
#   theme_bw() +
#   theme(plot.margin = margin(0.05,0.05,0.25,0.25, unit = "in"),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         axis.text.y = element_blank(),
#         axis.text.x = element_blank(),
#         legend.position = "bottom",
#         legend.key.width = unit(1, "in"),
#         legend.key.height = unit(0.25, "in"),
#         legend.title = element_blank()
#   ) + 
#   guides(fill = guide_colorbar(ticks.colour = "black")) +
#   xlab("\n\n Average GDD") +
#   ylab("NDVI \n\n")
# 
# ggplot() +
#   geom_raster(data = reshape2::melt(sp2.Pred), aes(x = Var1, y = Var2, fill = value)) +
#   scale_fill_gradientn(colors = pal5, limits = c(-100,100), na.value = RColorBrewer::brewer.pal(11, 'RdBu')[11]) +  
#   scale_y_continuous(expand = c(0,0)) +
#   scale_x_continuous(expand = c(0,0)) +
#   theme_bw() +
#   theme(plot.margin = margin(0.05,0.05,0.25,0.25, unit = "in"),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         axis.text.y = element_blank(),
#         axis.text.x = element_blank(),
#         legend.position = "bottom",
#         legend.key.width = unit(1, "in"),
#         legend.key.height = unit(0.25, "in"),
#         legend.title = element_blank()
#   ) + 
#   guides(fill = guide_colorbar(ticks.colour = "black")) +
#   xlab("\n\n Average GDD") +
#   ylab("NDVI \n\n")
# 
# ggplot() +
#   geom_raster(data = reshape2::melt(su.Pred), aes(x = Var1, y = Var2, fill = value)) +
#   scale_fill_gradientn(colors = pal5, limits = c(-100,100), na.value = RColorBrewer::brewer.pal(11, 'RdBu')[11]) + 
#   scale_y_continuous(expand = c(0,0)) +
#   scale_x_continuous(expand = c(0,0)) +
#   theme_bw() +
#   theme(plot.margin = margin(0.05,0.05,0.25,0.25, unit = "in"),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         axis.text.y = element_blank(),
#         axis.text.x = element_blank(),
#         legend.position = "bottom",
#         legend.key.width = unit(1, "in"),
#         legend.key.height = unit(0.25, "in"),
#         legend.title = element_blank()
#   ) + 
#   guides(fill = guide_colorbar(ticks.colour = "black")) +
#   xlab("\n\n Average GDD") +
#   ylab("NDVI \n\n")
# 
# sp1.points <- matrix(NA, nrow = 590, ncol = 2)
# for(i in 1:590){
#   sp1.points[i,] <- as.numeric(tst2[which.min(abs(tst2[,3] - tst3[i]) + abs(tst2[,4] - tst[i])),1:2])
# }
# 

#Figure 3
pred.GDD <- matrix(seq(0, 225, length.out = 100) , 100 , 100)
pred.NDVI <- t(matrix(seq(0, 0.8, length.out = 100) , 100 , 100))

mean.GDD <- mean(c(Data$gdd2, NodeDF$gdd.avg*14))
sd.GDD <- sd(c(Data$gdd2, NodeDF$gdd.avg*14))

mean.NDVI <- mean(c(Data$NDVI, NodeDF$NDVI))
sd.NDVI <- sd(c(Data$NDVI, NodeDF$NDVI))

pred.GDD <- (pred.GDD - mean.GDD)/sd.GDD
pred.NDVI <- (pred.NDVI - mean.NDVI)/sd.NDVI

sp1.Pred <- exp(beta[7] + beta1[1] * pred.NDVI + beta2[1] * pred.NDVI * pred.NDVI + beta3[1] * pred.GDD + beta4[1] * pred.GDD * pred.GDD)

sp2.Pred <- exp(beta[8] + beta1[2] * pred.NDVI + beta2[2] * pred.NDVI * pred.NDVI + beta3[2] * pred.GDD + beta4[2] * pred.GDD * pred.GDD)

su.Pred <- exp(beta[9] + beta1[3] * pred.NDVI + beta2[3] * pred.NDVI * pred.NDVI + beta3[3] * pred.GDD + beta4[3] * pred.GDD * pred.GDD)

pal4 <- colorRampPalette(RColorBrewer::brewer.pal(9, 'Greens'))(100)
pal5 <- colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(100)

gdd.temp <- (Data[Data$period == 1, "gdd2"]$gdd2 - mean.GDD)/sd.GDD
ndvi.temp <- (Data[Data$period == 1, "NDVI"]$NDVI - mean.NDVI)/sd.NDVI
grid.loc <- full_join(reshape2::melt(pred.NDVI), reshape2::melt(pred.GDD), by = c("Var1", "Var2"))
sp1.points <- matrix(NA, nrow = length(gdd.temp), ncol = 2)
for(i in 1:length(gdd.temp)){
  sp1.points[i,] <- as.numeric(grid.loc[which.min(abs(grid.loc[,3] - ndvi.temp[i]) + abs(grid.loc[,4] - gdd.temp[i])),1:2])
}

vline <- which.min(abs(daymet %>% st_drop_geometry(.) %>% 
                         filter(period == 1) %>% 
                         summarise(mu = mean(gdd2.avg)) %>% select(mu) %>% .$mu -
                         seq(0, 225, length.out = 100)))

hline <- which.min(abs(daymet %>% st_drop_geometry(.) %>% 
                         filter(period == 1) %>% 
                         summarise(mu = mean(NDVI)) %>% select(mu) %>% .$mu -
                         seq(0, 0.8, length.out = 100)))



Fig3D <- ggplotGrob(ggplot() +
                      geom_raster(data = reshape2::melt(sp1.Pred), aes(x = Var1, y = Var2, fill = value)) +
                      scale_fill_gradientn(colors = pal4, limits = c(0,0.05), na.value = RColorBrewer::brewer.pal(9, 'Greens')[9]) + 
                      geom_point(data = as.data.frame(sp1.points), aes(x = V1, y = V2), shape = 3, size = 1) +
                      scale_y_continuous(expand = c(0,0), 
                                         labels = round(seq(0, 0.8, length.out = 100)[c(1,25,50,75,100)], digits = 2)) +
                      scale_x_continuous(expand = c(0,0),
                                         labels = round(seq(0, 225, length.out = 100)[c(1,25,50,75,100)], digits = 2)) +
                      theme_bw() +
                      theme(plot.margin = margin(0.05,0.05,0.25,0.25, unit = "in"),
                            panel.background = element_blank(),
                            panel.grid = element_blank(),
                            # legend.position = "bottom",
                            legend.key.width = unit(0.125, "in"),
                            legend.key.height = unit(0.5, "in"),
                            legend.title = element_blank()
                      ) + 
                      guides(fill = guide_colorbar(ticks.colour = "black", frame.colour = "black"))+
                      xlab("GDD") +
                      ylab("NDVI") +
                      geom_vline(aes(xintercept = vline), linetype = "dashed") + 
                      geom_hline(aes(yintercept = hline), linetype = "dashed"))

tiff(filename = "~/Monarchs/PostAnalysis/Figure3D.tiff", width = 6.5, height = 6.5, units = "in", res = 600)
grid::grid.draw(Fig3D)
dev.off()


gdd.temp <- (Data[Data$period == 2, "gdd2"]$gdd2 - mean.GDD)/sd.GDD
ndvi.temp <- (Data[Data$period == 2, "NDVI"]$NDVI - mean.NDVI)/sd.NDVI
grid.loc <- full_join(reshape2::melt(pred.NDVI), reshape2::melt(pred.GDD), by = c("Var1", "Var2"))
sp2.points <- matrix(NA, nrow = length(gdd.temp), ncol = 2)
for(i in 1:length(gdd.temp)){
  sp2.points[i,] <- as.numeric(grid.loc[which.min(abs(grid.loc[,3] - ndvi.temp[i]) + abs(grid.loc[,4] - gdd.temp[i])),1:2])
}

vline <- which.min(abs(daymet %>% st_drop_geometry(.) %>% 
                         filter(period == 2) %>% 
                         summarise(mu = mean(gdd2.avg)) %>% select(mu) %>% .$mu -
                         seq(0, 225, length.out = 100)))

hline <- which.min(abs(daymet %>% st_drop_geometry(.) %>% 
                         filter(period == 2) %>% 
                         summarise(mu = mean(NDVI)) %>% select(mu) %>% .$mu -
                         seq(0, 0.8, length.out = 100)))


Fig3E <- ggplotGrob(ggplot() +
                      geom_raster(data = reshape2::melt(sp2.Pred), aes(x = Var1, y = Var2, fill = value)) +
                      scale_fill_gradientn(colors = pal4) +  
                      geom_point(data = as.data.frame(sp2.points), aes(x = V1, y = V2), shape = 3, size = 1) +
                      scale_y_continuous(expand = c(0,0), 
                                         labels = round(seq(0, 0.8, length.out = 100)[c(1,25,50,75,100)], digits = 2)) +
                      scale_x_continuous(expand = c(0,0),
                                         labels = round(seq(0, 225, length.out = 100)[c(1,25,50,75,100)], digits = 2)) +
                      theme_bw() +
                      theme(plot.margin = margin(0.05,0.05,0.25,0.25, unit = "in"),
                            panel.background = element_blank(),
                            panel.grid = element_blank(),
                            # legend.position = "bottom",
                            legend.key.width = unit(0.125, "in"),
                            legend.key.height = unit(0.5, "in"),
                            legend.title = element_blank()
                      ) + 
                      guides(fill = guide_colorbar(ticks.colour = "black", frame.colour = "black")) +
                      xlab("GDD") +
                      ylab("NDVI") +
                      geom_vline(aes(xintercept = vline), linetype = "dashed") + 
                      geom_hline(aes(yintercept = hline), linetype = "dashed"))

tiff(filename = "~/Monarchs/PostAnalysis/Figure3E.tiff", width = 6.5, height = 6.5, units = "in", res = 600)
grid::grid.draw(Fig3E)
dev.off()

gdd.temp <- (Data[Data$period == 3, "gdd2"]$gdd2 - mean.GDD)/sd.GDD
ndvi.temp <- (Data[Data$period == 3, "NDVI"]$NDVI - mean.NDVI)/sd.NDVI
grid.loc <- full_join(reshape2::melt(pred.NDVI), reshape2::melt(pred.GDD), by = c("Var1", "Var2"))
su.points <- matrix(NA, nrow = length(gdd.temp), ncol = 2)
for(i in 1:length(gdd.temp)){
  su.points[i,] <- as.numeric(grid.loc[which.min(abs(grid.loc[,3] - ndvi.temp[i]) + abs(grid.loc[,4] - gdd.temp[i])),1:2])
}

vline <- which.min(abs(daymet %>% st_drop_geometry(.) %>% 
                         filter(period == 3) %>% 
                         summarise(mu = mean(gdd2.avg)) %>% select(mu) %>% .$mu -
                         seq(0, 225, length.out = 100)))

hline <- which.min(abs(daymet %>% st_drop_geometry(.) %>% 
                         filter(period == 3) %>% 
                         summarise(mu = mean(NDVI)) %>% select(mu) %>% .$mu -
                         seq(0, 0.8, length.out = 100)))

Fig3F <- ggplotGrob(ggplot() +
                      geom_raster(data = reshape2::melt(su.Pred), aes(x = Var1, y = Var2, fill = value)) +
                      scale_fill_gradientn(colors = pal4) + 
                      geom_point(data = as.data.frame(su.points), aes(x = V1, y = V2), shape = 3, size = 1) +
                      # scale_y_continuous(labels=LETTERS, breaks=1:26, limits=c(1,26))
                      scale_y_continuous(expand = c(0,0), 
                                         labels = round(seq(0, 0.8, length.out = 100)[c(1,25,50,75,100)], digits = 2)) +
                      scale_x_continuous(expand = c(0,0),
                                         labels = round(seq(0, 255, length.out = 100)[c(1,25,50,75,100)], digits = 2)) +
                      theme_bw() +
                      theme(plot.margin = margin(0.05,0.05,0.25,0.25, unit = "in"),
                            panel.background = element_blank(),
                            panel.grid = element_blank(),
                            # legend.position = "bottom",
                            legend.key.width = unit(0.125, "in"),
                            legend.key.height = unit(0.5, "in"),
                            legend.title = element_blank()
                      ) + 
                      guides(fill = guide_colorbar(ticks.colour = "black", frame.colour = "black")) +
                      xlab("GDD") +
                      ylab("NDVI")  +
                      geom_vline(aes(xintercept = vline), linetype = "dashed") + 
                      geom_hline(aes(yintercept = hline), linetype = "dashed"))

tiff(filename = "~/Monarchs/PostAnalysis/Figure3F.tiff", width = 6.5, height = 6.5, units = "in", res = 600)
grid::grid.draw(Fig3F)
dev.off()

Effect <- data.frame(param = rep(c("beta1", "beta2", "beta3", "beta4"), each = 3),
                     period = factor(rep(1:3, 4)),
                     mean.den = c(beta1, beta2, beta3, beta4))

Effect$lower <- NA
Effect$upper <- NA

Effect[1,4:5] <- confint(tmbprofile(obj1, lincomb = c(0,1,rep(0,34))))
Effect[2,4:5] <- confint(tmbprofile(obj1, lincomb = c(0,0,1,rep(0,33))))
Effect[3,4:5] <- confint(tmbprofile(obj1, lincomb = c(0,0,0,1,rep(0,32))))
Effect[4,4:5] <- confint(tmbprofile(obj1, lincomb = c(0,0,0,0,1,rep(0,31))))
Effect[5,4:5] <- confint(tmbprofile(obj1, lincomb = c(0,0,0,0,0,1,rep(0,30))))
Effect[6,4:5] <- confint(tmbprofile(obj1, lincomb = c(rep(0,6),1,rep(0,29))))
Effect[7,4:5] <- confint(tmbprofile(obj1, lincomb = c(rep(0,7),1,rep(0,28))))
Effect[8,4:5] <- confint(tmbprofile(obj1, lincomb = c(rep(0,8),1,rep(0,27))))
Effect[9,4:5] <- confint(tmbprofile(obj1, lincomb = c(rep(0,9),1,rep(0,26))))
Effect[10,4:5] <- confint(tmbprofile(obj1, lincomb = c(rep(0,10),1,rep(0,25))))
Effect[11,4:5] <- confint(tmbprofile(obj1, lincomb = c(rep(0,11),1,rep(0,24))))
Effect[12,4:5] <- confint(tmbprofile(obj1, lincomb = c(rep(0,12),1,rep(0,23))))

Fig3A <- ggplotGrob(ggplot(Effect, aes(y = mean.den, x = param, col = period)) +
                      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
                      geom_point(position = position_dodge(0.5)) +
                      geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.5), width = 0.25) +
                      scale_color_manual(name = "Stage", values = RColorBrewer::brewer.pal(4, 'Greens')[2:4]) +
                      scale_x_discrete(labels = c("beta1" = "NDVI", "beta2" = expression(NDVI^2),
                                                  "beta3" = "GDD", "beta4" = expression(GDD^2))) +
                      theme_few() +
                      theme(legend.position = "top") +
                      xlab("") +
                      ylab("Covariate effect (log-scale)"))

tiff(filename = "~/Monarchs/PostAnalysis/Figure3A.tiff", width = 6.5, height = 6.5, units = "in", res = 600)
grid::grid.draw(Fig3A)
dev.off()


#Figure XX

pred.GDD <- seq(0, 225, length.out = 100)
pred.NDVI <- seq(0, 0.8, length.out = 100)

pred.GDD <- (pred.GDD - mean.GDD)/sd.GDD
pred.NDVI <- (pred.NDVI - mean.NDVI)/sd.NDVI

sp1.NDVI.mean <- exp(beta[7] + beta1[1] * pred.NDVI + beta2[1] * pred.NDVI * pred.NDVI)
sp2.NDVI.mean<- exp(beta[8] + beta1[2] * pred.NDVI + beta2[2] * pred.NDVI * pred.NDVI)
su.NDVI.mean <- exp(beta[9] + beta1[3] * pred.NDVI + beta2[3] * pred.NDVI * pred.NDVI)

sp1.NDVI.lower <- exp(beta[7] + Effect[1,4] * pred.NDVI + Effect[4,4] * pred.NDVI * pred.NDVI)
sp2.NDVI.lower <- exp(beta[8] + Effect[2,4] * pred.NDVI + Effect[5,4] * pred.NDVI * pred.NDVI)
su.NDVI.lower <- exp(beta[9] + Effect[3,4] * pred.NDVI + Effect[6,4] * pred.NDVI * pred.NDVI)

sp1.NDVI.upper <- exp(beta[7] + Effect[1,5] * pred.NDVI + Effect[4,5] * pred.NDVI * pred.NDVI)
sp2.NDVI.upper <- exp(beta[8] + Effect[2,5] * pred.NDVI + Effect[5,5] * pred.NDVI * pred.NDVI)
su.NDVI.upper <- exp(beta[9] + Effect[3,5] * pred.NDVI + Effect[6,5] * pred.NDVI * pred.NDVI)


#Avg GDD

sp1.GDD.mean <- exp(beta[7] + beta3[1] * pred.GDD + beta4[1] * pred.GDD * pred.GDD)
sp2.GDD.mean<- exp(beta[8] + beta3[2] * pred.GDD + beta4[2] * pred.GDD * pred.GDD)
su.GDD.mean <- exp(beta[9] + beta3[3] * pred.GDD + beta4[3] * pred.GDD * pred.GDD)

sp1.GDD.lower <- exp(beta[7] + Effect[7,4] * pred.GDD + Effect[10,4] * pred.GDD * pred.GDD)
sp2.GDD.lower <- exp(beta[8] + Effect[8,4] * pred.GDD + Effect[11,4] * pred.GDD * pred.GDD)
su.GDD.lower <- exp(beta[9] + Effect[9,4] * pred.GDD + Effect[12,4] * pred.GDD * pred.GDD)

sp1.GDD.upper <- exp(beta[7] + Effect[7,5] * pred.GDD + Effect[10,5] * pred.GDD * pred.GDD)
sp2.GDD.upper <- exp(beta[8] + Effect[8,5] * pred.GDD + Effect[11,5] * pred.GDD * pred.GDD)
su.GDD.upper <- exp(beta[9] + Effect[9,5] * pred.GDD + Effect[12,5] * pred.GDD * pred.GDD)

pred.GDD <- seq(0, 225, length.out = 100)
pred.NDVI <- seq(0, 0.8, length.out = 100)

Pred.df <- data.frame(NDVI.mean = c(sp1.NDVI.mean, sp2.NDVI.mean, su.NDVI.mean),
                      NDVI.lower = c(sp1.NDVI.lower, sp2.NDVI.lower, su.NDVI.lower),
                      NDVI.upper = c(sp1.NDVI.upper, sp2.NDVI.upper, su.NDVI.upper),
                      GDD.mean = c(sp1.GDD.mean, sp2.GDD.mean, su.GDD.mean),
                      GDD.lower = c(sp1.GDD.lower, sp2.GDD.lower, su.GDD.lower),
                      GDD.upper = c(sp1.GDD.upper, sp2.GDD.upper, su.GDD.upper),
                      period = rep(1:3, each = 100),
                      NDVI = pred.NDVI,
                      GDD = pred.GDD)

Fig3B <- ggplotGrob(ggplot() +
                      geom_ribbon(data = Pred.df %>% filter(period == 1), aes(x = NDVI, ymin = NDVI.lower, ymax = NDVI.upper), alpha = 0.75, fill = RColorBrewer::brewer.pal(4, 'Greens')[2]) +
                      geom_ribbon(data = Pred.df %>% filter(period == 2), aes(x = NDVI, ymin = NDVI.lower, ymax = NDVI.upper), alpha = 0.75, fill = RColorBrewer::brewer.pal(4, 'Greens')[3]) +
                      geom_ribbon(data = Pred.df %>% filter(period == 3), aes(x = NDVI, ymin = NDVI.lower, ymax = NDVI.upper), alpha = 0.75, fill = RColorBrewer::brewer.pal(4, 'Greens')[4]) +
                      geom_line(data = Pred.df %>% filter(period == 1), aes(x = NDVI, y = NDVI.mean)) +
                      geom_line(data = Pred.df %>% filter(period == 2), aes(x = NDVI, y = NDVI.mean), linetype = "dashed") +
                      geom_line(data = Pred.df %>% filter(period == 3), aes(x = NDVI, y = NDVI.mean), linetype = "dotted") +
                      theme_few() +
                      ylab(expression("Population Abundance / 100 m"^2)))

tiff(filename = "~/Monarchs/PostAnalysis/Figure3B.tiff", width = 6.5, height = 6.5, units = "in", res = 600)
grid::grid.draw(Fig3B)
dev.off()

Fig3C <- ggplotGrob(ggplot() +
                      geom_ribbon(data = Pred.df %>% filter(period == 1), aes(x = GDD, ymin = GDD.lower, ymax = GDD.upper), alpha = 0.75, fill = RColorBrewer::brewer.pal(4, 'Greens')[2]) +
                      geom_ribbon(data = Pred.df %>% filter(period == 2), aes(x = GDD, ymin = GDD.lower, ymax = GDD.upper), alpha = 0.75, fill = RColorBrewer::brewer.pal(4, 'Greens')[3]) +
                      geom_ribbon(data = Pred.df %>% filter(period == 3), aes(x = GDD, ymin = GDD.lower, ymax = GDD.upper), alpha = 0.75, fill = RColorBrewer::brewer.pal(4, 'Greens')[4]) +
                      geom_line(data = Pred.df %>% filter(period == 1), aes(x = GDD, y = GDD.mean)) +
                      geom_line(data = Pred.df %>% filter(period == 2), aes(x = GDD, y = GDD.mean), linetype = "dashed") +
                      geom_line(data = Pred.df %>% filter(period == 3), aes(x = GDD, y = GDD.mean), linetype = "dotted") +
                      theme_few() +
                      xlab("GDD") +
                      ylab(expression("Population Abundance / 100 m"^2)))

tiff(filename = "~/Monarchs/PostAnalysis/Figure3C.tiff", width = 6.5, height = 6.5, units = "in", res = 600)
grid::grid.draw(Fig3C)
dev.off()


Figure3A <- arrangeGrob(Fig3A, top = grid::textGrob("A", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0, hjust = 0,
                                                    gp=grid::gpar(fontsize=18, fontface = 2)))
Figure3B <- arrangeGrob(Fig3B, top = grid::textGrob("B", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0, hjust = 0,
                                                    gp=grid::gpar(fontsize=18, fontface = 2)))
Figure3C <- arrangeGrob(Fig3C, top = grid::textGrob("C", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0, hjust = 0,
                                                    gp=grid::gpar(fontsize=18, fontface = 2)))
Figure3D <- arrangeGrob(Fig3D, top = grid::textGrob("D", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0, hjust = 0,
                                                    gp=grid::gpar(fontsize=18, fontface = 2)))
Figure3E <- arrangeGrob(Fig3E, top = grid::textGrob("E", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0, hjust = 0,
                                                    gp=grid::gpar(fontsize=18, fontface = 2)))
Figure3F <- arrangeGrob(Fig3F, top = grid::textGrob("F", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0, hjust = 0,
                                                    gp=grid::gpar(fontsize=18, fontface = 2)))



tiff(filename = "~/Monarchs/PostAnalysis/Figure3.tiff", width = 8, height = 11.5, units = "in", res = 600)
grid.arrange(Figure3A, Figure3B, Figure3C, Figure3D, Figure3E, Figure3F, layout_matrix = matrix(c(1,2,3,4,5,6), nrow = 3))
dev.off()

# ggplot() + 
#   geom_point(data = Output, aes(y = exp(mean.den), x = period, color = year)) +
#   geom_errorbar(data = Output, aes(ymin = exp(lower.den), ymax = exp(upper.den), x = period, color = year))