##################################################
## Referential radial file construction
## Pipeline example
##
## vesrion date: 12.2020
## author: Estelle Noyer, noyer.estelle@gmail.com
##################################################

#####
## Ddl
#####

rm(list=ls())

library(grDevices)

# paths
#pathScript <-"RRF"
#pathR <- "RRF/dataset"
#pathS <- "RRF/output"
#pathG <- "RRF/graph"

# function
setwd(pathScript)
source("(Fct)RRF_estimation.R")

# dataframe
setwd(pathR)
coord  <- read.table("ash_cell_coordinates_clean.txt", header = TRUE, dec = ".", sep = " ")

######################################################## APPLICATION OF THE FUNCTIONS #####

######
# alpha shape #####
######
setwd(pathG)
pdf("ishape_plot_ash.pdf", width = 6, height = 15)

dataV  <- coord[which(coord$type == "Vessel"),]
dataNV <- coord[which(coord$type != "Vessel"),]

outAS <- ishape(dataset = dataNV, dataV = dataV, make.plot = TRUE)
dev.off()

rm(dataNV, dataV)

## export
setwd(pathS)
write.table(outAS, "ishape_output_ash.txt", row.names = FALSE)


#####
# RRF: 5 steps ####
#####

# 1. cambial line and layer characteristics
#widthL = 40 for ring-porous species
#widthL = 30 for diffuse-porous species
outI     <- ccls(dataset = coord, dataIS = outAS, widthL = 40)
carCfile <- as.data.frame(outI[1])
carCell  <- as.data.frame(outI[2])
rm(outI)

## export
setwd(pathS)
write.table(carCfile, "output_ccls_characteristics_ash.txt", row.names = TRUE)
write.table(carCell, "output_ccls_summary_ash.txt", row.names = TRUE)

# 2. polygon creation

setwd(pathG)
pdf("sectors_plot_ash.pdf", width = 6, height = 20)

SecEdges <- sectors(dataset = coord, cFile = carCfile, cell = carCell, dataIS = outAS, make.plot = TRUE)

dev.off()

## export
setwd(pathS)
write.table(SecEdges, "output_sectors_ash.txt", row.names = TRUE)


# 3. attribution idSector to each cell and abs/presence type cell per sector
setwd(pathG)
pdf("labelSec_plot_ash.pdf", width = 6, height = 20)

outM <- labelSec(dataset = coord, EdgeS = SecEdges, make.plot = TRUE)
outA <- as.data.frame(outM[1])
outP <- as.data.frame(outM[2])

dev.off()

## export
setwd(pathS)
write.table(outA, "output_labelSec_cell_sectorID_ash.txt", row.names = TRUE)
write.table(outP, "output_labelSec_sector_summary_ash.txt", row.names = TRUE)


# 4. vessel sectors ####
## 4.1 vessel neighbor cell id ####
outN   <- vNeigh(dataset = outA)
Vcar   <- as.data.frame(outN[1])
Vneigh <- as.data.frame(outN[2])

### export
setwd(pathS)
write.table(Vcar, "output_vNeigh_vessel_neighbors_ash.txt", row.names = TRUE)
write.table(Vneigh, "output_vNeigh_cell_characteristics_ash.txt", row.names = TRUE)


## 4.2 pool of sectors attributed to vessel  ####
setwd(pathG)
pdf("vesselS_plot_ash.pdf", width = 6, height = 20)
if(all(nrow(Vcar) != 0, nrow(Vneigh) != 0)){
  outV <- vesselS(dataset = outA, EdgeS = outP, CarV = Vcar, NeighV = Vneigh, make.plot = TRUE)
}
dev.off()

### export
setwd(pathS)
write.table(outV, "output_vesselS_vessel_sectors_ash.txt", row.names = TRUE)


# 5. estimation of the number of APF for one cambial file producing exclusively APF and their number per sector ####
setwd(pathG)
pdf("rrf_plot_ash.pdf", width = 10, height = 20)

outC  <- rrf(dataset = outA, EdgeS = outP, Vsec = outV, make.plot = TRUE)
cLcal <- as.data.frame(outC[1])
cLsum <- as.data.frame(outC[2])

dev.off()

## export
setwd(pathS)
write.table(cLcal, "output_rrf_produced_cell_number_ash.txt", row.names = TRUE)
write.table(cLsum, "output_rrf_produced_cell_image_ash.txt", row.names = TRUE)


# 6. Vessel indexation ########
if(nrow(outV) != 0){
  outI <- vesselID(dataV = outV, dataC = Vcar, pLab = outA, rrf = cLcal)
  outIndex <- as.data.frame(outI[1])
  newRRF   <- as.data.frame(outI[2])
}

# 7. yearly relative vessel index ########
outFinal <- merge(outIndex, cLsum[, c("idD", "id", "DOY", "RRF_avg")], by = "idD")
outFinal <- merge(outFinal, newRRF[, c("idD", "cumulCellLay")], by = "idD")
outFinal$vessel.Relindex <- (outFinal$vessel.index * 100) / outFinal$cumulCellLay

outFinal[which(outFinal$vessel.Relindex > 100), ]


# 8. export ########
setwd(pathS)
write.table(outIndex, "output_vesselID_vessel_features_ash.txt", row.names = TRUE)
write.table(outFinal, "output_pipeline_ash.txt", row.names = TRUE)


########################################################################################################## END #####