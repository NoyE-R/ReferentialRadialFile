##################################################
## Referential radial file construction
## Used Functions
##
## vesrion date: 12.2020
## author: Estelle Noyer, noyer.estelle@gmail.com
##################################################

#################### FUNCTION ccls: characteristics of cambial layers/lines #####
# function to estimate the number of layer for each cambial line
# return list of 2 objects in outC: cambial line characteristics dataframe and first + cambial cells features 

#' @dataset: dataframe with reoriented data (X.cor, Y.cor) and "type" column WITHOUT vessel data
#' @dataV: dataframe with reoriented data (X.cor, Y.cor) and "type" column ONLY WITH vessel data

ishape <- function(dataset, dataV, make.plot = TRUE){
  require(alphahull) || install.packages(alphahull)
  require(sp) || install.packages(sp)
  library(alphahull)
  library(reshape2)
  
  # to avoid effect of big vessel on the alpha shape: different alpha angle following maximal vessel area
  if(all(is.na(dataV$Area.v))){
    alpha <- 135
  }else if(max(dataV$Area.v, na.rm = TRUE) > 900000){
    alpha <- 655
  }else if(max(dataV$Area.v, na.rm = TRUE) > 200000){
    alpha <- 500
  }else if(max(dataV$Area.v, na.rm = TRUE) > 50000){
    alpha <- 300
  }else if(max(dataV$Area.v, na.rm = TRUE) > 9000){
    alpha <- 250
  }else if(max(dataV$Area.v, na.rm = TRUE) == 0){
    alpha <- 655
  }else{
    alpha <- 135
  }
  
  if(nrow(dataset[which(duplicated(paste(dataset$X.cor, dataset$Y.cor))),]) != 0){
    # conversion to spatial points dataframe and removing the duplicated data
    spaData <- sp::SpatialPointsDataFrame(coords = cbind(dataset$X.cor, dataset$Y.cor), data = dataset)
    spaData <- spaData[-sp::zerodist(spaData)[, 1], ]
    AS <- alphahull::ashape(coordinates(spaData)[, 1], coordinates(spaData)[, 2], alpha = alpha)
  }else{
    AS <- alphahull::ashape(dataset[, c("X.cor", "Y.cor")], alpha = alpha)
  }
  
  # plot
  plotI <- plot(AS, main = unique(dataset$idD))
  
  # vector of cell at the edges
  outIshape      <- data.frame(AS$edges[, 1:6])
  outIshape$idD  <- unique(dataset$idD)
  outIshape$id   <- unique(dataset$id)
  outIshape$tree <- unique(dataset$tree)
  outIshape$year <- unique(dataset$year)
  outIshape$DOY  <- unique(dataset$DOY)
  outIshape$date <- unique(dataset$date)
  
  # return
  return(outIshape)
}

#### end function ccls ##########

#################### FUNCTION ccls: characteristics of cambial layers/lines #####
# function to estimate the number of layer for each cambial line
# return list of 2 objects in outC: cambial line characteristics dataframe and first + cambial cells features 

#' @dataset: dataframe with reoriented data (X.cor, Y.cor) and "type" column
#' @dataIS: dataframe output of ishape()
#' @widthL: width of the layer. Default = 30 um

ccls <- function(dataset, dataIS, widthL = 30){
  j      <- NULL
  carCC  <- NULL # dataframe with cambial layer and cambial line characteristics
  cFile  <- NULL # number of cambial cell = cambial line index in the image
  nCc    <- NULL # length of the cambial layer
  cD     <- NULL # mean cambial cell diameter
  newX   <- NULL # corrected value of the x coordinate of the edges of the cambial cell
  fC     <- NULL # subset of the dataframe ThisId of the cell with Y.cor < 60
  fCmin  <- NULL # min X.cor of fC
  fCmax  <- NULL # max X.cor of fC
  dif    <- NULL # difference between the X max and min coordinates of fC
  newDif <- NULL # difference between the X max and min coordinates of new edge coordinates
  ThisC  <- NULL
  Ymin   <- NULL # min Y.cor value in fC
  Ycell  <- NULL # type of cell of min Y.cor value in fC
  cell   <- NULL # dataframe with cambial cell and first cells characteristics per image
  outC   <- NULL # output of the function
  
  # 1. calculation of the length of the cambial layer and edges of each cambial lines #####
  ## cambial cells position and their "theoretical transversal diameter"
  carCC <- data.frame(id = unique(dataset$idD), idC = unique(dataset[dataset$Cells == "Cambial Row", ]$idC),
                      cFile = NA, X.cor = NA, Y.cor = NA,
                      Xmin = NA, Xmax = NA, Ymin = NA, Ymax = NA,
                      TransvDiaTheo = NA, perLength = NA, cumperLength = NA,
                      Xinter.min = NA, Xinter.max = NA, Ymin.line = NA, layer = NA)
  
  ## X and Y coordinates
  for(j in 1:nrow(carCC)){
    carCC[j, ]$X.cor <- round(dataset[dataset$idC == carCC[j, ]$idC,]$X.cor, 3)
    carCC[j, ]$Y.cor <- round(dataset[dataset$idC == carCC[j, ]$idC,]$Y.cor, 3)
  }
  rm(j)
  
  ## order and indexation of the cambial line
  carCC       <- carCC[order(carCC$X.cor), ]  
  carCC$cFile <- c(1:nrow(carCC))
  
  ## calculation of the Xmin and Xmax
  for(j in 1:nrow(carCC)){
    ### Xmin and Xmax values
    if(j == 1){
      carCC[j, ]$Xmin <- round(carCC[j, ]$X.cor - (carCC[j + 1, ]$X.cor - carCC[j, ]$X.cor)/2, 3)
      carCC[j, ]$Xmax <- round((carCC[j, ]$X.cor + carCC[j + 1, ]$X.cor)/2, 3)
      
      carCC[j, ]$Ymin <- round(carCC[j, ]$Y.cor - (carCC[j + 1, ]$Y.cor - carCC[j, ]$Y.cor)/2, 3)
      carCC[j, ]$Ymax <- round((carCC[j, ]$Y.cor + carCC[j + 1, ]$Y.cor)/2, 3)
    }else{
      if(j == max(nrow(carCC))){
        carCC[j, ]$Xmin <- round((carCC[j, ]$X.cor + carCC[j - 1, ]$X.cor)/2, 3)
        carCC[j, ]$Xmax <- round(carCC[j, ]$X.cor + (carCC[j, ]$X.cor - carCC[j - 1, ]$X.cor)/2, 3)
        
        carCC[j, ]$Ymin <- round((carCC[j, ]$Y.cor + carCC[j - 1, ]$Y.cor)/2, 3)
        carCC[j, ]$Ymax <- round(carCC[j, ]$Y.cor + (carCC[j, ]$Y.cor - carCC[j - 1, ]$Y.cor)/2, 3)
      }else{
        carCC[j, ]$Xmin <- round((carCC[j, ]$X.cor + carCC[j - 1, ]$X.cor)/2, 3)
        carCC[j, ]$Xmax <- round((carCC[j, ]$X.cor + carCC[j + 1, ]$X.cor)/2, 3)
        
        carCC[j, ]$Ymin <- round((carCC[j, ]$Y.cor + carCC[j - 1, ]$Y.cor)/2, 3)
        carCC[j, ]$Ymax <- round((carCC[j, ]$Y.cor + carCC[j + 1, ]$Y.cor)/2, 3)
      }
    }
  }
  rm(j)
  
  ## homogenization of the edges values : need to be in other loop because need to have all the values to make the homogenization
  for(j in 1:nrow(carCC)){
    if(j != max(nrow(carCC))){
      if(carCC[j, ]$Xmax != carCC[j + 1, ]$Xmin){
        newX                <- NULL
        newX                <- round((carCC[j, ]$Xmax + carCC[j + 1, ]$Xmin) / 2, 3)
        carCC[j, ]$Xmax     <- newX
        carCC[j + 1, ]$Xmin <- newX
        rm(newX)
      }
    }
    
    ### theoretical transverse diameter and perLength values
    carCC[j, ]$TransvDiaTheo <- sqrt((carCC[j, ]$Xmin - carCC[j, ]$Xmax)^2 + (carCC[j, ]$Ymin - carCC[j, ]$Ymax)^2)
  }
  rm(j)
  
  ## number of cambial line
  cFile <- c(1:length(dataset[dataset$type == "Cambial",]$type))
  ## length cambial layer +/- equal to the width of the zone
  nCc  <- sum(carCC$TransvDiaTheo)
  ## mean transverse cell diameter
  cD  <- nCc / max(cFile)
  
  ## percentage of length for each radial line
  for(j in 1:nrow(carCC)){ carCC[j,]$perLength  <- round(carCC[j,]$TransvDiaTheo * 100 / nCc, 2) }
  rm(j)
  
  ## cumul perLength  
  carCC$cumperLength  <- cumsum(carCC$perLength)
  
  # 2. First produced cells coordinates #####
  if(unique(dataset$TWood == "RP")){
    if(nrow(dataset[dataset$type == "Vessel", ]) > 0){
      if(min(dataset[dataset$type == "Vessel", ]$Y.cor) < 100){
        fC <- dataset[dataset$Y.cor <= min(dataset[dataset$type == "Vessel", ]$Y.cor) + widthL,]
      }else{
        fC <- dataset[dataset$Y.cor <= 100,]
      }
    }else{
      fC <- dataset[dataset$Y.cor <= 100,]
    }
    
    fCAS <- dataIS[dataIS$y1 <= 300 & dataIS$y2 <= 300,]
    
  }else{
    if(nrow(dataset[dataset$type == "Vessel", ]) > 0){
      if(min(dataset[dataset$type == "Vessel", ]$Y.cor) < (widthL * 3)){
        fC <- dataset[dataset$Y.cor <= (widthL * 3),]
      }else{
        fC <- dataset[dataset$Y.cor <= (widthL * 1.5),]
      }
    }else{
      fC <- dataset[dataset$Y.cor <= (widthL * 1.5),]
    }
    
    fCAS <- dataIS[dataIS$y1 <= (widthL * 1.5) & dataIS$y2 <= (widthL * 1.5),]
  }
  
  ## X values of the alpha shape
  XcorAS <- unique(c(fCAS$x1, fCAS$x2))
  YcorAS <- unique(c(fCAS$y1, fCAS$y2))
  
  ## max X.cor of fC
  if(fC[which.max(fC$X.cor), ]$type == "Vessel"){
    ### if it is a vessel in the right bottom corner
    fCmax <- fC[which.max(fC$X.cor), ]$X.cor + ((fC[which.max(fC$X.cor),]$Feret.v + fC[which.max(fC$X.cor),]$MinFeret.v) / 2)
    if(fCmax > max(dataset$X.cor, na.rm = TRUE)){ fCmax <- max(dataset$X.cor, na.rm = TRUE) + cD/4 }else{}
  }else{
    fCmax <- fC[which.max(fC$X.cor), ]$X.cor
  }
  
  ## min X.cor of fC
  if(fC[which.min(fC$X.cor), ]$type == "Vessel"){
    ### if it is a vessel in the left bottom corner
    fCmin <- (fC[which.min(fC$X.cor), ]$X.cor - ((fC[which.min(fC$X.cor),]$Feret.v + fC[which.min(fC$X.cor),]$MinFeret.v) / 2))
    if(fCmin < 0){ fCmin <- 0 }
  }else{
    fCmin <- fC[which.min(fC$X.cor), ]$X.cor
  }
  
  ## extraction of the ashape
  if(unique(dataset$type == "DP")){
    if(all(fCmax %in% XcorAS, abs(fCmax - max(XcorAS)) > cD)){
      fCmaxAS <- fCmax
    }else{
      fCmaxAS <- max(XcorAS)
    }
    
    if(all(fCmin %in% XcorAS, abs(fCmin - min(XcorAS)) > cD)){
      fCminAS <- fCmin
    }else{
      fCminAS <- min(XcorAS)
    }
  }else{
    fCmaxAS <- max(XcorAS)
    fCminAS <- min(XcorAS)
  }
  
  
  ## X edges difference (raw and ashape values)
  dif <- fCmax - fCmin
  if(unique(dataset$type == "DP")){
    difAS <- fCmaxAS - fCminAS
  }else{
    difAS <- nCc / 1.8
  }
  
  ## comparison of dif (length of the first produced cells) to difAS (length of first layer AS), in the case of ring porous species and vessel dif < difAS : correction
  if(dif < difAS){ ## ## dif < difAS: RP + Vessel case ##
    ### X values
    if(fC[which.max(fC$X.cor), ]$Y.cor - fC[which.min(fC$X.cor), ]$Y.cor > 0){ ## ## If positive, going to the right so left part missing
      carCC[nrow(carCC), ]$Xinter.max <- fCmaxAS
      newDif                          <- carCC[nrow(carCC), ]$Xinter.max
      carCC[nrow(carCC), ]$Xinter.min <- carCC[nrow(carCC), ]$Xinter.max - (carCC[nrow(carCC), ]$perLength * newDif / 100)
      
      for(j in (nrow(carCC) - 1):1){
        carCC[j, ]$Xinter.max <- carCC[j + 1, ]$Xinter.min
        carCC[j, ]$Xinter.min <- carCC[j, ]$Xinter.max - (carCC[j, ]$perLength * newDif / 100)
      }
      rm(j)
      
    }else{ ## ## If negative, going to the left so right part missing
      carCC[1, ]$Xinter.min           <- fCminAS
      carCC[nrow(carCC), ]$Xinter.max <- fCmaxAS
      newDif                          <- carCC[nrow(carCC), ]$Xinter.max - carCC[1, ]$Xinter.min
      carCC[1, ]$Xinter.max           <- carCC[1, ]$Xinter.min + (carCC[1, ]$cumperLength * newDif / 100)
      
      for(j in 2:(nrow(carCC)-1)){
        carCC[j, ]$Xinter.min <- carCC[1, ]$Xinter.min + (carCC[j - 1, ]$cumperLength * newDif / 100)
        carCC[j, ]$Xinter.max <- carCC[1, ]$Xinter.min + (carCC[j, ]$cumperLength * newDif / 100)
      }
      rm(j)
      
      carCC[nrow(carCC), ]$Xinter.min <- carCC[nrow(carCC) - 1, ]$Xinter.max
      
    }
    
    ### Ymin.line: based on the ashape
    for(j in (1:nrow(carCC))){
      #### extraction of the ashape value included cFile j
      ThisC <- rbind(fCAS[which(fCAS$x1 < carCC[j, ]$Xinter.min & fCAS$x2 > carCC[j, ]$Xinter.min), ],
                     fCAS[which(fCAS$x1 > carCC[j, ]$Xinter.min & fCAS$x2 < carCC[j, ]$Xinter.min), ],
                     fCAS[which(fCAS$x1 < carCC[j, ]$Xinter.max & fCAS$x2 > carCC[j, ]$Xinter.max), ],
                     fCAS[which(fCAS$x1 > carCC[j, ]$Xinter.max & fCAS$x2 < carCC[j, ]$Xinter.max), ])
      if(nrow(ThisC) == 0){
        carCC[j, ]$Ymin.line <- 0
      }else{
        carCC[j, ]$Ymin.line <- min(ThisC$y1, ThisC$y2)
      }
    }
    rm(j)
    
  }else{ ## ## dif > difAS, so theoretically no edged vessel as cell 1 detected in fC ##
    ### X and Y values of the first cells row
    for(j in 1:nrow(carCC)){
      #### X values
      if(j == 1){
        carCC[j, ]$Xinter.min <- fCmin
        carCC[j, ]$Xinter.max <- fCmin + (carCC[j, ]$cumperLength * dif / 100)
      }else{
        if(j == max(nrow(carCC))){
          carCC[j, ]$Xinter.min <- fCmin + (carCC[j-1, ]$cumperLength * dif / 100)
          carCC[j, ]$Xinter.max <- fCmax
        }else{
          carCC[j, ]$Xinter.min <- fCmin + (carCC[j-1, ]$cumperLength * dif / 100)
          carCC[j, ]$Xinter.max <- fCmin + (carCC[j, ]$cumperLength * dif / 100)
        }
      }
      
      #### Ymin.line
      if(length(fC[fC$X.cor >= (carCC[j, ]$Xinter.min - 1) & fC$X.cor < (carCC[j, ]$Xinter.max + 1), ]$Y.cor) == 0){
        if(j == 1){
          carCC[j, ]$Ymin.line <- 0
        }else{
          carCC[j, ]$Ymin.line <- carCC[j - 1, ]$Ymin.line - 1
        }
        
      }else{
        ##### if the line starts by a vessel, the value is biased. Calculation of the correction
        Ymin  <- min(fC[fC$X.cor >= (carCC[j, ]$Xinter.min - 1) & fC$X.cor < (carCC[j, ]$Xinter.max + 1),]$Y.cor, na.rm = TRUE)
        Ycell <- fC[fC$Y.cor == Ymin, ]$type
        
        if(any(Ycell %in% "Vessel")){
          if(fC[fC$Y.cor == Ymin  & fC$type == "Vessel", ]$Feret.v == 0){
            if(j ==1){
              carCC[j, ]$Ymin.line <- 0
            }else{
              carCC[j, ]$Ymin.line <- carCC[j - 1, ]$Ymin.line - 1
            }
          }else if((Ymin - ((fC[fC$Y.cor == Ymin & fC$type == "Vessel", ]$Feret.v + fC[fC$Y.cor == Ymin & fC$type == "Vessel", ]$MinFeret.v) / 2)) < 0){
            if(j ==1){
              carCC[j, ]$Ymin.line <- 0
            }else{
              carCC[j, ]$Ymin.line <- carCC[j - 1, ]$Ymin.line - 1
            }
          }else{
            carCC[j, ]$Ymin.line <- Ymin - ((fC[fC$Y.cor == Ymin & fC$type == "Vessel", ]$Feret.v + fC[fC$Y.cor == Ymin & fC$type == "Vessel", ]$MinFeret.v) / 2)
          }
          
        }else{ ## test on the off set between Ymin.line of the j cFile and the previous cFile
          carCC[j, ]$Ymin.line <- Ymin
          
          if(all(j == 1, Ymin >= widthL)){
              carCC[1, ]$Ymin.line <- mean(c(Ymin, widthL))
          }else if(all(j != 1, Ymin > carCC[j-1, ]$Ymin.line + widthL/2)){
              carCC[j, ]$Ymin.line <- carCC[j - 1, ]$Ymin.line + 1
          }
        }
      }
    }
    rm(j)
  }
  
  ## check if the angle of the first row of produced cell is high: if yes: recalculation of the Xinter.min and Xinter.max
  ### calculation of the angle of the vector
  newDif <- sqrt((carCC[1, ]$Xinter.min - carCC[nrow(carCC), ]$Xinter.max)^2 + (max(carCC$Ymin.line))^2) # norme of the vector
  ### the angle in regard to the X axis. To obtain the angle of correction, remove to 90 deg (which is the total angle)
  angle.deg.X <- acos((carCC[nrow(carCC), ]$Xinter.max - carCC[1, ]$Xinter.min) / newDif) * 180 / pi
  
  if(angle.deg.X > 5){
    carCC[1, ]$Xinter.max <- carCC[1, ]$Xinter.min + ((carCC[1, ]$perLength * (carCC[nrow(carCC), ]$Xinter.max - carCC[1, ]$Xinter.min) / 100) * cos(angle.deg.X * pi/180))
    for(j in 2:(nrow(carCC) - 1)){
      carCC[j, ]$Xinter.min <- carCC[j - 1, ]$Xinter.max
      carCC[j, ]$Xinter.max <- carCC[1, ]$Xinter.min + ((carCC[j, ]$cumperLength * (carCC[nrow(carCC), ]$Xinter.max - carCC[1, ]$Xinter.min) / 100) * cos(angle.deg.X * pi/180))
    }
    rm(j)
    carCC[nrow(carCC), ]$Xinter.min <- carCC[nrow(carCC) - 1, ]$Xinter.max
  }
  
  
  rm(fC, dif, fCmax, fCmin, Ymin, Ycell)
  
  ## Number of layer of 30 µm of width (roughly 1.5 fibers)
  carCC$layer <- (((carCC$Y.cor + cD/4) - (carCC$Ymin.line - cD/4)) / widthL) 
  
  # 3. summary of the cambial cell and first cells characteristics per image #####
  ## cell dataframe
  cell <- data.frame(id = unique(dataset$idD),
                     cFile = max(cFile), cD = cD, LgcFile = nCc)
  
  rm(cFile, cD, nCc)
  
  ##### 4. to export #####
  outC <- list(carCC, cell)
  return(outC)
  
}

#### end function ccls ##########

#################### FUNCTION sectors #####

# function to estimate the edges of each horizontal layers and then sectors of the image
# return graphic with point coordinates of the image and sectors borders
# return list of 1 object: sector dataframe

#' @dataset: dataframe with reoriented data (X.cor, Y.cor) and "type" column
#' @carCC: dataframe output from ccls function with cambial lines and layers characteristics for each images
#' @cell: dataframe output from ccls function with number of cambial line (cFile), average cambial cell diameter (cD), length of cambial layer (LgcFile), min and max X.cor values (fCmin, fCmax) of the first part of the image for each images
#' @dataIS: dataframe with coordinates of the ashape sector

sectors <- function(dataset, cFile, cell, dataIS, make.plot = TRUE){
  j         <- NULL ; l <- NULL
  yInter    <- NULL # values of the Y coordinates for each layer
  SecEdges <- NULL # dataframe summarizing horizontal layers and sectors characteristics
  ThisL     <- NULL # susbet of the dataframe SecEdges (one cambial line)
  Call      <- NULL # selected cells within the layer l
  Cmin      <- NULL # selected cells for the left edge of the layer
  Cmax      <- NULL # selected cells for the right edge of the layer
  XcorAs    <- NULL # X.cor of the dataIS cell (ind1 and ind2 number are not matching with idC)
  CaSi      <- NULL ; CaSa <- NULL ; CaS <- NULL
  
  newC      <- NULL # new cell for xLeftC or xRightC
  estX      <- NULL # model for estimation of the xLeftL or xRightL coordinates
  outC      <- NULL # output dataframe after calculation of layer X edges
  outL      <- NULL # output dataframe after calculation of layer X edges
  idSector    <- NULL # vector of sectors ids
  c         <- NULL ; newLength <- NULL
  ThisP     <- NULL # subset of the dataframe outL (one sectors)
  outP      <- NULL # output dataframe after calculation of sectors characteristics
  title     <- NULL # title of the chart 
  
  # 1. creation of vector and object from inputed dataframes #####
  carCC  <- cFile[cFile$id == unique(dataset$idD), ]
  cFile  <- unique(carCC$cFile)
  cD     <- cell[cell$id == unique(dataset$idD), ]$cD
  layer  <- c(1:max(round(carCC$layer)))
  
  # 2. creation of horizontal layers and sectors id #####
  if(max(layer) > 1){
    SecEdges <- data.frame(idD = NA, layer = rep(layer, each = max(cFile)), cFile = rep(cFile, max(layer)), idSector = NA,
                            xLeftC = NA, xRightC = NA,
                            yInfL = NA, ySupL = NA, xLeftL = NA, xRightL = NA, yLeftL = NA, yRightL = NA,
                            xLeftLed = NA, xRightLed = NA,
                            xLengthL = NA, xLengthInfL = NA, xLengthSupL = NA,
                            xLeftInfP = NA, xLeftSupP = NA, xLeftEdgeS = NA, yLeftInfS = NA, yLeftSupS = NA, yLeftEdgeS = NA,
                            xRightInfS = NA, xRightSupS = NA, xRightEdgeS = NA, yRightInfS = NA, yRightSupS = NA, yRightEdgeS = NA)
  }else{
    SecEdges <- data.frame(idD = NA, layer = rep(c(1, 2), each = max(cFile)), cFile = rep(cFile, 2), idSector = NA,
                            xLeftC = NA, xRightC = NA,
                            yInfL = NA, ySupL = NA, xLeftL = NA, xRightL = NA, yLeftL = NA, yRightL = NA,
                            xLeftLed = NA, xRightLed = NA,
                            xLengthL = NA, xLengthInfL = NA, xLengthSupL = NA,
                            xLeftInfP = NA, xLeftSupP = NA, xLeftEdgeS = NA, yLeftInfS = NA, yLeftSupS = NA, yLeftEdgeS = NA,
                            xRightInfS = NA, xRightSupS = NA, xRightEdgeS = NA, yRightInfS = NA, yRightSupS = NA, yRightEdgeS = NA)
    
    
    if(layer == 1){ layer <- seq(1, max(SecEdges$layer), 1) }
  }
  
  
  SecEdges$idD    <- unique(as.character(dataset$idD))
  SecEdges$idSector <- paste("c", SecEdges$cFile, ".l", SecEdges$layer, sep = "")
  
  ## 2a. assessment yInfL et ySupL per horizontal layers for each cambial line ####
  for(j in 1:length(cFile)){
    ThisL <- SecEdges[SecEdges$cFile == cFile[j], ]
    
    #y interval based on roughly 1.5 fibers 30 or 40 um regarding the species) and varied to follow fluctuations of first cells position
    ### recalculation of the interval to have a homogenization of values
    if(max(layer) > 1){
      newInter <- ((carCC[carCC$cFile == cFile[j], ]$Y.cor + 3) - (carCC[carCC$cFile == cFile[j], ]$Ymin.line - 5)) / max(layer)
      yInter   <- seq(carCC[carCC$cFile == cFile[j], ]$Ymin.line - 5, carCC[carCC$cFile == cFile[j], ]$Y.cor + 3, newInter)
      
      #if the effective of yInter is not equal to layer + 1
      if((length(yInter)-1) != nrow(ThisL)){
        yInter <- c(yInter, carCC[carCC$cFile == cFile[j],]$Y.cor + 3)
      }
      
    }else{
      newInter <- ((carCC[carCC$cFile == cFile[j], ]$Y.cor + 3) - (carCC[carCC$cFile == cFile[j], ]$Ymin.line - 5)) / 2
      yInter   <- seq(carCC[carCC$cFile == cFile[j], ]$Ymin.line - 5, carCC[carCC$cFile == cFile[j], ]$Y.cor + 3, newInter)
    }
    
    ### saving data
    ThisL$yInfL <- round(yInter[c(1: length(yInter)-1)], 3)
    ThisL$ySupL <- round(yInter[c(2: length(yInter))], 3)
    
    ### saving
    outC <- rbind(outC, ThisL)
    
  } #end loop j
  
  rm(j, yInter, ThisL, newInter, newC, estX)
  
  ## 2b. assessment xLeft et xRight per horizontal layers for each cambial line ####
  for(l in 1:length(layer)){
    ThisL <- outC[outC$layer == layer[l], ]
    vxE   <- NULL
    Call  <- NULL
    
    ### extraction of the potential values
    #### cells within the layer l
    if(layer[l] != 1){
      Call <- dataset[dataset$Y.cor > min(ThisL$yInfL) & dataset$Y.cor <=  max(ThisL$ySupL), ]
      vxE  <- Call[Call$type == "Vessel", ]
    }else{
      Call <- dataset[dataset$Y.cor <= max(ThisL$ySupL), ]
      vxE  <- Call[Call$type == "Vessel", ]
    }
    
    ### coordinates of the ashape within the layer
    XcorAs <- NULL
    CaS    <- NULL
    CaSi   <- rbind(dataIS[dataIS$y1 > ThisL[ThisL$cFile == 1, ]$yInfL & dataIS$y1 <=  ThisL[ThisL$cFile == 1, ]$ySupL, ],
                    dataIS[dataIS$y2 > ThisL[ThisL$cFile == 1, ]$yInfL & dataIS$y2 <=  ThisL[ThisL$cFile == 1, ]$ySupL, ])
    CaSa   <- rbind(dataIS[dataIS$y1 > ThisL[ThisL$cFile == max(ThisL$cFile), ]$yInfL & dataIS$y1 <=  ThisL[ThisL$cFile == max(ThisL$cFile), ]$ySupL, ],
                    dataIS[dataIS$y2 > ThisL[ThisL$cFile == max(ThisL$cFile), ]$yInfL & dataIS$y2 <=  ThisL[ThisL$cFile == max(ThisL$cFile), ]$ySupL, ])
    if(layer[l] != 1){
      CaSi <- CaSi[c(which(abs(CaSi$x2 - unique(outL[outL$layer == (layer[l] - 1), ]$xLeftL)) <= 20),
                     which(abs(CaSi$x1 - unique(outL[outL$layer == (layer[l] - 1), ]$xLeftL)) <= 20)), ]
      CaSa <- CaSa[c(which(abs(CaSa$x2 - unique(outL[outL$layer == (layer[l] - 1), ]$xRightL)) <= 20),
                     which(abs(CaSa$x1 - unique(outL[outL$layer == (layer[l] - 1), ]$xRightL)) <= 20)), ]
    }
    if(!all(nrow(CaSi) == 0, nrow(CaSa) == 0)){
      CaS    <- rbind(CaSi, CaSa)
      CaS    <- CaS[-which(duplicated(CaS)),]
      XcorAs <- unique(c(CaS$x1, CaS$x2))
    }
    
    ### special cases
    if(nrow(Call) == 0){
      #yLeftL and yRightL
      ThisL$yLeftL  <- mean(c(ThisL$yInfL, ThisL$ySupL))
      ThisL$yRightL <- mean(c(ThisL$yInfL, ThisL$ySupL))
      #xLeftL and yRightL
      ThisL$xLeftL  <- unique(outL[outL$layer == (layer[l]-1), ]$xLeftL)
      ThisL$xLeftC  <- unique(outL[outL$layer == (layer[l]-1), ]$xLeftC)
      ThisL$xRightL <- unique(outL[outL$layer == (layer[l]-1), ]$xRightL)
      ThisL$xRightC <- unique(outL[outL$layer == (layer[l]-1), ]$xRightC)
      
    }else{ #nrow(Call) !=0
      #case when layer between cFile are switched
      if(any(max(ThisL$yInfL) - min(ThisL$yInfL) > cD, max(ThisL$ySupL) - min(ThisL$ySupL) > cD)){
        Cmin <- dataset[dataset$Y.cor > ThisL[ThisL$cFile == 1, ]$yInfL & dataset$Y.cor <=  ThisL[ThisL$cFile == 1,]$ySupL, ]
        Cmax <- dataset[dataset$Y.cor > ThisL[ThisL$cFile == max(ThisL$cFile), ]$yInfL & dataset$Y.cor <=  ThisL[ThisL$cFile == max(ThisL$cFile),]$ySupL, ]
        
        if(all(nrow(Cmax) != 0, nrow(Cmin) != 0)){ Call <- rbind(Cmax[Cmax$X.cor >= min(Cmin$X.cor), ], Cmin[Cmin$X.cor <= max(Cmax$X.cor), ])
        }else if(nrow(Cmax) == 0){ Call <- Cmin
        }else if(nrow(Cmin) == 0){ Call <- Cmax }
        
        if(length(which(duplicated(Call$idC))) != 0){
          Call <- Call[-which(duplicated(Call$idC)),]
        }
        vxE <- Call[Call$type == "Vessel", ]
      }
      
      if(all(!is.null(XcorAs), any(Call$X.cor %in% XcorAs))){ # if idC of Call are included in dataIS #
        #### extraction of the minimal and maximal X.cor value of the layer l
        Cmin <- Call[Call$X.cor == min(Call[Call$X.cor %in% XcorAs, ]$X.cor, na.rm = TRUE), ]
        Cmax <- Call[Call$X.cor == max(Call[Call$X.cor %in% XcorAs, ]$X.cor, na.rm = TRUE), ]
        Cmin <- Cmin[which.min(Cmin$Y.cor), ]
        Cmax <- Cmax[which.min(Cmax$Y.cor), ]
        #### yLeftL and yRightL
        ThisL$yLeftL  <- Cmin$Y.cor
        ThisL$yRightL <- Cmax$Y.cor
        
      }else if(layer[l] != max(layer)){ # if idC of Call are not included in dataIS #
        #### test if there is a coordinate in the layer + 1
        Call <- dataset[dataset$Y.cor > min(outC[outC$layer == layer[l] + 1, ]$yInfL) & dataset$Y.cor <= max(outC[outC$layer == layer[l] + 1, ]$ySupL), ]
        vxE  <- Call[Call$type == "Vessel", ]
        #### test if the new values are not kicking out the potential edges of the layer l
        if(all(!is.null(XcorAs), any(Call$X.cor %in% XcorAs))){
          if(Call[Call$X.cor == min(Call[Call$X.cor %in% XcorAs, ]$X.cor, na.rm = TRUE), ]$X.cor < min(Cmin$X.cor)){
            Cmin <- Call[Call$X.cor == min(Call[Call$X.cor %in% XcorAs, ]$X.cor, na.rm = TRUE), ]
          }else{
            Cmin <- Cmin[which.min(Cmin$X.cor), ]
          }
          if(all(!is.null(XcorAs), Call[Call$X.cor == max(Call[Call$X.cor %in% XcorAs, ]$X.cor, na.rm = TRUE), ]$X.cor > max(Cmax$X.cor))){
            Cmax <- Call[Call$X.cor == max(Call[Call$X.cor %in% XcorAs, ]$X.cor, na.rm = TRUE), ]
          }else{
            Cmax <- Cmax[which.max(Cmax$X.cor), ]
          }
          ##### yLeftL and yRightL
          ThisL$yLeftL  <- mean(c(ThisL$yInfL, ThisL$ySupL))
          ThisL$yRightL <- mean(c(ThisL$yInfL, ThisL$ySupL))
          
        }else{
          #### extraction of the minimal and maximal X.cor value of the layer l based on the cell coordinates and not the alpha shape
          Cmin <- dataset[dataset$Y.cor > ThisL[ThisL$cFile == 1,]$yInfL & dataset$Y.cor <=  ThisL[ThisL$cFile == 1,]$ySupL, ]
          Cmax <- dataset[dataset$Y.cor > ThisL[ThisL$cFile == max(ThisL$cFile),]$yInfL & dataset$Y.cor <=  ThisL[ThisL$cFile == max(ThisL$cFile),]$ySupL, ]
          
          if(all(nrow(Cmax) != 0, nrow(Cmin) != 0)){ Call <- rbind(Cmax[Cmax$X.cor >= min(Cmin$X.cor), ], Cmin[Cmin$X.cor <= max(Cmax$X.cor), ])
          }else if(nrow(Cmax) == 0){ Call <- Cmin
          }else if(nrow(Cmin) == 0){ Call <- Cmax }
          
          vxE  <- Call[Call$type == "Vessel", ]
          
          if(layer[l] != 1){ Call <- Call[Call$type != "Vessel", ] }
          if(length(which(duplicated(Call$idC))) != 0){
            Call <- Call[-which(duplicated(Call$idC)),]
          }
          Cmin <- Call[which.min(Call$X.cor), ]
          Cmax <- Call[which.max(Call$X.cor), ]
          #### yLeftL and yRightL
          ThisL$yLeftL  <- mean(c(ThisL$yInfL, ThisL$ySupL))
          ThisL$yRightL <- mean(c(ThisL$yInfL, ThisL$ySupL))
        }
      }else if(layer[l] == max(layer)){
        Cmin <- Call[which.min(Call$X.cor), ]
        Cmax <- Call[which.max(Call$X.cor), ]
        #yLeftL and yRightL
        ThisL$yLeftL  <- mean(c(ThisL$yInfL, ThisL$ySupL))
        ThisL$yRightL <- mean(c(ThisL$yInfL, ThisL$ySupL))
      }
    }
    
    ### comparison with referential values
    if(nrow(Call) == 0){
      #yLeftL and yRightL
      ThisL$yLeftL  <- mean(c(ThisL$yInfL, ThisL$ySupL))
      ThisL$yRightL <- mean(c(ThisL$yInfL, ThisL$ySupL))
      #xLeftL and yRightL
      ThisL$xLeftL  <- unique(outL[outL$layer == (layer[l]-1), ]$xLeftL)
      ThisL$xLeftC  <- unique(outL[outL$layer == (layer[l]-1), ]$xLeftC)
      ThisL$xRightL <- unique(outL[outL$layer == (layer[l]-1), ]$xRightL)
      ThisL$xRightC <- unique(outL[outL$layer == (layer[l]-1), ]$xRightC)
      
    }else{
      if(layer[l] == max(layer)){
        #first test on Cmin and Cmax from ashape
        refLa <- Cmin$X.cor - carCC[carCC$cFile == 1,]$Xmin                <= carCC[carCC$cFile == 1,]$TransvDiaTheo
        refRa <- carCC[carCC$cFile == max(carCC$cFile),]$Xmax - Cmax$X.cor <= carCC[carCC$cFile == max(carCC$cFile),]$TransvDiaTheo
        #second test on Cmin and Cmax from raw data
        refLr <- min(Call$X.cor) - carCC[carCC$cFile == 1,]$Xmin                <= 20
        refRr <- carCC[carCC$cFile == max(carCC$cFile),]$Xmax - max(Call$X.cor) <= 20  
        
      }else if(all(layer[l] != 1, layer[l] != max(layer))){
        #### conditions to validate X coordinates for layer edges ##
        if(unique(dataset$TWood) == "DP"){
          ##### first test on Cmin and Cmax from ashape
          refLa <- abs(Cmin$X.cor - unique(outL[outL$layer == (layer[l]-1), ]$xLeftL))  <= 20
          refRa <- abs(Cmax$X.cor - unique(outL[outL$layer == (layer[l]-1), ]$xRightL)) <= 20
          ##### second test on Cmin and Cmax from raw data
          refLr <- abs(min(Call[Call$type != "Vessel", ]$X.cor) - unique(outL[outL$layer == (layer[l]-1), ]$xLeftL))  <= mean(c(carCC[carCC$cFile == 1, ]$TransvDiaTheo, carCC[carCC$cFile == 2, ]$TransvDiaTheo))
          refRr <- abs(max(Call[Call$type != "Vessel", ]$X.cor) - unique(outL[outL$layer == (layer[l]-1), ]$xRightL)) <= mean(c(carCC[carCC$cFile == max(carCC$cFile), ]$TransvDiaTheo, carCC[carCC$cFile == max(carCC$cFile) - 1, ]$TransvDiaTheo))
        }else{
          ##### first test on Cmin and Cmax from ashape
          refLa <- Cmin$X.cor -  unique(outL[outL$layer == (layer[l]-1), ]$xLeftL) <= 45
          refRa <- unique(outL[outL$layer == (layer[l]-1), ]$xRightL) - Cmax$X.cor <= 45
          ##### second test on Cmin and Cmax from raw data
          refLr <- min(Call[Call$type != "Vessel", ]$X.cor) - unique(outL[outL$layer == (layer[l]-1), ]$xLeftL)  <= 40
          refRr <- unique(outL[outL$layer == (layer[l]-1), ]$xRightL) - max(Call[Call$type != "Vessel", ]$X.cor) <= 40
        }
      }else if (layer[l] == 1){
        if(unique(dataset$TWood) == "DP"){
          ##### first test on Cmin and Cmax from ashape
          refLa <- Cmin$X.cor - carCC[carCC$cFile == 1,]$Xinter.min                <= carCC[carCC$cFile == 1,]$TransvDiaTheo
          refRa <- carCC[carCC$cFile == max(carCC$cFile),]$Xinter.max - Cmax$X.cor <= carCC[carCC$cFile == max(carCC$cFile),]$TransvDiaTheo
          ##### second test on Cmin and Cmax from raw data
          refLr <- min(Call$X.cor) - carCC[carCC$cFile == 1,]$Xinter.min                <= 20
          refRr <- carCC[carCC$cFile == max(carCC$cFile),]$Xinter.max - max(Call$X.cor) <= 20
          
        }else{
          ##### first test on Cmin and Cmax from ashape
          refLa <- which.min(abs(Cmin$X.cor - sort(dataset[dataset$type == "Cambial",]$X.cor))) <= round((length(dataset[dataset$type == "Cambial",]$idC) * (2/5)))
          refRa <- which.min(abs(Cmax$X.cor - sort(dataset[dataset$type == "Cambial",]$X.cor))) >= round((length(dataset[dataset$type == "Cambial",]$idC) * (4/5)))
          ##### second test on Cmin and Cmax from raw data
          refLr <- which.min(abs(Call[which.min(Call$X.cor), ]$X.cor - sort(dataset[dataset$type == "Cambial",]$X.cor))) <= round((length(dataset[dataset$type == "Cambial",]$idC) * (2/5)))
          refRr <- which.min(abs(Call[which.max(Call$X.cor), ]$X.cor - sort(dataset[dataset$type == "Cambial",]$X.cor))) >= round((length(dataset[dataset$type == "Cambial",]$idC) * (4/5)))
        }
      }
    }
    
    if(nrow(Call) != 0){
      #### xLeftL ##
      if(refLa == TRUE){
        #if detected Cmin APF is in the first third of the cambial lines
        ThisL$xLeftC <- Cmin$idC
        ThisL$xLeftL <- Cmin$X.cor
        
        #several possibilities: alpha shape too convex, vessel in the bottom corner or spatial shift
      }else if(refLr == TRUE){
        #in the case of Cmin from alpha shape: test on the raw data
        #if the raw data is ok: too convex alpha shape
        ThisL$xLeftC <- Call[which.min(Call$X.cor), ]$idC
        ThisL$xLeftL <- Call[which.min(Call$X.cor), ]$X.cor
        
      }else if(layer[l] == 1){
        ##### case of the first layer with potentially a vessel as the first cell
        ThisL$xLeftL <- carCC[carCC$cFile == 1, ]$Xinter.min
        
      }else{
        #if the raw data is not ok: spatial shift
        if(all(any(!is.null(CaS), nrow(CaS) != 0),
               any(abs(unique(outL[outL$layer == (layer[l] - 1), ]$xLeftL) - CaS[which.min(abs(CaS$x1 - unique(outL[outL$layer == (layer[l] - 1), ]$xLeftL))), ]$x1) <= 20,
                   abs(unique(outL[outL$layer == (layer[l] - 1), ]$xLeftL) - CaS[which.min(abs(CaS$x2 - unique(outL[outL$layer == (layer[l] - 1), ]$xLeftL))), ]$x2) <= 20))){
          newC <- CaS[which.min(abs(CaS$x1 - unique(outL[outL$layer == (layer[l] - 1), ]$xLeftL))), ]
          estX <- lm(c(newC$y1, newC$y2) ~ c(newC$x1, newC$x2))
          ThisL$xLeftL <- (unique(ThisL$yLeftL) - as.numeric(estX$coefficients[1])) / as.numeric(estX$coefficients[2])
          ThisL$xLeftC <- unique(outL[outL$layer == (layer[l]-1), ]$xLeftC)
        }else{
          newC    <- rbind(dataIS[which(dataIS$y1 < unique(ThisL$yLeftL) & dataIS$y2 > unique(ThisL$yLeftL)), ],
                           dataIS[which(dataIS$y1 > unique(ThisL$yLeftL) & dataIS$y2 < unique(ThisL$yLeftL)), ])
          newC$n1 <- sqrt((newC$x1 - unique(outL[outL$layer == (layer[l] - 1), ]$xLeftL))^2 + (newC$y1 - unique(ThisL$yLeftL))^2)
          newC$n2 <- sqrt((newC$x2 - unique(outL[outL$layer == (layer[l] - 1), ]$xLeftL))^2 + (newC$y2 - unique(ThisL$yLeftL))^2)
          newC    <- rbind(newC[which.min(newC$n1), ], newC[which.min(newC$n2), ])
          
          if(newC[which.min(newC$n1), ]$n1 < newC[which.min(newC$n2), ]$n2){
            newC <- newC[which.min(newC$n1), ]
          }else{
            newC <- newC[which.min(newC$n2), ]
          }
          #model
          if(unique(dataset$TWood == "DP")){                                    #spatial shift in ashape are moderate in DP species
            #estimated x coordinate
            if(identical(round(newC$x1, 3), round(newC$x2, 3))){
              ThisL$xLeftL <- newC$x2
            }else{
              estX <- lm(c(newC$y1, newC$y2) ~ c(newC$x1, newC$x2))
              ThisL$xLeftL <- (unique(ThisL$yLeftL) - as.numeric(estX$coefficients[1])) / as.numeric(estX$coefficients[2]) }
          }else{                                                                #spatial shift in ashape can be big in RP species
            if(unique(outL[outL$layer == (layer[l]-1), ]$xLeftL) %in% c(newC$x1, newC$x2)){
              ThisL$xLeftL <- unique(outL[outL$layer == (layer[l]-1), ]$xLeftL)
            }else if(newC$y1 > newC$y2){
              estX <- lm(c(newC$y1, unique(outL[outL$layer == (layer[l]-1), ]$yLeftL)) ~ c(newC$x1, unique(outL[outL$layer == (layer[l]-1), ]$xLeftL)))
              #estimated x coordinate 
              ThisL$xLeftL <- (unique(ThisL$yLeftL) - as.numeric(estX$coefficients[1])) / as.numeric(estX$coefficients[2])
            }else if(newC$x1 == newC$x2){
              ThisL$xLeftL <- newC$x2
            }else{
              estX <- lm(c(newC$y2, unique(outL[outL$layer == (layer[l]-1), ]$yLeftL)) ~ c(newC$x2, unique(outL[outL$layer == (layer[l]-1), ]$xLeftL)))
              #estimated x coordinate 
              ThisL$xLeftL <- (unique(ThisL$yLeftL) - as.numeric(estX$coefficients[1])) / as.numeric(estX$coefficients[2])
            }
          }
          #idC
          ThisL$xLeftC <- unique(outC[outC$layer == (layer[l]-1), ]$xLeftC)
        }
      }
      
      #### xLeftL ##
      if(refRa == TRUE){
        #if detected Cmax APF is in the last third of the cambial lines
        ThisL$xRightC <- Cmax$idC
        ThisL$xRightL <- Cmax$X.cor
        
        #several possibilities: alpha shape too convex, vessel in the bottom corner or spatial shift
      }else if(refRr == TRUE){
        #in the case of Cmax from alpha shape: test on the Cmax from raw data
        #if the raw data is ok: too convex alpha shape
        ThisL$xRightC <- Call[which.max(Call$X.cor), ]$idC
        ThisL$xRightL <- Call[which.max(Call$X.cor), ]$X.cor
        
      }else if(layer[l] == 1){
        ##### case of the first layer with potentially a vessel as the first cell
        ThisL$xRightL <- carCC[carCC$cFile == max(carCC$cFile),]$Xinter.max 
        
      }else{
        #if the raw data is not ok: spatial shift
        if(all(any(!is.null(CaS), nrow(CaS) != 0),
               any(abs(unique(outL[outL$layer == (layer[l] - 1), ]$xRightL) - CaS[which.min(abs(CaS$x1 - unique(outL[outL$layer == (layer[l] - 1), ]$xRightL))), ]$x1) <= 20,
                   abs(unique(outL[outL$layer == (layer[l] - 1), ]$xRightL) - CaS[which.min(abs(CaS$x2 - unique(outL[outL$layer == (layer[l] - 1), ]$xRightL))), ]$x2) <= 20))){
          newC <- CaS[which.min(abs(CaS$x1 - unique(outL[outL$layer == (layer[l] - 1), ]$xRightL))), ]
          estX <- lm(c(newC$y1, newC$y2) ~ c(newC$x1, newC$x2))
          ThisL$xRightL <- ((unique(ThisL$yRightL) - as.numeric(estX$coefficients[1])) / as.numeric(estX$coefficients[2]))
          ThisL$xRightC <- unique(outL[outL$layer == (layer[l]-1), ]$xRightC)
        }else{
          newC    <- rbind(dataIS[which(dataIS$y1 < unique(ThisL$yRightL) & dataIS$y2 > unique(ThisL$yRightL)), ],
                           dataIS[which(dataIS$y1 > unique(ThisL$yRightL) & dataIS$y2 < unique(ThisL$yRightL)), ])
          newC$n1 <- sqrt((newC$x1 - unique(outL[outL$layer == (layer[l] - 1), ]$xRightL))^2 + (newC$y1 - unique(ThisL$yRightL))^2)
          newC$n2 <- sqrt((newC$x2 - unique(outL[outL$layer == (layer[l] - 1), ]$xRightL))^2 + (newC$y2 - unique(ThisL$yRightL))^2)
          newC    <- rbind(newC[which.min(newC$n1), ], newC[which.min(newC$n2), ])
          
          if(newC[which.min(newC$n1), ]$n1 < newC[which.min(newC$n2), ]$n2){
            newC <- newC[which.min(newC$n1), ]
          }else{
            newC <- newC[which.min(newC$n2), ]
          }
          if(unique(dataset$TWood == "DP")){                                    #spatial shift in ashape are moderate in DP species
            #estimated x coordinate
            if(identical(round(newC$x1, 3), round(newC$x2, 3))){
              ThisL$xRightL <- newC$x2
            }else{
              estX <- lm(c(newC$y1, newC$y2) ~ c(newC$x1, newC$x2))
              ThisL$xRightL <- ((unique(ThisL$yRightL) - as.numeric(estX$coefficients[1])) / as.numeric(estX$coefficients[2])) }
          }else{                                                      # spatial shift in ashape can be big in RP species
            if(unique(outL[outL$layer == (layer[l]-1), ]$xRightL) %in% c(newC$x1, newC$x2)){
              ThisL$xRightL <- unique(outL[outL$layer == (layer[l]-1), ]$xRightL)
            }else if(newC$y1 > newC$y2){
              estX <- lm(c(newC$y1, unique(outL[outL$layer == (layer[l]-1), ]$yRightL)) ~ c(newC$x1, unique(outL[outL$layer == (layer[l]-1), ]$xRightL)))
              #estimated x coordinate 
              ThisL$xRightL <- ((unique(ThisL$yRightL) - as.numeric(estX$coefficients[1])) / as.numeric(estX$coefficients[2]))
            }else if(newC$x1 == newC$x2){
              ThisL$xRightL <- newC$x2
            }else{
              estX <- lm(c(newC$y2, unique(outL[outL$layer == (layer[l]-1), ]$yRightL)) ~ c(newC$x2, unique(outL[outL$layer == (layer[l]-1), ]$xRightL)))
              ThisL$xRightL <- ((unique(ThisL$yRightL) - as.numeric(estX$coefficients[1])) / as.numeric(estX$coefficients[2]))
            }
          }
          ThisL$xRightC <- unique(outL[outL$layer == (layer[l]-1), ]$xRightC)
        }
      }
    }
    
    ### if there is an edged vessel, test on the edge X ##
    if(nrow(vxE) != 0){
      if(min(vxE$X.cor) < unique(ThisL$xLeftL)){
        ThisL$xLeftL <- min(vxE$X.cor)
        ThisL$xLeftC <- vxE[which.min(vxE$X.cor), ]$idC
      }
    }
    
    if(nrow(vxE) != 0){
      if(max(vxE$X.cor) > unique(ThisL$xRightL)){
        ThisL$xRightL <- max(vxE$X.cor)
        ThisL$xRightC <- vxE[which.max(vxE$X.cor), ]$idC
      }
    }
    
    ### comparison of results for the first layer ##
    if(layer[l] == 1){
      if(carCC[carCC$cFile == 1, ]$Xinter.min < unique(ThisL$xLeftL)){
        ThisL$xLeftL <- carCC[carCC$cFile == 1, ]$Xinter.min
      }
      
      if(carCC[carCC$cFile == max(carCC$cFile), ]$Xinter.max > unique(ThisL$xRightL)){
        ThisL$xRightL <-carCC[carCC$cFile == max(carCC$cFile), ]$Xinter.max
      }
    }
    
    ### saving
    outL <- rbind(outL, ThisL)
    
  } # end loop l
  
  rm(l)
  rm(ThisL, Call, Cmin, Cmax)
  
  # 3. coordinates of each sectors #####
  ## addition of the half of the tranversal diameter
  #test with only m um of margin around cells
  m <- ifelse(unique(dataset$TWood == "DP"), 10, 5)
  s <- ifelse(cD/2 > 10, 7, cD/2)
  outL$xLeftLed  <- outL$xLeftL - s
  outL$xRightLed <- outL$xRightL + s
  
  ## layer length/norm
  outL$xLengthL <- sqrt((outL$xLeftLed - outL$xRightLed)^2 + (outL$yLeftL - outL$yRightL)^2)
  
  for(j in 1:length(layer)){
    ThisP <- outL[outL$layer == layer[j],]
    
    ### ## case of the first layer ## L1 ####  
    if(unique(ThisP$layer) == 1){
      #### newLength of the layer used to recalculate the transversal layer width of each cFile ##
      #to measure newLentgh used for SupP values, need ThisP[ThisP$cFile == 1,]$x- and y-LeftSupP and ThisP[ThisP$cFile == max(cFile),]$x- and y-RightSupP
      ##### C1
      ThisP[ThisP$cFile == 1, ]$xLeftInfP   <- ThisP[ThisP$cFile == 1, ]$xLeftLed
      if(all(max(outL$layer) > 1, outL[outL$idSector == "c1.l2", ]$xLeftLed <= ThisP[ThisP$cFile == 1, ]$xLeftLed + carCC[carCC$cFile == 1, ]$TransvDiaTheo)){
        ThisP[ThisP$cFile == 1, ]$xLeftSupP <- min(outL[outL$idSector == "c1.l2", ]$xLeftLed, ThisP[ThisP$cFile == 1, ]$xLeftLed)
      }else{
        ThisP[ThisP$cFile == 1, ]$xLeftSupP <- ThisP[ThisP$cFile == 1, ]$xLeftLed
      }
      
      ThisP[ThisP$cFile == 1, ]$yLeftSupS   <- ThisP[ThisP$cFile == 1, ]$ySupL
      ThisP[ThisP$cFile == 1, ]$yLeftInfS   <- ThisP[ThisP$cFile == 1, ]$yInfL - m
      
      ##### Cmax
      ThisP[ThisP$cFile == max(cFile), ]$xRightInfS   <- ThisP[ThisP$cFile == max(cFile), ]$xRightLed
      if(all(max(outL$layer) > 1, outL[outL$idSector == paste("c", max(cFile), ".l2", sep = ""), ]$xRightLed >= ThisP[ThisP$cFile == max(cFile), ]$xRightInfS - carCC[carCC$cFile == max(carCC$cFile), ]$TransvDiaTheo)){
        ThisP[ThisP$cFile == max(cFile), ]$xRightSupS <- max(outL[outL$idSector == paste("c", max(cFile), ".l2", sep = ""), ]$xRightLed, ThisP[ThisP$cFile == max(cFile), ]$xRightLed)
      }else{
        ThisP[ThisP$cFile == max(cFile), ]$xRightSupS <- ThisP[ThisP$cFile == max(cFile), ]$xRightLed
      }
      ThisP[ThisP$cFile == max(cFile), ]$yRightSupS   <- ThisP[ThisP$cFile == max(cFile), ]$ySupL
      
      ##### newLength
      ThisP$xLengthSupL <- sqrt((ThisP[ThisP$cFile == 1, ]$xLeftSupP - ThisP[ThisP$cFile == max(cFile), ]$xRightSupS)^2 + (ThisP[ThisP$cFile == 1, ]$yLeftSupS - ThisP[ThisP$cFile == max(cFile),]$yRightSupS)^2)
      ThisP$xLengthInfL <- sqrt((ThisP[ThisP$cFile == 1, ]$xLeftInfP - ThisP[ThisP$cFile == max(cFile), ]$xRightInfS)^2 + (ThisP[ThisP$cFile == 1, ]$yInfL - ThisP[ThisP$cFile == max(cFile),]$yInfL)^2)
      
      #### y Right ##
      ##### yRightInfS and yRightSupS
      for(c in 1:(max(cFile)-1)){
        ThisP[ThisP$cFile == c, ]$yRightInfS <- min(c(outL[outL$cFile == (c + 1) & outL$layer == ThisP$layer, ]$yInfL, ThisP[ThisP$cFile == c, ]$yInfL)) - m
        ThisP[ThisP$cFile == c, ]$yRightSupS <- mean(c(ThisP[ThisP$cFile == c, ]$ySupL, ThisP[ThisP$cFile == c + 1, ]$ySupL))
      }
      rm(c)
      
      if(nrow(ThisP[ThisP$yRightL == 0,]) != 0 ){ ThisP[ThisP$yRightL == 0,]$yRightL <- mean(ThisP$yLeftL) } # in case of vessel bias and so yRightL = 0
      ThisP[ThisP$cFile == max(cFile),]$yRightInfS <- min(ThisP[ThisP$cFile == max(cFile),]$yInfL, ThisP[ThisP$cFile == max(cFile) - 1,]$yInf)
      
      ##### yRightEdgeS
      for(c in 1:(max(cFile) - 1)){
        ThisP[ThisP$cFile == c, ]$yRightEdgeS <-  mean(c(ThisP[ThisP$cFile == c + 1, ]$yInfL, ThisP[ThisP$cFile == c + 1, ]$ySupL))
      }
      rm(c)
      
      if(ThisP[ThisP$cFile == max(cFile), ]$yRightL >= ThisP[ThisP$cFile == max(cFile), ]$yRightSupS){
        ThisP[ThisP$cFile == max(cFile), ]$yRightEdgeS <- mean(c(ThisP[ThisP$cFile == max(cFile), ]$yInfL, ThisP[ThisP$cFile == max(cFile), ]$ySupL))
      }else{
        if(ThisP[ThisP$cFile == max(cFile), ]$yRightL <= ThisP[ThisP$cFile == max(cFile), ]$yRightInfS){
          ThisP[ThisP$cFile == max(cFile), ]$yRightEdgeS <- mean(c(ThisP[ThisP$cFile == max(cFile), ]$yInfL, ThisP[ThisP$cFile == max(cFile), ]$ySupL))
        }else{
          ThisP[ThisP$cFile == max(cFile), ]$yRightEdgeS <- ThisP[ThisP$cFile == max(cFile), ]$yRightL
        }
      } 
      
      #### y Left ##
      if(nrow(ThisP[ThisP$yLeftL == 0, ]) != 0 ){ ThisP[ThisP$yLeftL == 0, ]$yLeftL <- mean(ThisP$yRightL) }
      
      if(any(ThisP[ThisP$cFile == 1, ]$yLeftL <= ThisP[ThisP$cFile == 1, ]$yLeftInfS, ThisP[ThisP$cFile == 1, ]$yLeftL >= ThisP[ThisP$cFile == 1, ]$yLeftSupS)){
        ThisP[ThisP$cFile == 1, ]$yLeftEdgeS <- mean(c(ThisP[ThisP$cFile == 1, ]$yInfL, ThisP[ThisP$cFile == 1, ]$ySupL))
      }else{
        ThisP[ThisP$cFile == 1, ]$yLeftEdgeS <- ThisP[ThisP$cFile == 1, ]$yLeftL
      }
      
      for(c in 2:max(cFile)){
        ThisP[ThisP$cFile == c, ]$yLeftInfS  <- min(c(ThisP[ThisP$cFile == c - 1, ]$yRightInfS, ThisP[ThisP$cFile == c, ]$yInfL - m), na.rm = TRUE) 
        ThisP[ThisP$cFile == c, ]$yLeftEdgeS <- ThisP[ThisP$cFile == c - 1, ]$yRightEdgeS
        ThisP[ThisP$cFile == c, ]$yLeftSupS  <- ThisP[ThisP$cFile == c - 1, ]$yRightSup
      }
      rm(c)
      
      #### x ##
      ##### L1.C1
      ThisP[ThisP$cFile == 1, ]$xRightSupS  <- (ThisP[ThisP$cFile == 1, ]$xLeftSupP) + (carCC[carCC$cFile == 1, ]$cumperLength * ThisP[ThisP$cFile == 1, ]$xLengthSupL / 100)
      ThisP[ThisP$cFile == 1, ]$xRightInfS  <- (ThisP[ThisP$cFile == 1, ]$xLeftInfP) + (carCC[carCC$cFile == 1, ]$cumperLength * ThisP[ThisP$cFile == 1, ]$xLengthInfL / 100)
      
      if(all(ThisP[ThisP$cFile == 1, ]$xLeftLed >= ThisP[ThisP$cFile == 1, ]$xLeftSupP, ThisP[ThisP$cFile == 1, ]$xLeftLed >= ThisP[ThisP$cFile == 1, ]$xLeftInfP)){
        ThisP[ThisP$cFile == 1, ]$xLeftEdgeS <- mean(c(ThisP[ThisP$cFile == 1, ]$xLeftSupP, ThisP[ThisP$cFile == 1, ]$xLeftInfP))
      }else{
        ThisP[ThisP$cFile == 1, ]$xLeftEdgeS <- ThisP[ThisP$cFile == 1, ]$xLeftLed
      }
      
      ##### L1.C middle
      for(c in 2:(max(cFile)-1)){
        ThisP[ThisP$cFile == c, ]$xLeftSupP   <- ThisP[ThisP$cFile == c - 1, ]$xRightSupS
        ThisP[ThisP$cFile == c, ]$xLeftInfP   <- ThisP[ThisP$cFile == c - 1, ]$xRightInfS
        ThisP[ThisP$cFile == c, ]$xRightSupS  <- (ThisP[ThisP$cFile == 1, ]$xLeftSupP) + (carCC[carCC$cFile == c, ]$cumperLength * ThisP[ThisP$cFile == c, ]$xLengthSupL / 100)
        ThisP[ThisP$cFile == c, ]$xRightInfS  <- (ThisP[ThisP$cFile == 1, ]$xLeftInfP) + (carCC[carCC$cFile == c, ]$cumperLength * ThisP[ThisP$cFile == c, ]$xLengthInfL / 100)
      }
      rm(c)
      
      ##### L1.Cmax
      ThisP[ThisP$cFile == max(cFile), ]$xLeftSupP   <- ThisP[ThisP$cFile == (max(cFile) - 1), ]$xRightSupS
      ThisP[ThisP$cFile == max(cFile), ]$xLeftInfP   <- ThisP[ThisP$cFile == (max(cFile) - 1), ]$xRightInfS
      
      #check if the angle of the first row of produced cell is high: if yes: recalculation of the Xinter.min and Xinter.max
      newLength <- sqrt((ThisP[ThisP$cFile == max(cFile), ]$xRightSupS - ThisP[ThisP$cFile == 1, ]$xLeftSupP)^2 +
                          (ThisP[ThisP$cFile == max(cFile), ]$yRightSupS - ThisP[ThisP$cFile == 1, ]$yLeftSupS)^2)
      #the angle in regard to the X axis. To obtain the angle of correction, remove to 90 deg (which is the total angle)
      angle.deg.X <- acos((ThisP[ThisP$cFile == max(cFile), ]$xRightSupS - ThisP[ThisP$cFile == 1, ]$xLeftSupP) / newLength) * 180 / pi
      if(angle.deg.X > 5){
        ThisP[ThisP$cFile == 1, ]$xRightSupS   <- (ThisP[ThisP$cFile == 1, ]$xLeftSupP) + ((carCC[1, ]$perLength * (ThisP[ThisP$cFile == max(cFile), ]$xRightSupS - ThisP[ThisP$cFile == 1, ]$xLeftSupP) / 100) * cos(angle.deg.X * pi/180))
        for(c in 2:(max(cFile) - 1)){
          ThisP[ThisP$cFile == c, ]$xLeftSupP  <- ThisP[ThisP$cFile == c - 1, ]$xRightSupS
          ThisP[ThisP$cFile == c, ]$xRightSupS <- (ThisP[ThisP$cFile == 1, ]$xLeftSupP) + ((carCC[c, ]$cumperLength * (ThisP[ThisP$cFile == max(cFile), ]$xRightSupS - ThisP[ThisP$cFile == 1, ]$xLeftSupP) / 100) * cos(angle.deg.X * pi/180))
        }
        rm(c)
        ThisP[ThisP$cFile == max(cFile), ]$xLeftSupP   <- ThisP[ThisP$cFile == (max(cFile) - 1), ]$xRightSupS
      } #end test on angle.deg
      
      #### x Edges  ##
      ThisP[ThisP$cFile == 1, ]$xRightEdgeS <- mean(c(ThisP[ThisP$cFile == 1, ]$xRightSupS, ThisP[ThisP$cFile == 1, ]$xRightInfS))
      
      for(c in 2:(max(cFile) - 1)){
        ThisP[ThisP$cFile == c, ]$xLeftEdgeS  <- ThisP[ThisP$cFile == c - 1, ]$xRightEdgeS
        ThisP[ThisP$cFile == c, ]$xRightEdgeS <- (ThisP[ThisP$cFile == c, ]$xRightSupS + ThisP[ThisP$cFile == c, ]$xRightInfS) / 2
      }
      rm(c)
      
      ThisP[ThisP$cFile == max(cFile), ]$xLeftEdgeS  <- ThisP[ThisP$cFile == (max(cFile) - 1), ]$xRightEdgeS
      if(all(ThisP[ThisP$cFile == max(cFile), ]$xRightLed < ThisP[ThisP$cFile == max(cFile), ]$xRightSupS, ThisP[ThisP$cFile == max(cFile), ]$xRightLed < ThisP[ThisP$cFile == max(cFile), ]$RightInfP)){
        ThisP[ThisP$cFile == max(cFile), ]$xRightEdgeS <- mean(c(ThisP[ThisP$cFile == max(cFile), ]$xRightSupS, ThisP[ThisP$cFile == max(cFile), ]$xRightInfS))
      }else{
        ThisP[ThisP$cFile == max(cFile), ]$xRightEdgeS <- ThisP[ThisP$cFile == max(cFile), ]$xRightLed
      }
      
      ### ## case of the last layer : Lmax ####
    }else{
      if(unique(ThisP$layer) == max(layer)){
        #### y Right ##
        for(c in 1:max(cFile)){
          ThisP[ThisP$cFile == c, ]$yRightInfS <- outP[outP$idSector == paste("c", c, ".l", layer[j]-1, sep = ""), ]$yRightSupS
        }
        rm(c)
        
        for(c in 1:(max(cFile)-1)){
          ThisP[ThisP$cFile == c, ]$yRightEdgeS <- carCC[carCC$cFile == c + 1, ]$Y.cor
          ThisP[ThisP$cFile == c, ]$yRightSupS  <- max(ThisP[ThisP$cFile == c, ]$ySupL, ThisP[ThisP$cFile == c + 1, ]$ySupL) + m # + m to be sure that the cambial cell is included to the sectors
        }
        rm(c)
        
        ThisP[ThisP$cFile == max(cFile), ]$yRightEdgeS <- carCC[carCC$cFile == max(cFile), ]$Y.cor
        ThisP[ThisP$cFile == max(cFile), ]$yRightSupS  <- ThisP[ThisP$cFile == max(cFile), ]$ySupL + m
        
        #### y Left ##
        ThisP[ThisP$cFile == 1, ]$yLeftSupS  <- ThisP[ThisP$cFile == 1, ]$ySupL + m      
        ThisP[ThisP$cFile == 1, ]$yLeftInfS  <- outP[outP$idSector == paste("c1.l", layer[j]-1, sep = ""), ]$yLeftSupS
        ThisP[ThisP$cFile == 1, ]$yLeftEdgeS <- carCC[carCC$cFile == 1, ]$Y.cor
        
        for(c in 2:max(cFile)){
          ThisP[ThisP$cFile == c, ]$yLeftSupS  <- ThisP[ThisP$cFile == c - 1, ]$yRightSupS
          ThisP[ThisP$cFile == c, ]$yLeftInfS  <- ThisP[ThisP$cFile == c - 1, ]$yRightInfS
          ThisP[ThisP$cFile == c, ]$yLeftEdgeS <- ThisP[ThisP$cFile == c - 1, ]$yRightEdgeS
        }
        rm(c)
        
        #### x ##
        ##### Lmax.C1
        ThisP[ThisP$cFile == 1, ]$xLeftInfP   <- outP[outP$idSector == paste("c1.l", layer[j]-1, sep = ""), ]$xLeftSupP
        ThisP[ThisP$cFile == 1, ]$xLeftSupP   <- carCC[carCC$cFile == 1, ]$Xmin
        ThisP[ThisP$cFile == 1, ]$xLeftEdgeS  <- min(carCC[carCC$cFile == 1, ]$Xmin, ThisP[ThisP$cFile == 1, ]$xLeftLed)
        ThisP[ThisP$cFile == 1, ]$xRightInfS  <- outP[outP$idSector == paste("c1.l", layer[j]-1, sep = ""), ]$xRightSupS
        ThisP[ThisP$cFile == 1, ]$xRightSupS  <- carCC[carCC$cFile == 1, ]$Xmax
        ThisP[ThisP$cFile == 1, ]$xRightEdgeS <- carCC[carCC$cFile == 1, ]$Xmax
        
        #test on the x coordinates of the Right side in regard to the cambial cell position
        if(ThisP[ThisP$cFile == 1, ]$xRightInfS >= carCC[carCC$cFile == 2, ]$X.cor){
          ThisP[ThisP$cFile == 1, ]$xRightInfS <- carCC[carCC$cFile == 1, ]$Xmax
          outP[outP$idSector == paste("c1.l", layer[j]-1, sep = ""), ]$xRightSupS <- carCC[carCC$cFile == 1, ]$Xmax
          outP[outP$idSector == paste("c2.l", layer[j]-1, sep = ""), ]$xLeftSupP  <- carCC[carCC$cFile == 1, ]$Xmax
        }
        
        ##### Lmax.C middle
        for(c in 2:(max(cFile)-1)){
          ThisP[ThisP$cFile == c, ]$xLeftInfP   <- ThisP[ThisP$cFile == c - 1, ]$xRightInfS
          ThisP[ThisP$cFile == c, ]$xLeftSupP   <- ThisP[ThisP$cFile == c - 1, ]$xRightSupS
          ThisP[ThisP$cFile == c, ]$xLeftEdgeS  <- ThisP[ThisP$cFile == c - 1, ]$xRightEdgeS
          ThisP[ThisP$cFile == c, ]$xRightSupS  <- carCC[carCC$cFile == c, ]$Xmax
          ThisP[ThisP$cFile == c, ]$xRightInfS  <- outP[outP$idSector == paste("c", c, ".l", layer[j]-1, sep = ""), ]$xRightSupS
          ThisP[ThisP$cFile == c, ]$xRightEdgeS <- carCC[carCC$cFile == c, ]$Xmax
          
          if(ThisP[ThisP$cFile == c, ]$xRightInfS >= carCC[carCC$cFile == c + 1, ]$X.cor){
            ThisP[ThisP$cFile == c, ]$xRightInfS <- carCC[carCC$cFile == c, ]$Xmax
            outP[outP$idSector == paste("c", c, ".l", layer[j]-1, sep = ""), ]$xRightSupS    <- carCC[carCC$cFile == c, ]$Xmax
            outP[outP$idSector == paste("c", c + 1, ".l", layer[j]-1, sep = ""), ]$xLeftSupP <- carCC[carCC$cFile == c, ]$Xmax
          }
        }
        rm(c)
        
        ##### Lmax.Cmax
        ThisP[ThisP$cFile == max(cFile), ]$xLeftSupP   <- ThisP[ThisP$cFile == (max(cFile) - 1), ]$xRightSupS
        ThisP[ThisP$cFile == max(cFile), ]$xLeftEdgeS  <- ThisP[ThisP$cFile == (max(cFile) - 1), ]$xRightEdgeS
        ThisP[ThisP$cFile == max(cFile), ]$xLeftInfP   <- ThisP[ThisP$cFile == (max(cFile) - 1), ]$xRightInfS
        ThisP[ThisP$cFile == max(cFile), ]$xRightSupS  <- carCC[carCC$cFile == max(cFile), ]$Xmax
        ThisP[ThisP$cFile == max(cFile), ]$xRightInfS  <- outP[outP$idSector == paste("c", max(cFile), ".l", layer[j]-1, sep = ""), ]$xRightSupS
        ThisP[ThisP$cFile == max(cFile), ]$xRightEdgeS <- max(carCC[carCC$cFile == max(cFile), ]$Xmax, ThisP[ThisP$cFile == max(cFile), ]$xRightLed)
        
        #####
        ### ## case of middle layers : L ####
      }else{
        #### newLength ##
        #to measure newLentgh used for SupP values, need ThisP[ThisP$cFile == 1,]$x- and y-LeftSupP and ThisP[ThisP$cFile == max(cFile),]$x- and y-RightSupP
        ##### C1
        ThisP[ThisP$cFile == 1, ]$xLeftInfP   <- outP[outP$idSector == paste("c1.l", layer[j]-1, sep = ""), ]$xLeftSupP
        ThisP[ThisP$cFile == 1, ]$xRightInfS  <- outP[outP$idSector == paste("c1.l", layer[j]-1, sep = ""), ]$xRightSupS
        
        if(outL[outL$idSector == paste("c1.l", layer[j]+1, sep = ""), ]$xLeftLed <= ThisP[ThisP$cFile == 1, ]$xLeftInfP + carCC[carCC$cFile == 1, ]$TransvDiaTheo){
          ThisP[ThisP$cFile == 1, ]$xLeftSupP <- min(outL[outL$idSector == paste("c1.l", layer[j]+1, sep = ""), ]$xLeftLed, ThisP[ThisP$cFile == 1, ]$xLeftLed)
        }else{
          ThisP[ThisP$cFile == 1, ]$xLeftSupP <- mean(c(ThisP[ThisP$cFile == 1, ]$xLeftInfP, ThisP[ThisP$cFile == 1, ]$xLeftLed))
        }
        
        ThisP[ThisP$cFile == 1, ]$yLeftSupS   <- ThisP[ThisP$cFile == 1, ]$ySupL
        
        ##### Cmax
        ThisP[ThisP$cFile == max(cFile), ]$xRightEdgeS <- ThisP[ThisP$cFile == max(cFile), ]$xRightLed
        ThisP[ThisP$cFile == max(cFile), ]$xRightInfS  <- outP[outP$idSector == paste("c", max(cFile), ".l", layer[j]-1, sep = ""), ]$xRightSupS
        
        if(outL[outL$idSector == paste("c", max(cFile), ".l", layer[j]+1, sep = ""), ]$xRightLed > ThisP[ThisP$cFile == max(cFile), ]$xRightEdgeS - carCC[carCC$cFile == max(carCC$cFile), ]$TransvDiaTheo){
          ThisP[ThisP$cFile == max(cFile), ]$xRightSupS <- max(outL[outL$idSector == paste("c", max(cFile), ".l", layer[j]+1, sep = ""), ]$xRightLed,
                                                               ThisP[ThisP$cFile == max(cFile), ]$xRightLed)
        }else{
          ThisP[ThisP$cFile == max(cFile), ]$xRightSupS <- mean(c(ThisP[ThisP$cFile == max(cFile), ]$xRightEdgeS, ThisP[ThisP$cFile == max(cFile), ]$xRightLed))
        }
        
        ThisP[ThisP$cFile == max(cFile), ]$yRightSupS   <- ThisP[ThisP$cFile == max(cFile), ]$ySupL
        
        ##### y Right ##
        ThisP[ThisP$cFile != max(cFile), ]$yRightEdgeS <- (ThisP[ThisP$cFile != max(cFile), ]$yRightSupS + ThisP[ThisP$cFile != max(cFile), ]$yRightInfS) / 2
        ThisP[ThisP$cFile == max(cFile), ]$yRightInfS  <- ThisP[ThisP$cFile == max(cFile), ]$yInfL
        
        #when the detected cell for the edge coordinate is from the layer-1, the y value needs to be corrected
        if(any(ThisP[ThisP$cFile == max(cFile), ]$yRightL < ThisP[ThisP$cFile == max(cFile), ]$yRightInfS, ThisP[ThisP$cFile == max(cFile), ]$yRightL > ThisP[ThisP$cFile == max(cFile), ]$yRightSupS)){
          ThisP[ThisP$cFile == max(cFile), ]$yRightEdgeS <- mean(c(ThisP[ThisP$cFile == max(cFile), ]$yRightInfS, ThisP[ThisP$cFile == max(cFile), ]$yRightSupS))
        }else{
          ThisP[ThisP$cFile == max(cFile), ]$yRightEdgeS <- ThisP[ThisP$cFile == max(cFile), ]$yRightL
        }
        
        for(c in 1:(max(cFile)-1)){
          ThisP[ThisP$cFile == c, ]$yRightInfS <- outP[outP$idSector == paste("c", c,".l", layer[j]-1, sep = ""), ]$yRightSupS
          ThisP[ThisP$cFile == c, ]$yRightSupS <- outL[outL$idSector == paste("c", c+1,".l", layer[j], sep = ""), ]$ySupL
        }
        rm(c)
        
        ThisP[ThisP$cFile != max(cFile), ]$yRightEdgeS <- (ThisP[ThisP$cFile != max(cFile), ]$yRightSupS + ThisP[ThisP$cFile != max(cFile), ]$yRightInfS) / 2
        
        ##### y Left ##
        ThisP[ThisP$cFile == 1, ]$yLeftInfS    <- outP[outP$idSector == paste("c1.l", layer[j]-1, sep = ""),]$yLeftSupS
        
        #when the detected cell for the edge coordinate is from the layer-1, the y value needs to be corrected
        if(any(ThisP[ThisP$cFile == 1, ]$yLeftL < ThisP[ThisP$cFile == 1, ]$yLeftInfS, ThisP[ThisP$cFile == 1, ]$yLeftL > ThisP[ThisP$cFile == 1, ]$yLeftSupS)){
          ThisP[ThisP$cFile == 1, ]$yLeftEdgeS <- mean(c(ThisP[ThisP$cFile == 1, ]$yLeftInfS, ThisP[ThisP$cFile == 1, ]$yLeftSupS))
        }else{
          ThisP[ThisP$cFile == 1, ]$yLeftEdgeS <- ThisP[ThisP$cFile == 1, ]$yLeftL
        }
        
        for(c in 2:max(cFile)){
          ThisP[ThisP$cFile == c, ]$yLeftSupS  <- ThisP[ThisP$cFile == c-1, ]$yRightSupS
          ThisP[ThisP$cFile == c, ]$yLeftEdgeS <- ThisP[ThisP$cFile == c-1, ]$yRightEdgeS
          ThisP[ThisP$cFile == c, ]$yLeftInfS  <- outP[outP$idSector == paste("c", c,".l", layer[j]-1, sep = ""),]$yLeftSupS
        }
        rm(c)
        
        #### newLength ##
        ThisP$xLengthSupL <- sqrt((ThisP[ThisP$cFile == 1, ]$xLeftSupP - ThisP[ThisP$cFile == max(cFile), ]$xRightSupS)^2 + (ThisP[ThisP$cFile == 1, ]$yLeftSupS - ThisP[ThisP$cFile == max(cFile), ]$yRightSupS)^2)
        ThisP$xLengthInfL <- sqrt((ThisP[ThisP$cFile == 1, ]$xLeftInfP - ThisP[ThisP$cFile == max(cFile), ]$xRightInfS)^2 + (ThisP[ThisP$cFile == 1, ]$yLeftInfS - ThisP[ThisP$cFile == max(cFile), ]$yRightInfS)^2)
        
        ##### x ##
        ###### L.C1
        ThisP[ThisP$cFile == 1, ]$xLeftEdgeS   <- ThisP[ThisP$cFile == 1, ]$xLeftLed
        ThisP[ThisP$cFile == 1, ]$xRightSupS   <- (ThisP[ThisP$cFile == 1, ]$xLeftSupP) + (carCC[carCC$cFile == 1, ]$cumperLength * ThisP[ThisP$cFile == 1, ]$xLengthSupL / 100)
        
        #test on the x coordinates of the Left side in regard to x coordinates of the Right side
        if(all(ThisP[ThisP$cFile == 1, ]$xLeftLed >= ThisP[ThisP$cFile == 1, ]$xLeftSupP, ThisP[ThisP$cFile == 1, ]$xLeftLed >= ThisP[ThisP$cFile == 1, ]$xLeftInfP)){
          ThisP[ThisP$cFile == 1, ]$xLeftEdgeS <- min(c(ThisP[ThisP$cFile == 1, ]$xLeftSupP, ThisP[ThisP$cFile == 1, ]$xLeftInfP))
        }else{
          ThisP[ThisP$cFile == 1, ]$xLeftEdgeS <- ThisP[ThisP$cFile == 1, ]$xLeftLed
        }
        
        ###### L.C middle
        for(c in 2:(max(cFile)-1)){
          ThisP[ThisP$cFile == c, ]$xLeftSupP   <- ThisP[ThisP$cFile == c - 1, ]$xRightSupS
          ThisP[ThisP$cFile == c, ]$xLeftInfP   <- ThisP[ThisP$cFile == c - 1, ]$xRightInfS
          ThisP[ThisP$cFile == c, ]$xRightSupS  <- (ThisP[ThisP$cFile == 1, ]$xLeftSupP) + (carCC[carCC$cFile == c, ]$cumperLength * ThisP[ThisP$cFile == c, ]$xLengthSupL / 100)
          ThisP[ThisP$cFile == c, ]$xRightInfS  <- outP[outP$idSector == paste("c", c,".l", layer[j]-1, sep = ""), ]$xRightSupS
        }
        rm(c)
        
        ###### L.Cmax
        ThisP[ThisP$cFile == max(cFile), ]$xLeftSupP   <- ThisP[ThisP$cFile == (max(cFile) - 1), ]$xRightSupS
        ThisP[ThisP$cFile == max(cFile), ]$xLeftInfP   <- ThisP[ThisP$cFile == (max(cFile) - 1), ]$xRightInfS
        
        #check if the angle of the first row of produced cell is high: if yes: recalculation of the Xinter.min and Xinter.max
        newLength <- sqrt((ThisP[ThisP$cFile == max(cFile), ]$xRightSupS - ThisP[ThisP$cFile == 1, ]$xLeftSupP)^2 +
                            (ThisP[ThisP$cFile == max(cFile), ]$yRightSupS - ThisP[ThisP$cFile == 1, ]$yLeftSupS)^2)
        #the angle in regard to the X axis. To obtain the angle of correction, remove to 90 deg (which is the total angle)
        angle.deg.X <- acos((ThisP[ThisP$cFile == max(cFile), ]$xRightSupS - ThisP[ThisP$cFile == 1, ]$xLeftSupP) / newLength) * 180 / pi
        if(angle.deg.X > 5){
          ThisP[ThisP$cFile == 1, ]$xRightSupS <- (ThisP[ThisP$cFile == 1, ]$xLeftSupP) + ((carCC[1, ]$perLength * (ThisP[ThisP$cFile == max(cFile), ]$xRightSupS - ThisP[ThisP$cFile == 1, ]$xLeftSupP) / 100) * cos(angle.deg.X * pi/180))
          for(c in 2:(max(cFile) - 1)){
            ThisP[ThisP$cFile == c, ]$xLeftSupP  <- ThisP[ThisP$cFile == c - 1, ]$xRightSupS
            ThisP[ThisP$cFile == c, ]$xRightSupS <- (ThisP[ThisP$cFile == 1, ]$xLeftSupP) + ((carCC[c, ]$cumperLength * (ThisP[ThisP$cFile == max(cFile), ]$xRightSupS - ThisP[ThisP$cFile == 1, ]$xLeftSupP) / 100) * cos(angle.deg.X * pi/180))
          }
          rm(c)
          ThisP[ThisP$cFile == max(cFile), ]$xLeftSupP   <- ThisP[ThisP$cFile == (max(cFile) - 1), ]$xRightSupS
        } #end test on angle.deg
        
        
        ##### x Edges ##
        ThisP[ThisP$cFile == 1, ]$xRightEdgeS  <- (ThisP[ThisP$cFile == 1, ]$xRightSupS + ThisP[ThisP$cFile == 1, ]$xRightInfS) / 2
        for(c in 2:(max(cFile) - 1)){
          ThisP[ThisP$cFile == c, ]$xLeftEdgeS  <- ThisP[ThisP$cFile == c - 1, ]$xRightEdgeS
          ThisP[ThisP$cFile == c, ]$xRightEdgeS <- (ThisP[ThisP$cFile == c, ]$xRightSupS + ThisP[ThisP$cFile == c, ]$xRightInfS) / 2
        }
        rm(c)
        
        ThisP[ThisP$cFile == max(cFile), ]$xLeftEdgeS  <- ThisP[ThisP$cFile == (max(cFile) - 1), ]$xRightEdgeS
        #test on the x coordinates of the Right side in regard to x coordinates of the Left side
        if(all(ThisP[ThisP$cFile == max(cFile), ]$xRightLed <= ThisP[ThisP$cFile == max(cFile), ]$xRightSupS, ThisP[ThisP$cFile == max(cFile), ]$xRightLed < ThisP[ThisP$cFile == max(cFile), ]$RightInfP)){
          ThisP[ThisP$cFile == max(cFile), ]$xRightEdgeS <- mean(c(ThisP[ThisP$cFile == max(cFile), ]$xRightSupS, ThisP[ThisP$cFile == max(cFile), ]$xRightInfS))
        }else{
          ThisP[ThisP$cFile == max(cFile), ]$xRightEdgeS <- ThisP[ThisP$cFile == max(cFile), ]$xRightLed
        }
        
      } #end L.middle
    }
    
    outP <- rbind(outP, ThisP)
  }
  rm(j)  

  # 4. plot for checking #####
  if(make.plot == TRUE){
    title <- paste(as.character(unique(dataset$idD)), " sectorss", sep = "")
    
    plotP <- plot(Y.cor ~ X.cor, dataset, main =  title, cex.main = 0.8, type = "n",
                  xlim = c(-20, max(dataset$X.cor, na.rm = TRUE) + 5),
                  ylim = c(-25, max(outP$yRightSupS, na.rm = TRUE) * 1.01)
    )
    abline(v = seq(-10, max(dataset$X.cor, na.rm = TRUE) + 10, 10), col = "grey", lty = 3)
    abline(v = seq(0, max(dataset$X.cor, na.rm = TRUE) + 10, 50), col = "grey20", lty = 2)
    abline(h = seq(0, max(dataset$Y.cor, na.rm = TRUE), 10), col = "grey", lty = 3)
    abline(h = seq(0, max(dataset$Y.cor, na.rm = TRUE), 50), col = "grey20", lty = 2)
    
    #with sectors edges
    idSector <- unique(SecEdges$idSector)
    for(j in 1:length(idSector)){
      polygon(c(outP[outP$idSector == idSector[j],]$xLeftInfP, outP[outP$idSector == idSector[j],]$xLeftEdgeS, outP[outP$idSector == idSector[j],]$xLeftSupP,
                outP[outP$idSector == idSector[j],]$xRightSupS, outP[outP$idSector == idSector[j],]$xRightEdgeS, outP[outP$idSector == idSector[j],]$xRightInfS),
              c(outP[outP$idSector == idSector[j],]$yLeftInfS, outP[outP$idSector == idSector[j],]$yLeftEdgeS, outP[outP$idSector == idSector[j],]$yLeftSupS, 
                outP[outP$idSector == idSector[j],]$yRightSupS, outP[outP$idSector == idSector[j],]$yRightEdgeS, outP[outP$idSector == idSector[j],]$yRightInfS),
              border = "wheat2", lwd = 2)
    }
    rm(j)
    
    points(Y.cor ~ X.cor, dataset[dataset$type == "Vessel",], col = "darkturquoise", pch = 16)
    points(Y.cor ~ X.cor, dataset[dataset$type == "Fibers",], col = "brown", pch = 16)  
    points(Y.cor ~ X.cor, dataset[dataset$type == "Axial Parenchyma",], col = "purple2", pch = 15, cex = 0.8)
    points(Y.cor ~ X.cor, dataset[dataset$type == "APF" & dataset$phase == "E",], col = "maroon2", pch = 1)
    points(Y.cor ~ X.cor, dataset[dataset$type == "APF" & dataset$phase == "S",], col = "maroon2", pch = 16)
    points(Y.cor ~ X.cor, dataset[dataset$type == "Cambial",], col = "seagreen3", pch = 15)
    
    #layer label
    for(j in 1:length(layer)){
      text(-19, (outP[outP$layer == layer[j] & outP$cFile == 1,]$yLeftSupS + outP[outP$layer == layer[j] & outP$cFile == 1,]$yLeftInfS) / 2,
           layer[j], cex = 1, col = "grey60", font = 2)
      
    }
    rm(j)
    
    #cambial line label
    Line <- sort(unique(outP$cFile))
    for(j in 1:length(Line)){
      text((outP[outP$layer == 1 & outP$cFile == Line[j], ]$xLeftInfP + outP[outP$layer == 1 & outP$cFile == Line[j], ]$xRightInfS) / 2, -18,
           Line[j], cex = 1, col = "grey60", font = 2)
    }
    rm(j)
    
    legend(-22, -26, pch = c(15, 1, 16, 16, 15, 16),
           bty = "n", horiz = TRUE,
           y.intersp = 0.8, cex = 0.7,
           col = c("seagreen3", "maroon2", "maroon2", "brown", "purple2", "darkturquoise"),
           legend = c("Cc", "APF.e", "APF.l", "F", "AP", "V"))
    
    rm(plotP, title)
  }
  
  # 5. to export #####
  return(outP)
}

##### END OF THE FUNCTION sectors #####

#################### FUNCTION labelSec #####
# function to attribute to each cell a sectors id a calculation of sectors cell characteristic
# return list of 2 objects: coordinates of each cell and id sectors and sectors edges and number of cells

#' @dataset: dataframe with reoriented data (X.cor, Y.cor) and "type" column
#' @EdgeS: dataframe of the sectors edges

labelSec <- function(dataset, EdgeS, make.plot = TRUE){
  require(sp) || install.packages(sp)
  require(reshape2) || install.packages(reshape2)
  library(sp)
  library(reshape2)
  
  # 1. attribution of sectors id for each cell ####
  SecList  <- NULL ; secX <- NULL ; secY <- NULL ; newsecY <- NULL
  idSector <- unique(EdgeS$idSector)
  
  ## loop to create the list of sectors in the right format
  for(j in 1:length(idSector)){
    #extraction of the coordinates, the first and the last have to be the same to obtain a closed sectors
    secX <- c(EdgeS[EdgeS$idSector == idSector[j],]$xLeftInfP,  EdgeS[EdgeS$idSector == idSector[j],]$xLeftEdgeS,  EdgeS[EdgeS$idSector == idSector[j],]$xLeftSupP,
              EdgeS[EdgeS$idSector == idSector[j],]$xRightSupS, EdgeS[EdgeS$idSector == idSector[j],]$xRightEdgeS, EdgeS[EdgeS$idSector == idSector[j],]$xRightInfS,
              EdgeS[EdgeS$idSector == idSector[j],]$xLeftInfP)
    secY <- c(EdgeS[EdgeS$idSector == idSector[j],]$yLeftInfS,  EdgeS[EdgeS$idSector == idSector[j],]$yLeftEdgeS,  EdgeS[EdgeS$idSector == idSector[j],]$yLeftSupS,
              EdgeS[EdgeS$idSector == idSector[j],]$yRightSupS, EdgeS[EdgeS$idSector == idSector[j],]$yRightEdgeS, EdgeS[EdgeS$idSector == idSector[j],]$yRightInfS,
              EdgeS[EdgeS$idSector == idSector[j],]$yLeftInfS)
    #creation of a sector sector class with id
    newSec <- sp::Polygons(list(sp::Polygon(cbind(secX, secY))), ID = idSector[j])
    #save and compilation of all sector of the image on a list with each element is a sector
    SecList <- append(SecList, list(newSec))
  }
  rm(j)
  rm(secX, secY, newSec)
  
  ## creation of spatial sector
  secSpat <- sp::SpatialPolygons(SecList)
  
  ## creation of spatial points dataframe class
  row.names(dataset) <- dataset$idC
  cells <- sp::SpatialPointsDataFrame(coords = cbind(dataset$X.cor, dataset$Y.cor), data = dataset)
  
  ## attribution of sector id per cell
  testBis           <- data.frame(sp::over(cells, secSpat, returnList = FALSE))
  names(testBis)    <- "idSector.raw"
  testBis$idSector  <- idSector[testBis$idSector.raw]
  testBis$idC       <- row.names(testBis)
  
  cellPol <- testBis[, c("idC", "idSector")]
  
  ## addition if the sector id to the dataset dataframe
  outD <- merge(dataset, cellPol, by = "idC")
  outD$cFile <- gsub("\\D+", "", substr(outD$idSector, start = 1, stop = 3))
  outD$layer <- gsub("\\D+", "", substr(outD$idSector, start = nchar(outD$idSector) - 1, stop = nchar(outD$idSector)))
  
  rm(cellPol, cells, secSpat)
  rm(testBis)

  # 2. distance of APF cell to its sector edges ####
  idC   <- unique(outD$idC)
  ThisC <- NULL ; outDi <- NULL ; ThisP <- NULL
  
  outD$LeftInf   <- NA
  outD$LeftSup   <- NA
  outD$LeftEdge  <- NA
  outD$RightInf  <- NA
  outD$RightSup  <- NA
  outD$RightEdge <- NA
  outD$PolySide  <- NA

  for(d in 1:length(idC)){
    ThisC <- outD[outD$idC == idC[d], ]
    
    if(is.na(ThisC$idSector)){
      print(paste("No sector detected for image", ThisC$idD, "cell", idC[d], sep = " "))
    }else{
      ThisP <- EdgeS[which(EdgeS$idSector == ThisC$idSector), ]
      
      #changes of the x- and y- -Right- and -Left- -EdgeS values to be at the middle of the sector edge and not at the cell location
      ThisP$xRightEdgeS <- (ThisP$xRightInfS + ThisP$xRightSupS) / 2
      ThisP$yRightEdgeS <- (ThisP$yRightInfS + ThisP$yRightSupS) / 2
      ThisP$xLeftEdgeS  <- (ThisP$xLeftInfP + ThisP$xLeftSupP) / 2
      ThisP$yLeftEdgeS  <- (ThisP$yLeftInfS + ThisP$yLeftSupS) / 2
      
      ThisC$RightEdge   <- sqrt((ThisC$X.cor - ThisP$xRightEdgeS)^2 + (ThisC$Y.cor - ThisP$yRightEdgeS)^2)
      ThisC$LeftEdge    <- sqrt((ThisC$X.cor - ThisP$xLeftEdgeS)^2  + (ThisC$Y.cor - ThisP$yLeftEdgeS)^2)
      ThisC$LeftInf     <- sqrt((ThisC$X.cor - ThisP$xLeftInfP)^2   + (ThisC$Y.cor - ThisP$yLeftInfS)^2)
      ThisC$LeftSup     <- sqrt((ThisC$X.cor - ThisP$xLeftSupP)^2   + (ThisC$Y.cor - ThisP$yLeftSupS)^2)
      ThisC$RightInf    <- sqrt((ThisC$X.cor - ThisP$xRightInfS)^2  + (ThisC$Y.cor - ThisP$yRightInfS)^2)
      ThisC$RightSup    <- sqrt((ThisC$X.cor - ThisP$xRightSupS)^2  + (ThisC$Y.cor - ThisP$yRightSupS)^2)
      ThisC$PolySide    <- names(ThisC[, c("LeftSup", "LeftInf", "LeftEdge", "RightSup", "RightInf", "RightEdge")])[which.min(apply(ThisC[, c("LeftSup", "LeftInf", "LeftEdge", "RightSup", "RightInf", "RightEdge")], MARGIN = 2, min))]
      
      outDi <- rbind(ThisC, outDi)
    }
  }
  rm(d, ThisC, ThisP)
  
  outD <- outDi
  rm(outDi)
  
  # 3. type of cell per sector ####
  ## cell pooling
  outD$pool <- NA
  outD[outD$type %in% c("APF", "Fibers", "Axial Parenchyma"), ]$pool <- "APF"
  outD[outD$type == "Cambial", ]$pool <- "Cambial"
  if(nrow(outD[outD$type == "Vessel", ]) != 0){
    outD[outD$type == "Vessel", ]$pool  <- "Vessel"
  }
  
  ## presence or absence of cells within sector
  EdgeS$APF      <- NA ; EdgeS$vessel   <- NA
  EdgeS$idVessel <- NA ; EdgeS$cell     <- NA
  
    #loop to attribute the type of cell per sector
  typeC <-  NULL ; maxV <- NULL
  for(j in 1:length(idSector)){
    typeC <- unique(outD[outD$idSector == idSector[j],]$pool)
    
    if(length(typeC) != 0){
      # APF
      if(any(typeC == "APF")){ EdgeS[EdgeS$idSector == idSector[j],]$APF <- 1 }
      else{ EdgeS[EdgeS$idSector == idSector[j],]$APF <- 0 }
      # Vessel
      if(any(typeC == "Vessel")){
        EdgeS[EdgeS$idSector == idSector[j],]$vessel <- 1
        if(length(outD[outD$idSector == idSector[j] & outD$type == "Vessel",]$idC) > 1){ # if there is several vessel within the sector, to choose the widest
          maxV <- outD[outD$idSector == idSector[j] & outD$type == "Vessel",]
          EdgeS[EdgeS$idSector == idSector[j],]$idVessel <- maxV[which.max(maxV$Feret.v),]$idC
        }else{
          EdgeS[EdgeS$idSector == idSector[j],]$idVessel <- outD[outD$idSector == idSector[j] & outD$type == "Vessel",]$idC
        }
        }else{ EdgeS[EdgeS$idSector == idSector[j],]$vessel <- 0 }
      # Cell
      EdgeS[EdgeS$idSector == idSector[j],]$cell <- 1
      
    }else{
      EdgeS[EdgeS$idSector == idSector[j],]$APF    <- 0
      EdgeS[EdgeS$idSector == idSector[j],]$vessel <- 0
      EdgeS[EdgeS$idSector == idSector[j],]$cell   <- 0
    }
  }
  rm(j, typeC, maxV)
  
  # 4. number of APF and cambial cells in each sector #####
  cDens    <- stats::aggregate(data.frame(nbcell = outD[outD$pool == "APF",]$idSector), list(idSector = outD[outD$pool == "APF",]$idSector), length)
  template <- data.frame(idSector = EdgeS$idSector, cFile = EdgeS$cFile, layer = EdgeS$layer)
  densI    <- merge(template, cDens, by = "idSector", all.x = TRUE)
  if(any(is.na(densI$nbcell))){
    densI[is.na(densI$nbcell), ]$nbcell <- 0
  }
  rm(cDens, template)
  
  secD <- merge(EdgeS, densI, by = c("idSector", "layer", "cFile"), all.x = TRUE)

  # 5. plot for checking #####
  if(make.plot == TRUE){
    title <- paste(as.character(unique(dataset$idD)), " labelled sectors", sep = "")
    
    plotL <- plot(Y.cor ~ X.cor, dataset, main =  title, cex.main = 0.8, type = "n",
                  xlim = c(-10, max(EdgeS$xRightSupS, EdgeS$xRightInfS, na.rm = TRUE) + 5),
                  ylim = c(-10, max(EdgeS$ySupL, na.rm = TRUE) * 1.01))
    abline(v = seq(-10, max(EdgeS$xRightSupS, EdgeS$xRightInfS, na.rm = TRUE) + 10, 10), col = "grey", lty = 3)
    abline(h = seq(0, max(dataset$Y.cor, na.rm = TRUE), 50), col = "grey", lty = 3)
    
    for(j in 1:length(idSector)){
      if(EdgeS[EdgeS$idSector == idSector[j],]$APF == 1){
        polygon(c(EdgeS[EdgeS$idSector == idSector[j],]$xLeftInfP,  EdgeS[EdgeS$idSector == idSector[j],]$xLeftEdgeS,  EdgeS[EdgeS$idSector == idSector[j],]$xLeftSupP,
                  EdgeS[EdgeS$idSector == idSector[j],]$xRightSupS, EdgeS[EdgeS$idSector == idSector[j],]$xRightEdgeS, EdgeS[EdgeS$idSector == idSector[j],]$xRightInfS),
                c(EdgeS[EdgeS$idSector == idSector[j],]$yLeftInfS,  EdgeS[EdgeS$idSector == idSector[j],]$yLeftEdgeS,  EdgeS[EdgeS$idSector == idSector[j],]$yLeftSupS, 
                  EdgeS[EdgeS$idSector == idSector[j],]$yRightSupS, EdgeS[EdgeS$idSector == idSector[j],]$yRightEdgeS, EdgeS[EdgeS$idSector == idSector[j],]$yRightInfS),
                col = adjustcolor(col = "maroon2", alpha.f = 0.1))
      }
      
      if(EdgeS[EdgeS$idSector == idSector[j],]$vessel == 1){
        polygon(c(EdgeS[EdgeS$idSector == idSector[j],]$xLeftInfP,  EdgeS[EdgeS$idSector == idSector[j],]$xLeftEdgeS,  EdgeS[EdgeS$idSector == idSector[j],]$xLeftSupP,
                  EdgeS[EdgeS$idSector == idSector[j],]$xRightSupS, EdgeS[EdgeS$idSector == idSector[j],]$xRightEdgeS, EdgeS[EdgeS$idSector == idSector[j],]$xRightInfS),
                c(EdgeS[EdgeS$idSector == idSector[j],]$yLeftInfS,  EdgeS[EdgeS$idSector == idSector[j],]$yLeftEdgeS,  EdgeS[EdgeS$idSector == idSector[j],]$yLeftSupS, 
                  EdgeS[EdgeS$idSector == idSector[j],]$yRightSupS, EdgeS[EdgeS$idSector == idSector[j],]$yRightEdgeS, EdgeS[EdgeS$idSector == idSector[j],]$yRightInfS),
                col = adjustcolor(col = "darkturquoise", alpha.f = 0.1))
      }
      
      polygon(c(EdgeS[EdgeS$idSector == idSector[j],]$xLeftInfP,  EdgeS[EdgeS$idSector == idSector[j],]$xLeftEdgeS,  EdgeS[EdgeS$idSector == idSector[j],]$xLeftSupP,
                EdgeS[EdgeS$idSector == idSector[j],]$xRightSupS, EdgeS[EdgeS$idSector == idSector[j],]$xRightEdgeS, EdgeS[EdgeS$idSector == idSector[j],]$xRightInfS),
              c(EdgeS[EdgeS$idSector == idSector[j],]$yLeftInfS,  EdgeS[EdgeS$idSector == idSector[j],]$yLeftEdgeS,  EdgeS[EdgeS$idSector == idSector[j],]$yLeftSupS, 
                EdgeS[EdgeS$idSector == idSector[j],]$yRightSupS, EdgeS[EdgeS$idSector == idSector[j],]$yRightEdgeS, EdgeS[EdgeS$idSector == idSector[j],]$yRightInfS),
              border = "wheat2", lwd = 2)
    }
    rm(j)
    
    points(Y.cor ~ X.cor, dataset[dataset$type == "Vessel",], col = "darkturquoise", pch = 16)
    points(Y.cor ~ X.cor, dataset[dataset$type == "Fibers",], col = "brown", pch = 16)  
    points(Y.cor ~ X.cor, dataset[dataset$type == "Axial Parenchyma",], col = "purple2", pch = 15, cex = 0.8)
    points(Y.cor ~ X.cor, dataset[dataset$type == "APF" & dataset$phase == "E",], col = "maroon2", pch = 1)
    points(Y.cor ~ X.cor, dataset[dataset$type == "APF" & dataset$phase == "S",], col = "maroon2", pch = 16)
    points(Y.cor ~ X.cor, dataset[dataset$type == "Cambial",], col = "seagreen3", pch = 15)
    
    layer <- unique(EdgeS$layer)
    for(j in 1:length(layer)){
      text(-11, (EdgeS[EdgeS$layer == layer[j] & EdgeS$cFile == 1,]$yLeftSupS + EdgeS[EdgeS$layer == layer[j] & EdgeS$cFile == 1,]$yLeftInfS) / 2,
           layer[j], cex = 1)
      
    }
    rm(j, plotL, title)
  }

  # 6. to export ####
  outA <- list(outD, secD)
  return(outA)
  
}

##### END OF THE FUNCTION labelSec #####

#################### FUNCTION vNeigh #####
# function to attribute to each cell a sector id a calculation of sector cell characteristic
# return list of 2 objects: Vessel spatial characteristics and Vessel neighbor cells

#' @dataset: dataframe with reoriented data (X.cor, Y.cor) and "type" column

vNeigh <- function(dataset){
  require(deldir) || install.packages(deldir)
  library(deldir)
  
  # 1. data and delaunay triangulation ####
  ThisE <- dataset[, c("idC", "Y.cor", "X.cor", "idSector", "cFile", "layer", "type", "idD")]
  ThisE <- ThisE[order(ThisE$idC),]
  DMat  <- deldir::deldir(ThisE$X.cor, ThisE$Y.cor)
  
  ## calculation of the norm
  norme <- function(x1, x2, y1, y2) {sqrt((x1 - x2)^2 + (y1 - y2)^2)}
  
  DMat$delsgs$dist <- with(DMat$delsgs, norme(DMat$delsgs$x1, DMat$delsgs$x2, DMat$delsgs$y1, DMat$delsgs$y2))
  
  ## create empty adjacency matrix, fill in distances
  adj <- matrix(NA, nrow = nrow(ThisE), ncol = nrow(ThisE))
  adj[as.matrix(DMat$delsgs[c("ind1", "ind2")])] <- DMat$delsgs$dist
  
  if(nrow(ThisE[ThisE$type == "Vessel", ]) != 0){
    # 2. Vessel characteristics and distances ####
    ## vessel characteristics
    ThoseVessel           <- ThisE[ThisE$type == "Vessel", c("idC", "idSector", "cFile", "layer")]
    names(ThoseVessel)[2] <- "idSectorV"
    ThoseVessel           <- ThoseVessel[order(ThoseVessel$idC),]
    ThoseVessel$nbNeigh   <- DMat$summary[ThoseVessel$idC,]$n.tri
    
    ## vessel neighbors id, idSector of neighbor and distance between vessel and neighbor
    outN      <- NULL
    ThoseCell <- NULL
    
    for(j in 1:nrow(ThoseVessel)){
      ThoseCell               <- data.frame(idV = ThoseVessel[j, ]$idC, idCNeigh = c(1:nrow(ThisE)), distance.From = adj[ThoseVessel[j, ]$idC,])
      ThoseCell$distance.To   <- adj[ ,ThoseVessel[j, ]$idC]
      ThoseCell$distance      <- NA
      ThoseCell[which(!is.na(ThoseCell$distance.To)),]$distance   <- ThoseCell[which(!is.na(ThoseCell$distance.To)),]$distance.To
      ThoseCell[which(!is.na(ThoseCell$distance.From)),]$distance <- ThoseCell[which(!is.na(ThoseCell$distance.From)),]$distance.From
      ThoseCell$distance.From <- NULL
      ThoseCell$distance.To   <- NULL
      
      ThoseCell <- na.omit(ThoseCell)
      outN      <- rbind(ThoseCell, outN)
    }
    rm(j, ThoseCell)
    
    outN$idSectorCN <- NA
    for(p in 1:nrow(outN)){ outN[p,]$idSectorCN <- dataset[dataset$idC == outN[p,]$idCNeigh,]$idSector }
    rm(p)
    
    # 3. to export ####
    ThoseVessel$idD <- unique(ThisE$idD)
    outN$idD        <- unique(ThisE$idD)
    
    ThoseVessel <- ThoseVessel[, c("idD", "idC", "idSectorV", "cFile", "layer", "nbNeigh")]
    outN        <- outN[, c("idD", "idV", "idCNeigh", "distance", "idSectorCN")]
  }else{
    ThoseVessel <- data.frame(idD = unique(ThisE$idD), idC = NA, idSectorV = NA, cFile = NA, layer = NA, nbNeigh = NA)
    outN        <- data.frame(idD = unique(ThisE$idD), idV = NA, idCNeigh = NA, distance = NA, idSectorCN = NA)
  }
  outV <- list(ThoseVessel, outN)
  return(outV)
}

##### END OF THE FUNCTION vNeigh #####

#################### FUNCTION vesselS #####
# function to attribute to each cell a sector id a calculation of sector cell characteristic
# return 1 object: outVtot (coordinates of each vessel and id sector)

#' @dataset: dataframe with reoriented data (X.cor, Y.cor) and "type" column
#' @EdgeS: dataframe of the sector edges
#' @CarV: dataframe with vessel characteristic (id image, idC of the vessel, idSector of the origin sector of the vessel, cFile and layer of idSector, number of adjacent cells)
#' @NeighV: dataframe with idsector of the neighbor cell, id image, idC of the Vessel, idC of the neighbor cell, distance cell-vessel

vesselS <- function(dataset, EdgeS, CarV, NeighV, make.plot = TRUE){
  if(all(!is.na(CarV$idSectorV))){
    # vector and dataframe
    idSector   <- unique(EdgeS$idSector)
    outD       <- dataset
    outD$layer <- as.numeric(outD$layer)
    outD$cFile <- as.numeric(outD$cFile)
    
    CarV$layer <- as.numeric(gsub("\\D+", "", substr(CarV$idSectorV, start = nchar(CarV$idSectorV) - 2, stop = nchar(CarV$idSectorV))))
    
    # 1. loop to manage vessel (vessel per vessel): test adjacent sectors, if empty or 1 APF only = put 1 in vessel column of EdgeS #####
    ## dataframe for disqualifiers criteria
    sec      <- c("secU", "secUR", "secR", "secDR", "secD", "secDL", "secL", "secUL")
    location <- c("LeftInf",  "RightInf",                                                      # secU
                  "LeftInf",  "LeftEdge",                                                      # secUR
                  "LeftSup",  "LeftEdge",  "LeftInf",                                          # secR
                  "LeftSup",  "LeftEdge",                                                      # secDR
                  "LeftSup",  "RightSup",                                                      # secD
                  "RightSup", "RightEdge",                                                     # secDL
                  "RightSup", "RightEdge", "RightInf",                                         # secL
                  "RightInf", "RightEdge")                                                     # secUL
    
    direction <- c("up", "up right", "right", "down right", " down", "down left", "left", "up left")
    qualif    <- data.frame(sec = rep(sec, c(2,2,3,2,2,2,3,2)), location = location, direction = rep(direction, c(2,2,3,2,2,2,3,2)))
    rm(direction, location, sec)
    
    ## empty objects
    ThisN   <- NULL ; ThisNr  <- NULL ; ThisP  <- NULL ; ThisC  <- NULL ; ThisL  <- NULL ; ThisVc <- NULL ; secRef <- NULL
    secAll  <- NULL ; Ctest   <- NULL ; CtestR <- NULL ; CtestL <- NULL ; Ltest  <- NULL ; Ptest  <- NULL ; Rtest   <- NULL ; secV <- NULL
    up      <- NULL ; down    <- NULL ; right  <- NULL ; left   <- NULL
    cLayerU <- NULL ; cLayerD <- NULL ; test   <- NULL ; testL  <- NULL ; testLu <- NULL ; testLd <- NULL
    outV    <- NULL ; outVtot <- NULL
    pPopo   <- NULL ; title   <- NULL ; layer  <- NULL ; Line   <- NULL ; secV  <- NULL
    
    ## loop j : j = vessel id
    for(j in 1:nrow(CarV)){
      ### 1.1.extraction of the vessel neighbors group
      ThisN       <- NeighV[NeighV$idV == CarV[j, ]$idC, ]
      ThisN$cFile <- as.numeric(gsub("\\D+", "", substr(ThisN$idSectorCN, start = 2, stop = 4)))
      ThisN$layer <- as.numeric(gsub("\\D+", "", substr(ThisN$idSectorCN, start = nchar(ThisN$idSectorCN) - 2, stop = nchar(ThisN$idSectorCN))))
      
      ### 1.2. creation of the dataframe of sectors between ThisN and vessel center
      ThisN$cFile <- as.numeric(ThisN$cFile)
      ThisN$layer <- as.numeric(ThisN$layer)
      #extraction of sector origin cFile and layer
      ThisC   <- as.numeric(CarV[j, ]$cFile)
      ThisL   <- as.numeric(CarV[j, ]$layer)
      secRef <- CarV[j, ]$idSectorV
      #creation of ThisP dataframe and addition of the column pool for each idCNeigh to ThisN
      ThisP           <- outD[, c("idSector", "PolySide", "pool", "idC")]
      names(ThisP)[4] <- "idCNeigh"
      ThisN           <- base::merge(ThisN, ThisP[, c("idCNeigh", "pool")], by = "idCNeigh")
      ThisP           <- NULL
      
      ThisN$cLdif  <- ThisC - ThisN$cFile
      ThisN$laydif <- ThisL - ThisN$layer
      
      if(any(ThisN$idSectorCN == secRef)){
        secAll <- ThisN[ThisN$idSectorCN == secRef, ]
        #dataframe with position of cell within sector
        Ptest           <- outD[which(outD$idSector %in% secAll$idSectorCN), c("idSector", "PolySide", "idC")]
        names(Ptest)[1] <- "idSectorCN"
        Ptest           <- base::merge(secAll, Ptest[Ptest$idC != unique(ThisN$idV), ], by = "idSectorCN", all.x = TRUE)
        #dataframe with validated sector: the ones which are not validated will be deleted with test
        secV          <- data.frame(idSectorCN = NA,
                                     cFile = c(rep(ThisC, 3), ThisC + 1, ThisC - 1),
                                     layer = c(ThisL, ThisL + 1, ThisL - 1, rep(ThisL, 2)),
                                     side = c("origin", "secU", "secD", "secR", "secL"),
                                     tested = 0, validated = 0)
        secV$idSectorCN <- paste("c", secV$cFile, ".l", secV$layer, sep = "")
        secV[1, ]$tested    <- 1
        secV[1, ]$validated <- 1
        
        #test 1: position of cell within sector to know which adjacent sector to test
        if(length(grep("Inf",   Ptest$PolySide)) == 0){ secV[secV$side == "secD", ]$tested <- 1 }
        if(length(grep("Sup",   Ptest$PolySide)) == 0){ secV[secV$side == "secU", ]$tested <- 1 }
        if(length(grep("Left",  Ptest$PolySide)) == 0){ secV[secV$side == "secL", ]$tested <- 1 }
        if(length(grep("Right", Ptest$PolySide)) == 0){ secV[secV$side == "secR", ]$tested <- 1 }
        
        #keeping only the tested sector
        secV <- secV[secV$tested == 1, ]
        
        #dataframe with tested sector
        ThisP           <- outD[which(outD$idSector %in% secV[secV$side != "origin", ]$idSectorCN), c("idSector", "PolySide", "idC")]
        names(ThisP)[1] <- "idSectorCN"
        ThisP           <- base::merge(secV, ThisP, by = "idSectorCN")
        
        #empty tested sectors are validated
        if(nrow(ThisP) != 0){
          secV[secV$tested == 1 & !(secV$idSectorCN %in% ThisP$idSectorCN), ]$validated <- 1
          
          #test 2 on cell position on adjacent sectors
          for(p in 1:length(unique(ThisP$idSectorCN))){
            if(!any(any(ThisP[ThisP$idSectorCN == unique(ThisP$idSectorCN)[p],]$PolySide %in% qualif[qualif$sec %in% unique(ThisP[ThisP$idSectorCN == unique(ThisP$idSectorCN)[p],]$side), ]$location))){
              secV[secV$idSectorCN == unique(ThisP$idSectorCN)[p], ]$validated <- 1
            }
          }
          rm(p)
        }else{
          secV$validated <- 1
        }
        
        #keeping only the validated sector
        secV <- secV[secV$validated == 1, ]
        
        
      }else{
        #### case of adjacent vessel
        Cref <- NULL ; Lref <- NULL
        if(any(ThisN$pool == "Vessel")){
          ThisVc <- ThisN[which(ThisN$pool == "Vessel"), ]
          for(v in 1:nrow(ThisVc)){
            if(abs(ThisVc[v, ]$cLdif) > abs(ThisVc[v, ]$laydif)){               #adj vessel on the left or the right side
              if(abs(ThisVc[v, ]$cLdif) == 1){
                Cref  <- rbind(Cref, data.frame(cFile = ThisC, layer = ThisVc[v, ]$layer))
              }else{
                Cref  <- rbind(Cref, data.frame(cFile = round(mean(c(ThisVc[v, ]$cFile, ThisC)), 0), layer = ThisVc[v, ]$layer))
              }
              
            }else if(abs(ThisVc[v, ]$cLdif) < abs(ThisVc[v, ]$laydif)){         #adj vessel on the top or the bottom side
              if(abs(ThisVc[v, ]$laydif) == 1){
                Lref  <- rbind(Lref, data.frame(cFile = ThisVc[v, ]$cFile, layer = ThisL))
              }else{
                Lref  <- rbind(Lref, data.frame(cFile = ThisVc[v, ]$cFile, layer = round(mean(c(ThisVc[v, ]$layer, ThisL)), 0)))
              }
              
            }else if(abs(ThisVc[v, ]$cLdif) == abs(ThisVc[v, ]$laydif)){        #adj vessel with same cFile and layer differences
              if(all(ThisVc[v, ]$cLdif == 0, ThisVc[v, ]$laydif == 0)){         #case of 2 vessels within the same sector
                Cref  <- rbind(Cref, data.frame(cFile = ThisC, layer = ThisVc[v, ]$layer))
                Lref  <- rbind(Lref, data.frame(cFile = ThisVc[v, ]$cFile, layer = ThisL))
              }else{                                                            #case of vessel in different sector
                Cref  <- rbind(Cref, data.frame(cFile = round(mean(c(ThisVc[v, ]$cFile, ThisC)), 0), layer = ThisVc[v, ]$layer))
                Lref  <- rbind(Lref, data.frame(cFile = ThisVc[v, ]$cFile, layer = round(mean(c(ThisVc[v, ]$layer, ThisL)), 0)))
              }
            }
          }
          rm(v, ThisVc)
        }
        
        #extraction of Layer to test
        if(ThisL == 1){ testLd <- ThisL                           }else{ testLd <- ThisL - 1}
        if(ThisL == max(as.numeric(outD$layer))){ testLu <- ThisL }else{ testLu <- ThisL + 1}
        
        if(is.null(Lref)){
          if(max(ThisN$layer) == ThisL){ up   <- testLu }else if(max(ThisN$layer) < (ThisL + 1)){ up <- ThisL + 1 }else{ up   <- seq(testLu, max(ThisN$layer), 1) }
          
          if(all(!(ThisN[which.min(ThisN$layer), ]$cFile %in% c((ThisC - 1):(ThisC + 1))), min(ThisN$layer) > 1)){
            if(min(ThisN$layer) == ThisL){ down <- testLd - 1 }else{ down <- sort(seq(min(ThisN$layer), testLd, 1), decreasing = TRUE) }
          }else{
            if(min(ThisN$layer) == ThisL){ down <- testLd }else{ down <- sort(seq(min(ThisN$layer), testLd, 1), decreasing = TRUE) }
          } 
          
        }else{
          if(max(Lref$layer) == ThisL){
            if(all(max(ThisN[which(ThisN$pool == "Vessel"), ]$layer) < ThisL,
                   max(ThisN$layer) >= testLu)){   up <- seq(ThisL, max(ThisN$layer), 1)
            }else{                                 up <- ThisL }
            
          }else if(max(Lref$layer) > ThisL){       up <- seq(testLu, max(Lref$layer), 1)
          }else{                                   up <- seq(ThisL, max(ThisN$layer), 1) }
          
          if(min(Lref$layer) == ThisL){
            if(all(min(ThisN[which(ThisN$pool == "Vessel"), ]$layer) > ThisL, 
               min(ThisN$layer) <= testLd)){                               down <- sort(seq(min(ThisN$layer), testLd, 1), decreasing = TRUE)
            }else{                                                         down <- ThisL }
          }else if(min(Lref$layer) < ThisL){                               down <- sort(seq(min(Lref$layer),  testLd, 1), decreasing = TRUE)
          }else{                                                           down <- sort(seq(min(ThisN$layer), ThisL, 1), decreasing = TRUE) }
        }

        #max and min cFile per layer
        layer <- sort(unique(c(ThisL, up, down)))
        Rtest <- data.frame(layer = seq(min(layer), max(layer), 1), maxC = NA, minC = NA)
        layer <- NULL
        #importation of Cref data
        if(!is.null(Cref)){
          for(l in 1:nrow(Cref)){
            if(all(Cref[l, ]$cFile < ThisC, Cref[l, ]$layer %in% Rtest$layer,               
                   any(is.na(Rtest[Rtest$layer ==  Cref[l, ]$layer, ]$minC),              
                       Rtest[Rtest$layer ==  Cref[l, ]$layer, ]$minC < Cref[l, ]$cFile))){         # Cref$cFile at the left side + Cref$layer included in Rtest + empty Rtest$layer OR ref$cFile is inferior to the Rtest$layer one
              Rtest[Rtest$layer == Cref[l, ]$layer, ]$minC <- Cref[l, ]$cFile
            }else if(all(Cref[l, ]$cFile < ThisC, Cref[l, ]$layer > max(Rtest$layer),
                         any(is.na(Rtest[Rtest$layer ==  max(Rtest$layer), ]$minC),
                             Rtest[Rtest$layer ==  max(Rtest$layer), ]$minC < Cref[l, ]$cFile))){  # Cref$cFile at the left side + Cref$layer superior to Rtest$layer + empty Rtest$layer OR ref$cFile is inferior to the Rtest$layer one
              Rtest[Rtest$layer ==  max(Rtest$layer), ]$minC <- Cref[l, ]$cFile
            }else if(all(Cref[l, ]$cFile < ThisC, Cref[l, ]$layer < min(Rtest$layer),
                         any(is.na(Rtest[Rtest$layer ==  min(Rtest$layer), ]$minC),
                             Rtest[Rtest$layer == min(Rtest$layer), ]$minC < Cref[l, ]$cFile))){   # Cref$cFile at the left side + Cref$layer inferior to Rtest$layer + empty Rtest$layer OR ref$cFile is inferior to the Rtest$layer one
              Rtest[Rtest$layer ==  min(Rtest$layer), ]$minC <- Cref[l, ]$cFile
            }
            
            if(all(Cref[l, ]$cFile > ThisC, Cref[l, ]$layer %in% Rtest$layer,               
                   any(is.na(Rtest[Rtest$layer ==  Cref[l, ]$layer, ]$maxC),              
                       Rtest[Rtest$layer ==  Cref[l, ]$layer, ]$maxC > Cref[l, ]$cFile))){         # Cref$cFile at the right side + Cref$layer included in Rtest + empty Rtest$layer OR ref$cFile is inferior to the Rtest$layer one
              Rtest[Rtest$layer == Cref[l, ]$layer, ]$maxC <- Cref[l, ]$cFile
            }else if(all(Cref[l, ]$cFile > ThisC, Cref[l, ]$layer > max(Rtest$layer),
                         any(is.na(Rtest[Rtest$layer ==  max(Rtest$layer), ]$maxC),
                             Rtest[Rtest$layer ==  max(Rtest$layer), ]$maxC > Cref[l, ]$cFile))){  # Cref$cFile at the right side + Cref$layer superior to Rtest$layer + empty Rtest$layer OR ref$cFile is inferior to the Rtest$layer one
              Rtest[Rtest$layer ==  max(Rtest$layer), ]$maxC <- Cref[l, ]$cFile
            }else if(all(Cref[l, ]$cFile > ThisC, Cref[l, ]$layer < min(Rtest$layer),
                         any(is.na(Rtest[Rtest$layer ==  min(Rtest$layer), ]$maxC),
                             Rtest[Rtest$layer == min(Rtest$layer), ]$maxC > Cref[l, ]$cFile))){   # Cref$cFile at the right side + Cref$layer inferior to Rtest$layer + empty Rtest$layer OR ref$cFile is inferior to the Rtest$layer one
              Rtest[Rtest$layer ==  min(Rtest$layer), ]$maxC <- Cref[l, ]$cFile
            }
          }
          rm(l)
        }
        #minC ThisL
        if(is.na(Rtest[Rtest$layer == ThisL, ]$minC)){
          if(nrow(ThisN[ThisN$layer == ThisL & ThisN$cFile <= ThisC, ]) != 0){
            Rtest[Rtest$layer == ThisL, ]$minC <- max(ThisN[ThisN$layer == ThisL & ThisN$cFile <= ThisC, ]$cFile)
          }else if(all(nrow(ThisN[ThisN$cFile <= ThisC, ]) != 0,
                       nrow(ThisN[ThisN$layer == (ThisL - 1) & ThisN$cFile <= ThisC, ]) == 0,
                       nrow(ThisN[ThisN$layer == (ThisL + 1) & ThisN$cFile <= ThisC, ]) == 0)){
            Rtest[Rtest$layer == ThisL, ]$minC <- min(ThisN[ThisN$cFile <= ThisC, ]$cFile)
          }else if(all(nrow(ThisN[ThisN$layer == (ThisL - 1) & ThisN$cFile <= ThisC, ]) == 0,
                       nrow(ThisN[ThisN$layer == (ThisL + 1) & ThisN$cFile <= ThisC, ]) == 0)){
            Rtest[Rtest$layer == ThisL, ]$minC <- 1
          }else{
            Rtest[Rtest$layer == ThisL, ]$minC <- min(c(ThisN[ThisN$layer == (ThisL - 1) & ThisN$cFile <= ThisC, ]$cFile,
                                                        ThisN[ThisN$layer == (ThisL + 1) & ThisN$cFile <= ThisC, ]$cFile))
          }
        }
        #maxC ThisL
        if(is.na(Rtest[Rtest$layer == ThisL, ]$maxC)){
          if(nrow(ThisN[ThisN$layer == ThisL & ThisN$cFile > ThisC, ]) != 0){
            Rtest[Rtest$layer == ThisL, ]$maxC <- min(ThisN[ThisN$layer == ThisL & ThisN$cFile > ThisC, ]$cFile)
          }else if(all(nrow(ThisN[ThisN$cFile > ThisC, ]) != 0,
                       nrow(ThisN[ThisN$layer == (ThisL - 1) & ThisN$cFile > ThisC, ]) == 0,
                       nrow(ThisN[ThisN$layer == (ThisL + 1) & ThisN$cFile > ThisC, ]) == 0)){
            Rtest[Rtest$layer == ThisL, ]$maxC <- max(ThisN[ThisN$cFile > ThisC, ]$cFile)
          }else if(all(nrow(ThisN[ThisN$layer == (ThisL - 1) & ThisN$cFile > ThisC, ]) == 0,
                       nrow(ThisN[ThisN$layer == (ThisL + 1) & ThisN$cFile > ThisC, ]) == 0)){
            Rtest[Rtest$layer == ThisL, ]$maxC <- max(outD$cFile)
          }else{
            Rtest[Rtest$layer == ThisL, ]$maxC <- max(c(ThisN[ThisN$layer == (ThisL - 1) & ThisN$cFile > ThisC, ]$cFile, ThisN[ThisN$layer == (ThisL + 1) & ThisN$cFile > ThisC, ]$cFile))
          }
        }
        
        #another layers
        if(nrow(Rtest) > 1){
          if(length(seq(1:nrow(Rtest))) > 2){
            Ltest <- NULL
            for(l in seq(2, nrow(Rtest)-1, 1)){
              Ltest <- Rtest[l, ]$layer
              #minC
              if(all(is.na(Rtest[l, ]$minC), nrow(ThisN[ThisN$layer == Ltest & ThisN$cFile < ThisC, ]) != 0)){
                Rtest[Rtest$layer == Ltest, ]$minC <- max(ThisN[ThisN$layer == Ltest & ThisN$cFile < ThisC, ]$cFile)
              }else if(all(is.na(Rtest[l, ]$minC), nrow(ThisN[ThisN$layer == Ltest & ThisN$cFile < ThisC, ]) == 0)){
                if(all(nrow(ThisN[ThisN$layer == (Ltest - 1) & ThisN$cFile < ThisC, ]) == 0,
                       nrow(ThisN[ThisN$layer == (Ltest + 1) & ThisN$cFile < ThisC, ]) == 0,
                       l != 2)){
                  Rtest[Rtest$layer == Ltest, ]$minC <- Rtest[Rtest$layer == (Ltest - 1), ]$minC
                }else{
                  Rtest[Rtest$layer == Ltest, ]$minC <- Rtest[Rtest$layer == ThisL, ]$minC
                }
                
              }
              #maxC
              if(all(is.na(Rtest[l, ]$maxC), nrow(ThisN[ThisN$layer == Ltest & ThisN$cFile > ThisC, ]) != 0)){
                Rtest[Rtest$layer == Ltest, ]$maxC <- min(ThisN[ThisN$layer == Ltest & ThisN$cFile > ThisC, ]$cFile)
              }else if(all(is.na(Rtest[l, ]$maxC), nrow(ThisN[ThisN$layer == Ltest & ThisN$cFile > ThisC, ]) == 0)){
                if(all(nrow(ThisN[ThisN$layer == (Ltest - 1) & ThisN$cFile > ThisC, ]) == 0,
                       nrow(ThisN[ThisN$layer == (Ltest + 1) & ThisN$cFile > ThisC, ]) == 0,
                       l != 2)){
                  Rtest[Rtest$layer == Ltest, ]$maxC <- Rtest[Rtest$layer == (Ltest - 1), ]$maxC
                }else{
                  Rtest[Rtest$layer == Ltest, ]$maxC <- Rtest[Rtest$layer == ThisL, ]$maxC
                }
              }
            }
            rm(l)
          }
          
          #first layer
          if(all(is.na(Rtest[1, ]$minC), Rtest[1, ]$layer != ThisL)){
            Rtest[1, ]$minC <- min(ThisN[ThisN$layer == Rtest[1, ]$layer & ThisN$cFile < ThisC, ]$cFile, Rtest[2, ]$minC) }
          if(all(is.na(Rtest[1, ]$maxC), Rtest[1, ]$layer != ThisL)){
            Rtest[1, ]$maxC <- max(ThisN[ThisN$layer == Rtest[1, ]$layer & ThisN$cFile > ThisC, ]$cFile,  Rtest[2, ]$maxC) }
          #last layer
          if(all(is.na(Rtest[nrow(Rtest), ]$minC), Rtest[nrow(Rtest), ]$layer != ThisL)){
            Rtest[nrow(Rtest), ]$minC <- min(ThisN[ThisN$layer == Rtest[nrow(Rtest), ]$layer & ThisN$cFile < ThisC, ]$cFile, Rtest[(nrow(Rtest) - 1), ]$minC) }
          if(all(is.na(Rtest[nrow(Rtest), ]$maxC), Rtest[nrow(Rtest), ]$layer != ThisL)){
            Rtest[nrow(Rtest), ]$maxC <- max(ThisN[ThisN$layer == Rtest[nrow(Rtest), ]$layer & ThisN$cFile > ThisC, ]$cFile, Rtest[(nrow(Rtest) - 1), ]$maxC) }
        }
        
        #number of radial file per layer
        Rtest$interC <- Rtest$maxC - Rtest$minC + 1
        
        #### creation the dataframe with all sectors to test ##
        layer <- NULL
        Line  <- NULL
        for(v in 1:nrow(Rtest)){
          layer <- c(layer, rep(Rtest[v, ]$layer, Rtest[v, ]$interC))
          Line  <- c(Line, c(Rtest[v, ]$minC:Rtest[v, ]$maxC))
        }
        rm(v)
        
        secAll <- data.frame(idSectorCN = NA, layer = layer, cFile = Line)
        layer <- NULL
        Line  <- NULL
        
        #addition of the sector id
        secAll$idSectorCN <- paste("c", secAll$cFile, ".l", secAll$layer, sep = "")
        
        #addition of the circle
        secAll$cLdif     <- ThisC - secAll$cFile
        secAll$laydif    <- ThisL - secAll$layer
        secAll$circle    <- NA
        secAll[abs(secAll$cLdif) > abs(secAll$laydif), ]$circle  <- abs(secAll[abs(secAll$cLdif) > abs(secAll$laydif), ]$cLdif)
        secAll[abs(secAll$cLdif) <= abs(secAll$laydif), ]$circle <- abs(secAll[abs(secAll$cLdif) <= abs(secAll$laydif), ]$laydif)
        secAll[secAll$cLdif  == 0, ]$circle <- abs(secAll[secAll$cLdif == 0, ]$laydif)
        secAll[secAll$laydif == 0, ]$circle <- abs(secAll[secAll$laydif == 0, ]$cLdif)
        
        #addition of the side in regard to vessel center position
        secAll$sec <- NA
        #straight
        if(nrow(secAll[which(secAll$cLdif == 0 & secAll$laydif < 0), ]) != 0){ secAll[which(secAll$cLdif == 0 & secAll$laydif < 0), ]$sec <- "secU" }
        if(nrow(secAll[which(secAll$cLdif == 0 & secAll$laydif > 0), ]) != 0){ secAll[which(secAll$cLdif == 0 & secAll$laydif > 0), ]$sec <- "secD" }
        if(nrow(secAll[which(secAll$cLdif > 0 & secAll$laydif == 0), ]) != 0){ secAll[which(secAll$cLdif > 0 & secAll$laydif == 0), ]$sec <- "secL" }
        if(nrow(secAll[which(secAll$cLdif < 0 & secAll$laydif == 0), ]) != 0){ secAll[which(secAll$cLdif < 0 & secAll$laydif == 0), ]$sec <- "secR" }
        #diagonal
        if(nrow(secAll[which(secAll$cLdif > 0 & secAll$laydif > 0), ]) != 0){ secAll[which(secAll$cLdif > 0 & secAll$laydif > 0), ]$sec <- "secDL" }
        if(nrow(secAll[which(secAll$cLdif < 0 & secAll$laydif > 0), ]) != 0){ secAll[which(secAll$cLdif < 0 & secAll$laydif > 0), ]$sec <- "secDR" }
        if(nrow(secAll[which(secAll$cLdif > 0 & secAll$laydif < 0), ]) != 0){ secAll[which(secAll$cLdif > 0 & secAll$laydif < 0), ]$sec <- "secUL" }
        if(nrow(secAll[which(secAll$cLdif < 0 & secAll$laydif < 0), ]) != 0){ secAll[which(secAll$cLdif < 0 & secAll$laydif < 0), ]$sec <- "secUR" }
        
        #merge dataframe with position of cell within each sector # when empty sectors = NA
        ThisP           <- outD[which(outD$idSector %in% secAll$idSectorCN), c("idSector", "PolySide", "pool", "idC")]
        names(ThisP)[1] <- "idSectorCN"
        ThisP           <- base::merge(secAll, ThisP, by = "idSectorCN", all.x = TRUE)
        ThisP$tested    <- 0
        ThisP$validated <- 0
        
        ### 1.3. test of sectors
        secV <- NULL
        if(nrow(ThisP[ThisP$cLdif == 0 & ThisP$laydif == 0, ]) != 0){ ThisP[ThisP$cLdif == 0 & ThisP$laydif == 0, ]$validated <- 1 }
        
        if(any(max(ThisP$layer) == ThisL, ThisL == max(outD$layer))){ up   <- ThisL }else{ up    <- seq(ThisL + 1, max(ThisP$layer), 1) }
        if(any(min(ThisP$layer) == ThisL, ThisL == 1)){               down <- ThisL }else{ down  <- sort(seq(min(ThisP$layer), ThisL - 1, 1), decreasing = TRUE) }
        
        #ThisC will be tested during down and up loops
        if(all(max(ThisP$cFile) <= ThisC + 1, ThisC != max(outD$cFile))){
          right <- ThisC + 1 
        }else if(all(max(ThisP$cFile) <= ThisC + 1, ThisC == max(outD$cFile))){
          right <- NA
        }else{
          right <- seq(ThisC + 1, max(ThisP$cFile), 1) }
        
        if(all(min(ThisP$cFile) >= ThisC - 1, ThisC != 1)){
          left <- ThisC - 1
        }else if(all(min(ThisP$cFile) >= ThisC - 1, ThisC == 1)){
          left <- NA
        }else{
          left <- sort(seq(min(ThisP$cFile), ThisC - 1, 1), decreasing = TRUE) }
        
        #up and down sectors on ThisC
        ThisP[ThisP$cFile == ThisC, ]$tested <- 1
        secV <- NULL ; Ptest <- NULL
        
        #### down
        for(l in 1:length(down)){
          Ptest <- ThisP[ThisP$cFile == ThisC & ThisP$layer == down[l], ]
          if(all(is.na(Ptest$PolySide),
                 any(is.na(ThisP[ThisP$cFile == ThisC & ThisP$layer == down[l] + 1, c("PolySide", "sec")])),
                 unique(ThisP[ThisP$cFile == ThisC & ThisP$layer == down[l] + 1, ]$validated) == 1)){
            ThisP[ThisP$cFile == ThisC & ThisP$layer == down[l], ]$validated <- 1
            Ptest$validated <- 1
            secV           <- rbind(secV, Ptest)
            
          }else if(all(!any(qualif[qualif$sec %in% Ptest$sec, ]$location %in% Ptest$PolySide),
                       any(is.na(ThisP[ThisP$cFile == ThisC & ThisP$layer == down[l] + 1, c("PolySide", "sec")])))){
            ThisP[ThisP$cFile == ThisC & ThisP$layer == down[l], ]$validated <- 1
            Ptest$validated <- 1
            secV           <- rbind(secV, Ptest)
          }else{
            ThisP[ThisP$cFile == ThisC & ThisP$layer >= down[l], ]$tested <- 1
            break
          }
        }#end of down
        rm(l, Ptest)
        
        #### up
        Ptest <- NULL
        for(l in 1:length(up)){
          Ptest <- ThisP[ThisP$cFile == ThisC & ThisP$layer == up[l], ]
          if(all(is.na(Ptest$PolySide), 
                 any(is.na(ThisP[ThisP$cFile == ThisC & ThisP$layer == up[l] - 1, c("PolySide", "sec")])),
                 unique(ThisP[ThisP$cFile == ThisC & ThisP$layer == up[l] - 1, ]$validated) == 1)){
            ThisP[ThisP$cFile == ThisC & ThisP$layer == up[l], ]$validated <- 1
            Ptest$validated <- 1
            secV           <- rbind(secV, Ptest)
            
          }else if(all(!any(qualif[qualif$sec %in% Ptest$sec, ]$location %in% Ptest$PolySide), 
                       any(is.na(ThisP[ThisP$cFile == ThisC & ThisP$layer == up[l] - 1, c("PolySide", "sec")])))){
            ThisP[ThisP$cFile == ThisC & ThisP$layer == up[l], ]$validated <- 1
            Ptest$validated <- 1
            secV           <- rbind(secV, Ptest)
          }else{
            ThisP[ThisP$cFile == ThisC & ThisP$layer <= up[l], ]$tested <- 1
          }
        }#end of up
        rm(l, Ptest)
        
        #### right and left : different cFile
        Ptest <- NULL ; test <- NULL ; cLayerU <- NULL ; cLayerD <- NULL
        ##### right
        options(warn = 2)
        if(any(!is.na(right))){
          for(r in 1:length(right)){
            if(all(nrow(ThisP[ThisP$cFile == right[r] - 1 & ThisP$validated == 1, ]) != 0,
                   nrow(ThisP[ThisP$cFile == right[r], ]) != 0)){
              #right up
              if(ThisL <= max(ThisP[ThisP$cFile == right[r] - 1 & ThisP$validated == 1, ]$layer)){
                cLayerU <- seq(ThisL, max(ThisP[ThisP$cFile == right[r] - 1 & ThisP$validated == 1, ]$layer), 1)
                Ptest   <- NULL ; test    <- NULL
                for(l in 1:length(cLayerU)){
                  if(nrow(ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerU[l], ]) != 0){
                    ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerU[l], ]$tested  <- 1
                    Ptest <- ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerU[l], ]
                    if(all(r == 1, l == 1)){                                                                                      ## first cFile and layer
                      test <- all(unique(ThisP[ThisP$cFile == right[r] - 1 & ThisP$layer == cLayerU[l], ]$validated) == 1)                                         # sector of the previous cFile has to be validated
                      
                    }else if(all(r != 1, l == 1)){                                                                                ## first layer but not cFile
                      test <- all(unique(ThisP[ThisP$cFile == right[r] - 1 & ThisP$layer == cLayerU[l], ]$validated) == 1,                                         # sector of the previous cFile has to be validated
                                  is.na(unique(ThisP[ThisP$cFile == right[r] - 1 & ThisP$layer == cLayerU[l], ]$pool)))                                            # sector of the previous cFile should be empty
                      
                    }else if(all(right[r] == max(right), l != 1)){                                                                ## max cFile but not layer: right side
                      test <- all(unique(ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerU[l] - 1, ]$validated) == 1,                                         # below sector has to be validated
                                  unique(ThisP[ThisP$cFile == right[r] - 1 & ThisP$layer == cLayerU[l], ]$validated) == 1,                                         # sector of the previous cFile has to be validated
                                  is.na(unique(ThisP[ThisP$cFile == right[r] - 1 & ThisP$layer == cLayerU[l], ]$pool)))                                            # sector of the previous cFile should be empty
                      
                    }else if(all(right[r] != max(right), l == length(cLayerU))){                                                  ## Top side
                      test <- all(unique(ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerU[l] - 1, ]$validated) == 1,                                         # above sector has to be validated
                                  unique(ThisP[ThisP$cFile == right[r] - 1 & ThisP$layer == cLayerU[l], ]$validated) == 1,
                                  any(is.na(ThisP[ThisP$cFile == right[r] - 1 & ThisP$layer == cLayerU[l], ]$PolySide),
                                      !any(ThisP[ThisP$cFile == right[r] - 1 & ThisP$layer == cLayerU[l], ]$PolySide %in% c("LeftInf", "RightInf"))))              # cells of sector of the previous cFile should be located correctly
                      
                    }else if(all(right[r] == max(right), l == length(cLayerU))){                                                  ## Top right corner
                      test <- all(unique(ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerU[l] - 1, ]$validated) == 1,                                         # above sector has to be validated
                                  unique(ThisP[ThisP$cFile == right[r] - 1 & ThisP$layer == cLayerU[l], ]$validated) == 1,                                         # sector of the previous cFile has to be validated
                                  any(is.na(ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerU[l] - 1, ]$PolySide),
                                      all(ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerU[l] - 1, ]$PolySide %in% c("RightSup", "RightInf", "RightEdge"))), # cells of sector of the previous layer should be located correctly
                                  any(is.na(ThisP[ThisP$cFile == right[r] - 1 & ThisP$layer == cLayerU[l], ]$PolySide),
                                      !any(ThisP[ThisP$cFile == right[r] - 1 & ThisP$layer == cLayerU[l], ]$PolySide %in% c("LeftInf", "RightInf"))))              # cells of sector of the previous cFile should be located correctly
                      
                    }else{                                                                                                      ## other sectors
                      test <- all(unique(ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerU[l] - 1, ]$validated) == 1,                                         # below sector has to be validated
                                  unique(ThisP[ThisP$cFile == right[r] - 1 & ThisP$layer == cLayerU[l], ]$validated) == 1,                                         # sector of the previous cFile has to be validated
                                  is.na(unique(ThisP[ThisP$cFile == right[r] - 1 & ThisP$layer == cLayerU[l], ]$pool)),
                                  is.na(unique(ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerU[l] - 1, ]$pool)))                                            # sector of the previous layer should be empty
                    }
                    
                    if(all(is.na(Ptest$PolySide), test)){
                      ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerU[l], ]$validated <- 1
                      Ptest$validated <- 1
                      secV <- rbind(secV, Ptest)
                      
                    }else if(all(!any(qualif[qualif$sec %in% Ptest$sec, ]$location %in% Ptest$PolySide), test)){
                      ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerU[l], ]$validated <- 1
                      Ptest$validated <- 1
                      secV <- rbind(secV, Ptest)
                    }else{
                      ThisP[ThisP$cFile == right[r] & ThisP$layer <= cLayerU[l], ]$tested <- 1
                      break
                    }
                  }
                }
                rm(l, Ptest, test)
              }
              
              #right down
              Ptest <- NULL ; test <- NULL
              if(min(ThisP[ThisP$cFile == right[r] - 1 & ThisP$validated == 1, ]$layer) < ThisL){
                cLayerD <- sort(seq(min(ThisP[ThisP$cFile == right[r] - 1 & ThisP$validated == 1, ]$layer), ThisL - 1,  1), decreasing = TRUE)
                Ptest   <- NULL
                for(l in 1:length(cLayerD)){
                  if(nrow(ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerD[l], ]) != 0){
                    ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerD[l], ]$tested  <- 1
                    Ptest <- ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerD[l], ]
                    
                    if(all(right[r] == max(right), l != length(cLayerD))){                                                       ## max cFile: right side
                      test <- all(unique(ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerD[l] + 1, ]$validated) == 1,                                         # above sector has to be validated
                                  unique(ThisP[ThisP$cFile == right[r] - 1 & ThisP$layer == cLayerD[l], ]$validated) == 1,                                         # sector of the previous cFile has to be validated
                                  is.na(unique(ThisP[ThisP$cFile == right[r] - 1 & ThisP$layer == cLayerD[l], ]$pool)))                                            # sector of the previous cFile should be empty
                      
                    }else if(all(right[r] != max(right), l == length(cLayerD))){                                                 ## last layer: bottom layer     
                      test <- all(unique(ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerD[l] + 1, ]$validated) == 1,                                         # above sector has to be validated
                                  unique(ThisP[ThisP$cFile == right[r] - 1 & ThisP$layer == cLayerD[l], ]$validated) == 1,                                            # sector of the previous layer should be empty
                                  any(is.na(ThisP[ThisP$cFile == right[r] - 1 & ThisP$layer == cLayerD[l], ]$PolySide),
                                      !any(ThisP[ThisP$cFile == right[r] - 1 & ThisP$layer == cLayerD[l], ]$PolySide %in% c("LeftSup", "RightSup"))))              # cells of sector of the previous cFile should be located correctly
                      
                    }else if(all(right[r] == max(right), l == length(cLayerD))){                                                 ## Bottom right corner     
                      test <- all(unique(ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerD[l] + 1, ]$validated) == 1,                                         # above sector has to be validated
                                  unique(ThisP[ThisP$cFile == right[r] - 1 & ThisP$layer == cLayerD[l], ]$validated) == 1,                                         # sector of the previous cFile has to be validated
                                  any(is.na(ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerD[l] + 1, ]$PolySide),
                                      all(ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerD[l] + 1, ]$PolySide %in% c("RightSup", "RightInf", "RightEdge"))), # cells of sector of the previous layer should be located correctly
                                  any(is.na(ThisP[ThisP$cFile == right[r] - 1 & ThisP$layer == cLayerD[l], ]$PolySide),
                                      !any(ThisP[ThisP$cFile == right[r] - 1 & ThisP$layer == cLayerD[l], ]$PolySide %in% c("LeftSup", "RightSup"))))              # cells of sector of the previous cFile should be located correctly
                      
                    }else{                                                                                                       ## other sectors
                      test <- all(unique(ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerD[l] + 1, ]$validated) == 1,                                         # above sector has to be validated
                                  unique(ThisP[ThisP$cFile == right[r] - 1 & ThisP$layer == cLayerD[l], ]$validated) == 1,                                         # sector of the previous cFile has to be validated
                                  is.na(unique(ThisP[ThisP$cFile == right[r] - 1 & ThisP$layer == cLayerD[l], ]$pool)),
                                  is.na(unique(ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerD[l] + 1, ]$pool)))                                            # sector of the previous layer should be empty
                    }
                    
                    if(all(is.na(Ptest$PolySide), test)){
                      ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerD[l], ]$validated <- 1
                      Ptest$validated <- 1
                      secV <- rbind(secV, Ptest)
                      
                    }else if(all(!any(qualif[qualif$sec %in% Ptest$sec, ]$location %in% Ptest$PolySide), test)){
                      ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerD[l], ]$validated <- 1
                      Ptest$validated <- 1
                      secV <- rbind(secV, Ptest)
                    }else{
                      ThisP[ThisP$cFile == right[r] & ThisP$layer >= cLayerD[l], ]$tested <- 1
                      break
                    }
                  }
                }
                rm(l, Ptest)
              }#end of right down
              
            }else{
              break
            }
          } #end of right
          rm(r)
        }
        
        ##### left
        if(any(!is.na(left))){
          Ptest <- NULL ; test <- NULL ; cLayerU <- NULL ; cLayerD <- NULL
          for(r in 1:length(left)){
            if(all(nrow(ThisP[ThisP$cFile == left[r] + 1 & ThisP$validated == 1, ]) != 0,
                   nrow(ThisP[ThisP$cFile == left[r], ]) != 0)){
              #left up
              if(ThisL <= max(ThisP[ThisP$cFile == left[r] + 1 & ThisP$validated == 1, ]$layer)){
                cLayerU <- seq(ThisL, max(ThisP[ThisP$cFile == left[r] + 1 & ThisP$validated == 1, ]$layer), 1)
                Ptest   <- NULL ; test <- NULL
                
                for(l in 1:length(cLayerU)){
                  if(nrow(ThisP[ThisP$cFile == left[r] & ThisP$layer == cLayerU[l], ]) != 0){
                    ThisP[ThisP$cFile == left[r] & ThisP$layer == cLayerU[l], ]$tested  <- 1
                    Ptest <- ThisP[ThisP$cFile == left[r] & ThisP$layer == cLayerU[l], ]
                    
                    if(all(r == 1, l == 1)){                                                                                    ## first cFile and layer
                      test <- all(unique(ThisP[ThisP$cFile == left[r] + 1 & ThisP$layer == cLayerU[l], ]$validated) == 1)                                         # sector of the previous cFile has to be validated
                      
                    }else if(all(r != 1, l == 1)){                                                                              ## first layer but not cFile
                      test <- all(unique(ThisP[ThisP$cFile == left[r] + 1 & ThisP$layer == cLayerU[l], ]$validated) == 1,                                         # sector of the previous cFile has to be validated
                                  is.na(unique(ThisP[ThisP$cFile == left[r] + 1 & ThisP$layer == cLayerU[l], ]$pool)))                                            # sector of the previous cFile should be empty
                      
                    }else if(all(left[r] %in% min(left), l != 1)){                                                              ## min cFile but not layer: left side
                      test <- all(unique(ThisP[ThisP$cFile == left[r] & ThisP$layer == cLayerU[l] - 1, ]$validated) == 1,                                         # below sector has to be validated
                                  unique(ThisP[ThisP$cFile == left[r] + 1 & ThisP$layer == cLayerU[l], ]$validated) == 1,                                         # sector of the previous cFile has to be validated
                                  is.na(unique(ThisP[ThisP$cFile == left[r] + 1 & ThisP$layer == cLayerU[l], ]$pool)))                                            # sector of the previous cFile should be empty
                      
                    }else if(all(left[r] != min(left), l == length(cLayerU))){                                                  ## Top side
                      test <- all(unique(ThisP[ThisP$cFile == left[r] & ThisP$layer == cLayerU[l] - 1, ]$validated) == 1,                                         # above sector has to be validated
                                  unique(ThisP[ThisP$cFile == left[r] + 1 & ThisP$layer == cLayerU[l], ]$validated) == 1,
                                  any(is.na(ThisP[ThisP$cFile == left[r] + 1 & ThisP$layer == cLayerU[l], ]$PolySide),
                                      !any(ThisP[ThisP$cFile == left[r] + 1 & ThisP$layer == cLayerU[l], ]$PolySide %in% c("LeftInf", "RightInf"))))              # cells of sector of the previous cFile should be located correctly
                      
                    }else if(all(left[r] == min(left), l == length(cLayerU))){                                                 ## Top left corner
                      test <- all(unique(ThisP[ThisP$cFile == left[r] & ThisP$layer == cLayerU[l] - 1, ]$validated) == 1,                                         # above sector has to be validated
                                  unique(ThisP[ThisP$cFile == left[r] + 1 & ThisP$layer == cLayerU[l], ]$validated) == 1,                                         # sector of the previous cFile has to be validated
                                  any(is.na(ThisP[ThisP$cFile == left[r] & ThisP$layer == cLayerU[l] - 1, ]$PolySide),
                                      all(ThisP[ThisP$cFile == left[r] & ThisP$layer == cLayerU[l] - 1, ]$PolySide %in% c("LeftSup", "LeftInf", "LeftEdge"))),    # cells of sector of the previous layer should be located correctly
                                  any(is.na(ThisP[ThisP$cFile == left[r] + 1 & ThisP$layer == cLayerU[l], ]$PolySide),
                                      !any(ThisP[ThisP$cFile == left[r] + 1 & ThisP$layer == cLayerU[l], ]$PolySide %in% c("LeftInf", "RightInf"))))              # cells of sector of the previous cFile should be located correctly
                      
                    }else{                                                                                                      ## other sectors
                      test <- all(unique(ThisP[ThisP$cFile == left[r] & ThisP$layer == cLayerU[l] - 1, ]$validated) == 1,                                         # below sector has to be validated
                                  unique(ThisP[ThisP$cFile == left[r] + 1 & ThisP$layer == cLayerU[l], ]$validated) == 1,                                         # sector of the previous cFile has to be validated
                                  is.na(unique(ThisP[ThisP$cFile == left[r] + 1 & ThisP$layer == cLayerU[l], ]$pool)),
                                  is.na(unique(ThisP[ThisP$cFile == left[r] & ThisP$layer == cLayerU[l] - 1, ]$pool)))                                            # sector of the previous layer should be empty
                    }
                    
                    if(all(is.na(Ptest$PolySide), test)){
                      ThisP[ThisP$cFile == left[r] & ThisP$layer == cLayerU[l], ]$validated <- 1
                      Ptest$validated <- 1
                      secV <- rbind(secV, Ptest)
                      
                    }else if(all(!any(qualif[qualif$sec %in% Ptest$sec, ]$location %in% Ptest$PolySide), test)){
                      ThisP[ThisP$cFile == left[r] & ThisP$layer == cLayerU[l], ]$validated <- 1
                      Ptest$validated <- 1
                      secV <- rbind(secV, Ptest)
                    }else{
                      ThisP[ThisP$cFile == left[r] & ThisP$layer <= cLayerU[l], ]$tested <- 1
                      break
                    }
                  }
                }
                rm(l, Ptest, test)
              }
              
              #left down
              testL <- NULL ; testC <- NULL
              if(min(ThisP[ThisP$cFile == left[r] + 1 & ThisP$validated == 1, ]$layer) < ThisL){
                cLayerD <- sort(seq(min(ThisP[ThisP$cFile == left[r] + 1 & ThisP$validated == 1, ]$layer), ThisL - 1,  1), decreasing = TRUE)
                Ptest   <- NULL
                for(l in 1:length(cLayerD)){
                  if(nrow(ThisP[ThisP$cFile == left[r] & ThisP$layer == cLayerD[l], ]) != 0){
                    ThisP[ThisP$cFile == left[r] & ThisP$layer == cLayerD[l], ]$tested  <- 1
                    Ptest <- ThisP[ThisP$cFile == left[r] & ThisP$layer == cLayerD[l], ]
                    
                    if(all(left[r] %in% c(min(left), min(left) + 1), l != length(cLayerD))){                                    ## min cFile: left side
                      test <- all(unique(ThisP[ThisP$cFile == left[r] & ThisP$layer == cLayerD[l] + 1, ]$validated) == 1,                                          # above sector has to be validated
                                  unique(ThisP[ThisP$cFile == left[r] + 1 & ThisP$layer == cLayerD[l], ]$validated) == 1,                                          # sector of the previous cFile has to be validated
                                  is.na(unique(ThisP[ThisP$cFile == left[r] + 1 & ThisP$layer == cLayerD[l], ]$pool)))                                             # sector of the previous cFile should be empty
                      
                    }else if(all(left[r] != min(left), l == length(cLayerD))){                                                  ## last layer: bottom layer     
                      test <- all(unique(ThisP[ThisP$cFile == left[r] & ThisP$layer == cLayerD[l] + 1, ]$validated) == 1,                                          # above sector has to be validated
                                  unique(ThisP[ThisP$cFile == left[r] + 1 & ThisP$layer == cLayerD[l], ]$validated) == 1,
                                  any(is.na(ThisP[ThisP$cFile == left[r] + 1 & ThisP$layer == cLayerD[l], ]$PolySide),
                                      !any(ThisP[ThisP$cFile == left[r] + 1 & ThisP$layer == cLayerD[l], ]$PolySide %in% c("LeftSup", "RightSup"))))               # cells of sector of the previous cFile should be located correctly
                      
                    }else if(all(left[r] == min(left), l == length(cLayerD))){                                                  ## Bottom left corner
                      test <- all(unique(ThisP[ThisP$cFile == left[r] & ThisP$layer == cLayerD[l] + 1, ]$validated) == 1,                                          # above sector has to be validated
                                  unique(ThisP[ThisP$cFile == left[r] + 1 & ThisP$layer == cLayerD[l], ]$validated) == 1,                                          # sector of the previous cFile has to be validated
                                  any(is.na(ThisP[ThisP$cFile == right[r] & ThisP$layer == cLayerD[l] + 1, ]$PolySide),
                                      all(ThisP[ThisP$cFile == left[r] & ThisP$layer == cLayerD[l] + 1, ]$PolySide %in% c("LeftSup", "LeftInf", "LeftEdge"))),     # cells of sector of the previous layer should be located correctly
                                  any(is.na(ThisP[ThisP$cFile == left[r] + 1 & ThisP$layer == cLayerD[l], ]$PolySide),
                                      !any(ThisP[ThisP$cFile == left[r] + 1 & ThisP$layer == cLayerD[l], ]$PolySide %in% c("LeftSup", "RightSup"))))               # cells of sector of the previous cFile should be located correctly
                      
                    }else{                                                                                                      ## other sectors
                      test <- all(unique(ThisP[ThisP$cFile == left[r] & ThisP$layer == cLayerD[l] + 1, ]$validated) == 1,                                          # below sector has to be validated
                                  unique(ThisP[ThisP$cFile == left[r] + 1 & ThisP$layer == cLayerD[l], ]$validated) == 1,                                          # sector of the previous cFile has to be validated
                                  is.na(unique(ThisP[ThisP$cFile == left[r] & ThisP$layer == cLayerD[l] + 1, ]$pool)))                                             # sector of the previous cFile should be empty
                    }
                    
                    if(all(is.na(Ptest$PolySide), test)){
                      ThisP[ThisP$cFile == left[r] & ThisP$layer == cLayerD[l], ]$validated <- 1
                      Ptest$validated <- 1
                      secV <- rbind(secV, Ptest)
                      
                    }else if(all(!any(qualif[qualif$sec %in% Ptest$sec, ]$location %in% Ptest$PolySide), test)){
                      ThisP[ThisP$cFile == left[r] & ThisP$layer == cLayerD[l], ]$validated <- 1
                      Ptest$validated <- 1
                      secV <- rbind(secV, Ptest)
                    }else{
                      ThisP[ThisP$cFile == left[r] & ThisP$layer >= cLayerD[l], ]$tested <- 1
                      break
                    }
                  }
                }
                rm(l, Ptest)
              }#end of left down
            }else{
              break
            }
          } #end of left
          rm(r)
        }
        
      }
    
    # 2. formatting results ####
    ## extraction of the identified vessel sectors
    outV <- EdgeS[EdgeS$idSector %in% sort(c(paste("c", ThisC, ".l", ThisL, sep = ""), unique(secV$idSectorCN))), ]
    
    ## formatting
    outV$tested       <- 1
    outV$validated    <- 1
    outV$Vorigin      <- NA
    outV[outV$idSector == paste("c", ThisC, ".l", ThisL, sep = ""), ]$Vorigin <- 1
    
    outV$idVessel <- unique(ThisN$idV)
    outV$cell     <- 1
    outV$vessel   <- 1 
    outV$cLdif    <- ThisC - outV$cFile
    outV$laydif   <- ThisL - outV$layer
    outV$circle   <- NA
    outV[abs(outV$cLdif) > abs(outV$laydif), ]$circle  <- abs(outV[abs(outV$cLdif) > abs(outV$laydif), ]$cLdif)
    outV[abs(outV$cLdif) <= abs(outV$laydif), ]$circle <- abs(outV[abs(outV$cLdif) <= abs(outV$laydif), ]$laydif)
    outV[outV$cLdif == 0, ]$circle  <- abs(outV[outV$cLdif == 0, ]$laydif)
    outV[outV$laydif == 0, ]$circle <- abs(outV[outV$laydif == 0, ]$cLdif)
    
    ## saving
    outVtot <- rbind(outVtot, outV)
  } ## end loop j
  
  #last formatting
  outVtot$cell      <- NULL
  outVtot$tested    <- NULL
  outVtot$validated <- NULL
  
  rm(j)
  #end of the decision tree cell or not cell
  
  # 3. plot for checking #####
  if(make.plot == TRUE){
    title <- paste(as.character(unique(dataset$idD)), " vessels are in the place", sep = "")
    
    plotV <- plot(Y.cor ~ X.cor, dataset, main =  title, cex.main = 0.8, type = "n",
                  xlim = c(-30, max(EdgeS$xRightSupS, EdgeS$xRightInfS, na.rm = TRUE) + 5),
                  ylim = c(-30, max(dataset$Y.cor, na.rm = TRUE) * 1.02))
    abline(v = seq(-10, max(EdgeS$xRightSupS, EdgeS$xRightInfS, na.rm = TRUE) + 10, 10), col = "grey", lty = 3)
    abline(h = seq(0, max(dataset$Y.cor, na.rm = TRUE), 50), col = "grey", lty = 3)
    
    for(j in 1:length(idSector)){
      if(EdgeS[EdgeS$idSector == idSector[j],]$APF == 1){
        polygon(c(EdgeS[EdgeS$idSector == idSector[j], ]$xLeftInfP,  EdgeS[EdgeS$idSector == idSector[j], ]$xLeftEdgeS,  EdgeS[EdgeS$idSector == idSector[j], ]$xLeftSupP,
                  EdgeS[EdgeS$idSector == idSector[j], ]$xRightSupS, EdgeS[EdgeS$idSector == idSector[j], ]$xRightEdgeS, EdgeS[EdgeS$idSector == idSector[j], ]$xRightInfS),
                c(EdgeS[EdgeS$idSector == idSector[j], ]$yLeftInfS,  EdgeS[EdgeS$idSector == idSector[j], ]$yLeftEdgeS,  EdgeS[EdgeS$idSector == idSector[j], ]$yLeftSupS, 
                  EdgeS[EdgeS$idSector == idSector[j], ]$yRightSupS, EdgeS[EdgeS$idSector == idSector[j], ]$yRightEdgeS, EdgeS[EdgeS$idSector == idSector[j], ]$yRightInfS),
                col = adjustcolor(col = "maroon2", alpha.f = 0.1))
      }
      
      if(EdgeS[EdgeS$idSector == idSector[j],]$vessel == 1){
        polygon(c(EdgeS[EdgeS$idSector == idSector[j], ]$xLeftInfP,  EdgeS[EdgeS$idSector == idSector[j], ]$xLeftEdgeS,  EdgeS[EdgeS$idSector == idSector[j], ]$xLeftSupP,
                  EdgeS[EdgeS$idSector == idSector[j], ]$xRightSupS, EdgeS[EdgeS$idSector == idSector[j], ]$xRightEdgeS, EdgeS[EdgeS$idSector == idSector[j], ]$xRightInfS),
                c(EdgeS[EdgeS$idSector == idSector[j], ]$yLeftInfS,  EdgeS[EdgeS$idSector == idSector[j], ]$yLeftEdgeS,  EdgeS[EdgeS$idSector == idSector[j], ]$yLeftSupS, 
                  EdgeS[EdgeS$idSector == idSector[j], ]$yRightSupS, EdgeS[EdgeS$idSector == idSector[j], ]$yRightEdgeS, EdgeS[EdgeS$idSector == idSector[j], ]$yRightInfS),
                col = adjustcolor(col = "turquoise4", alpha.f = 0.1))
      }
      
      polygon(c(EdgeS[EdgeS$idSector == idSector[j], ]$xLeftInfP,  EdgeS[EdgeS$idSector == idSector[j], ]$xLeftEdgeS,  EdgeS[EdgeS$idSector == idSector[j], ]$xLeftSupP,
                EdgeS[EdgeS$idSector == idSector[j], ]$xRightSupS, EdgeS[EdgeS$idSector == idSector[j], ]$xRightEdgeS, EdgeS[EdgeS$idSector == idSector[j], ]$xRightInfS),
              c(EdgeS[EdgeS$idSector == idSector[j], ]$yLeftInfS,  EdgeS[EdgeS$idSector == idSector[j], ]$yLeftEdgeS,  EdgeS[EdgeS$idSector == idSector[j], ]$yLeftSupS, 
                EdgeS[EdgeS$idSector == idSector[j], ]$yRightSupS, EdgeS[EdgeS$idSector == idSector[j], ]$yRightEdgeS, EdgeS[EdgeS$idSector == idSector[j], ]$yRightInfS),
              border = "wheat2", lwd = 2)
    }
    rm(j)
    
    for(v in 1:nrow(outVtot)){
      polygon(c(outVtot[v, ]$xLeftInfP,  outVtot[v, ]$xLeftEdgeS,  outVtot[v, ]$xLeftSupP,
                outVtot[v, ]$xRightSupS, outVtot[v, ]$xRightEdgeS, outVtot[v, ]$xRightInfS),
              c(outVtot[v, ]$yLeftInfS,  outVtot[v, ]$yLeftEdgeS,  outVtot[v, ]$yLeftSupS, 
                outVtot[v, ]$yRightSupS, outVtot[v, ]$yRightEdgeS, outVtot[v, ]$yRightInfS),
              col = adjustcolor(col = "darkturquoise", alpha.f = 0.1),
              border = "wheat2", lwd = 2)
    }
    rm(v)
    
    points(Y.cor ~ X.cor, dataset[dataset$type == "Vessel", ], col = "darkturquoise", pch = 16)
    points(Y.cor ~ X.cor, dataset[dataset$type == "Fibers", ], col = "brown", pch = 16)  
    points(Y.cor ~ X.cor, dataset[dataset$type == "Axial Parenchyma", ], col = "purple2", pch = 15, cex = 0.8)
    points(Y.cor ~ X.cor, dataset[dataset$type == "APF", ], col = "maroon2", pch = 16)
    points(Y.cor ~ X.cor, dataset[dataset$type == "Cambial", ], col = "seagreen3", pch = 15)
    
    # idVessel label
    secV <- aggregate(data.frame(idVessel = outVtot$idVessel), list(idSector = outVtot$idSector), unique)
    
    for(v in 1:nrow(secV)){
      text((outVtot[outVtot$idSector == secV[v, ]$idSector, ]$xRightInfS + outVtot[outVtot$idSector == secV[v, ]$idSector, ]$xLeftInfP + outVtot[outVtot$idSector == secV[v, ]$idSector, ]$xRightSupS + outVtot[outVtot$idSector == secV[v, ]$idSector, ]$xLeftSupP) / 4,
           (outVtot[outVtot$idSector == secV[v, ]$idSector, ]$yRightInfS + outVtot[outVtot$idSector == secV[v, ]$idSector, ]$yLeftInfS + outVtot[outVtot$idSector == secV[v, ]$idSector, ]$yRightSupS + outVtot[outVtot$idSector == secV[v, ]$idSector, ]$yLeftSupS) / 4,
           secV[v, ]$idVessel, cex = 0.6)
    }
    
    # layer label
    layer <- unique(EdgeS$layer)
    for(j in 1:length(layer)){
      text(-29, (EdgeS[EdgeS$layer == layer[j] & EdgeS$cFile == 1, ]$yLeftSupS + EdgeS[EdgeS$layer == layer[j] & EdgeS$cFile == 1, ]$yLeftInfS) / 2,
           layer[j], cex = 0.6, font = 4, col = "black")
      
    }
    rm(j)
    
    # cambial line label
    Line <- unique(EdgeS$cFile)
    for(j in 1:length(Line)){
      text((EdgeS[EdgeS$layer == 1 & EdgeS$cFile == Line[j], ]$xLeftInfP + EdgeS[EdgeS$layer == 1 & EdgeS$cFile == Line[j], ]$xRightInfS) / 2, -26,
           Line[j], cex = 0.6, font = 4, col = "black")
      
    }
    rm(j)
    
    rm(plotV, title)
  }
  
  }else{
    outVtot <- NA
  }
  
    # 4. to export ####
    return(outVtot)
}

##### END OF THE FUNCTION vesselS #####


#################### FUNCTION rrf #####
# function to attribute to each cell a sector id a calculation of sector cell characteristic
# return list of 1 object: dfCal

#' @dataset: dataframe with reoriented data (X.cor, Y.cor) and "type" column
#' @EdgeS: dataframe of the sector edges
#' @Vsec: dataframe of vessel sectors

rrf <- function(dataset, EdgeS, Vsec, make.plot = TRUE){
  require(reshape2) || install.packages(reshape2)
  library(reshape2)
  require(dplyr) || install.packages(dplyr)
  library(dplyr)
  
  # vector and dataframe
  work     <- EdgeS
  secIm    <- unique(work$idSector)
  idSector <- unique(work$idSector)
  
  ## Phase per layer
  # phase of the sector
  ThisP           <- aggregate(dataset$phase, list(idSector = dataset$idSector, layer = dataset$layer, type = dataset$Cells, phase = dataset$phase), length)
  names(ThisP)[5] <- "nbCell"
  ThisP           <- ThisP[grep("Fibers", ThisP$type), ]
  secPh           <- reshape2::dcast(ThisP, idSector + layer ~ phase, value.var = "nbCell")
  if(any(ThisP$phase %in% "E")){ if(any(is.na(secPh$E))){ secPh[is.na(secPh$E), ]$E <- 0 }}else{ secPh$E <- 0 }
  if(any(ThisP$phase %in% "S")){ if(any(is.na(secPh$S))){ secPh[is.na(secPh$S), ]$S <- 0 }}else{ secPh$S <- 0 }
  if(any(ThisP$phase %in% "M")){ if(any(is.na(secPh$M))){ secPh[is.na(secPh$M), ]$M <- 0 }}else{ secPh$M <- 0 }
  rm(ThisP)
  
  secPh$phaseF <- NA
  if(any(secPh$E > (secPh$S + secPh$M))){ secPh[secPh$E > (secPh$S + secPh$M),]$phaseF  <- "E" }
  if(any(secPh$E <= (secPh$S + secPh$M))){ secPh[secPh$E <= (secPh$S + secPh$M),]$phaseF <- "SM" }
  
  work <- merge(work, secPh[, c("idSector", "phaseF")], by = "idSector", all.x = TRUE)
  
  
  # phase per layer
  ThisL <- reshape2::dcast(aggregate(data.frame(nbPoly = secPh$phaseF), list(layer = secPh$layer, phaseF = secPh$phaseF), length),
                           layer ~ phaseF, value.var = "nbPoly")
  if(any(colnames(ThisL) %in% "E")){ if(any(is.na(ThisL$E))){ ThisL[is.na(ThisL$E), ]$E <- 0 }}else{ ThisL$E <- 0 }
  if(any(colnames(ThisL) %in% "SM")){ if(any(is.na(ThisL$SM))){ ThisL[is.na(ThisL$SM), ]$SM <- 0 }}else{ ThisL$SM <- 0 }
  
  ThisL$phaseLayF <- NA
  if(any(ThisL$E > ThisL$SM)){ ThisL[ThisL$E > ThisL$SM, ]$phaseLayF  <- "E" }
  if(any(ThisL$E <= ThisL$SM)){ ThisL[ThisL$E <= ThisL$SM,]$phaseLayF <- "SM" }
  
  work <- merge(work, ThisL[, c("layer", "phaseLayF")], by = "layer", all.x = TRUE)
  rm(secPh)
  
  ## Calculation
  # sector with only APF
  if(all(is.na(Vsec))){
    secAPF <- secIm
    secVA  <- NULL
    secVs  <- NULL
    secVl  <- NULL
  }else{
    secVs  <- unique(Vsec$idSector)                  # all vessel sectors
    secAPF <- secIm[!(secIm %in% secVs)]            # sectors with only APF
    secVA  <- unique(Vsec[Vsec$APF == 1,]$idSector) # sectors with vessels and APF
    if(length(secVA) != 0){                         # sector with only vessel
      secVl <- secVs[!(secVs %in% secVA)]
    } else {
      secVl <- secVs
    }
  }
  
  # test to add weight to sector regarding type of cell. APF = 1, APF+Vessel = 0.5, Vessel = 0
  work$weight <- NA
  work[work$idSector %in% secAPF,]$weight <- 1
  if(!any(!is.null(secVl), length(secVl) == 0)){ work[work$idSector %in% secVl,]$weight  <- 0 }
  if(length(secVA) != 0){ work[work$idSector %in% secVA,]$weight <- 0.5 }
  
  work$nbcellW <- work$nbcell^work$weight
  
  # dataframe without sectors with vessel
  withoutV  <- work[work$idSector %in% secIm[!(secIm %in% secVs)],] # without all vessel sectors (even with APF)
  withoutVw <- work[work$idSector %in% secIm[!(secIm %in% secVl)],] # without vessel sectors (not removing APF+Vessel sectors)
  
  # 1. number of cells #####
  ## sum of cell per cambial line
  sumClCL <- aggregate(data.frame(CellCLsum = work$nbcell), list(cFile = work$cFile), sum)
  
  ## weight mean (with weight column) = RRF
  wavgLay <- as.data.frame(withoutVw %>% group_by(layer) %>% summarise(RRF = weighted.mean(nbcell, weight), .groups = 'drop'))
  
  # 2. formatting output dataframes #####
  ## addition of the phase
  dfCal <- data.frame(idD = unique(dataset$idD), id = unique(dataset$id),
                      year = unique(dataset$year), DOY = unique(dataset$DOY), date = unique(dataset$date),
                      layer = c(unique(wavgLay$layer), rep(NA, length(sumClCL$cFile))),
                      cFile = c(rep(NA, length(wavgLay$layer)), unique(sumClCL$cFile)))
  
  dfCal <- merge(dfCal, wavgLay, by = "layer", all = TRUE)
  dfCal <- merge(dfCal, sumClCL, by = "cFile", all = TRUE)
  dfCal <- merge(dfCal, ThisL[, c("layer", "phaseLayF")], by = "layer", all = TRUE)
  
  dfCal <- dfCal[, c("idD", "id", "year", "DOY", "date", "layer", "cFile", "phaseLayF",
                     "CellCLsum", "RRF")]
  
  ## dataframe number of producted cell per method
  dfSummary <- data.frame(idD = unique(dataset$idD), id = unique(dataset$id),
                          year = unique(dataset$year), DOY = unique(dataset$DOY), date = unique(dataset$date),
                          CellCLsum_avg = NA, RRF_avg = NA)
  
  # number of producted APF per cFile following the different methods
  dfSummary$CellCLsum_avg <- ceiling(mean(sumClCL$CellCLsum))
  dfSummary$RRF_avg       <- ceiling(sum(wavgLay$RRF))
  
  # 3. plot for checking #####
  if(make.plot == TRUE){
    #### plot sector ####
    title <- paste(as.character(unique(dataset$idD)), " RRF", sep = "")
    layer <- sort(unique(work$layer))
    Line  <- sort(unique(work$cFile))
    
    pRRF <- plot(Y.cor ~ X.cor, dataset, main =  title, cex.main = 0.8, type = "n",
                  xlim = c(-10, max(work$xRightSupS, work$xRightInfS, na.rm = TRUE) + 25),
                  ylim = c(-10, max(work$yRightSupS, na.rm = TRUE) * 1.05))
    abline(v = seq(-10, max(work$xRightSupS, work$xRightInfS, na.rm = TRUE) + 30, 10), col = "grey", lty = 3)
    abline(h = seq(0, max(dataset$Y.cor, na.rm = TRUE), 50), col = "grey", lty = 3)
    
    for(j in 1:length(idSector)){
      if(work[work$idSector == idSector[j],]$APF == 1){
        polygon(c(work[work$idSector == idSector[j], ]$xLeftInfP,  work[work$idSector == idSector[j], ]$xLeftwork,  work[work$idSector == idSector[j], ]$xLeftSupP,
                  work[work$idSector == idSector[j], ]$xRightSupS, work[work$idSector == idSector[j], ]$xRightwork, work[work$idSector == idSector[j], ]$xRightInfS),
                c(work[work$idSector == idSector[j], ]$yLeftInfS,  work[work$idSector == idSector[j], ]$yLeftwork,  work[work$idSector == idSector[j], ]$yLeftSupS, 
                  work[work$idSector == idSector[j], ]$yRightSupS, work[work$idSector == idSector[j], ]$yRightwork, work[work$idSector == idSector[j], ]$yRightInfS),
                col = adjustcolor(col = "maroon2", alpha.f = 0.1))
      }
      
      if(work[work$idSector == idSector[j],]$vessel == 1){
        polygon(c(work[work$idSector == idSector[j], ]$xLeftInfP,  work[work$idSector == idSector[j], ]$xLeftwork,  work[work$idSector == idSector[j], ]$xLeftSupP,
                  work[work$idSector == idSector[j], ]$xRightSupS, work[work$idSector == idSector[j], ]$xRightwork, work[work$idSector == idSector[j], ]$xRightInfS),
                c(work[work$idSector == idSector[j], ]$yLeftInfS,  work[work$idSector == idSector[j], ]$yLeftwork,  work[work$idSector == idSector[j], ]$yLeftSupS, 
                  work[work$idSector == idSector[j], ]$yRightSupS, work[work$idSector == idSector[j], ]$yRightwork, work[work$idSector == idSector[j], ]$yRightInfS),
                col = adjustcolor(col = "turquoise4", alpha.f = 0.1))
      }
      
      polygon(c(work[work$idSector == idSector[j], ]$xLeftInfP,  work[work$idSector == idSector[j], ]$xLeftwork,  work[work$idSector == idSector[j], ]$xLeftSupP,
                work[work$idSector == idSector[j], ]$xRightSupS, work[work$idSector == idSector[j], ]$xRightwork, work[work$idSector == idSector[j], ]$xRightInfS),
              c(work[work$idSector == idSector[j], ]$yLeftInfS,  work[work$idSector == idSector[j], ]$yLeftwork,  work[work$idSector == idSector[j], ]$yLeftSupS, 
                work[work$idSector == idSector[j], ]$yRightSupS, work[work$idSector == idSector[j], ]$yRightwork, work[work$idSector == idSector[j], ]$yRightInfS),
              border = "wheat2", lwd = 2)
    }
    rm(j)
    
    ## Vessels
    if(!all(is.na(Vsec))){
      # Vessel sector
      for(v in 1:nrow(Vsec)){
        polygon(c(Vsec[v, ]$xLeftInfP,  Vsec[v, ]$xLeftwork,  Vsec[v, ]$xLeftSupP,
                  Vsec[v, ]$xRightSupS, Vsec[v, ]$xRightwork, Vsec[v, ]$xRightInfS),
                c(Vsec[v, ]$yLeftInfS,  Vsec[v, ]$yLeftwork,  Vsec[v, ]$yLeftSupS, 
                  Vsec[v, ]$yRightSupS, Vsec[v, ]$yRightwork, Vsec[v, ]$yRightInfS),
                col = adjustcolor(col = "darkturquoise", alpha.f = 0.1),
                border = "wheat2", lwd = 2)
      }
      rm(v)
      
      # idVessel label
      labV <- aggregate(data.frame(idVessel = Vsec$idVessel), list(idSector = Vsec$idSector), unique)
      for(v in 1:nrow(labV)){
        text((Vsec[Vsec$idSector == labV[v, ]$idSector, ]$xRightInfS + Vsec[Vsec$idSector == labV[v, ]$idSector, ]$xLeftInfP + Vsec[Vsec$idSector == labV[v, ]$idSector, ]$xRightSupS + Vsec[Vsec$idSector == labV[v, ]$idSector, ]$xLeftSupP) / 4,
             (Vsec[Vsec$idSector == labV[v, ]$idSector, ]$yRightInfS + Vsec[Vsec$idSector == labV[v, ]$idSector, ]$yLeftInfS + Vsec[Vsec$idSector == labV[v, ]$idSector, ]$yRightSupS + Vsec[Vsec$idSector == labV[v, ]$idSector, ]$yLeftSupS) / 4,
             labV[v, ]$idVessel, cex = 0.6)
      }
      rm(v, labV) 
    } 
    
    ## cells
    points(Y.cor ~ X.cor, dataset[dataset$type == "Vessel",], col = "darkturquoise", pch = 16)
    points(Y.cor ~ X.cor, dataset[dataset$type == "Fibers",], col = "brown", pch = 16)  
    points(Y.cor ~ X.cor, dataset[dataset$type == "Axial Parenchyma",], col = "purple2", pch = 15, cex = 0.8)
    points(Y.cor ~ X.cor, dataset[dataset$type == "APF" & dataset$phase == "E",], col = "maroon2", pch = 1)
    points(Y.cor ~ X.cor, dataset[dataset$type == "APF" & dataset$phase == "S",], col = "maroon2", pch = 16)
    points(Y.cor ~ X.cor, dataset[dataset$type == "Cambial",], col = "seagreen3", pch = 15)
    
    ## layer label and number of cell of the referential radial line
    for(j in 1:length(layer)){
      # layer label
      text(-11, (work[work$layer == layer[j] & work$cFile == 1, ]$yLeftSupS + work[work$layer == layer[j] & work$cFile == 1, ]$yLeftInfS) / 2,
           layer[j], cex = 1, col = "grey60")
      
      # weighted average number of APF
      text(max(work$xRightSupS, work$xRightInfS, na.rm = TRUE) + 5,
           (work[work$layer == layer[j] & work$cFile == max(Line), ]$yRightSupS + work[work$layer == layer[j] & work$cFile == max(Line), ]$yRightInfS) / 2,
           round(dfCal[dfCal$layer == layer[j],]$RRF, 1), cex = 0.9, font = 4, col = "hotpink3")

    }
    rm(j)

    text(max(work$xRightSupS, work$xRightInfS, na.rm = TRUE) + 15,
         max(dataset$Y.cor, na.rm = TRUE) * 1.025,
         paste("RFF", dfSummary$RRF), cex = 0.9,
         font = 4, col = "hotpink3")
    
  }
  
  ##### 4. export #####
  outVcL <- list(dfCal, dfSummary)
  return(outVcL)
}

##### END OF THE FUNCTION rrf #####

#################### FUNCTION vesselID #####
# function to index each vessel
# return list of 2 objects: outIndex and ThisW

#' @dataV: dataframe of vessel sector (output of vesselS)
#' @dataC: dataframe of vessel characteristic (idSectorV, idC, cFile, Layer, nbNeigh) (output of vNeigh[1])
#' @pLab:  dataframe with all cell and their position within sector 
#' @rrf:   dataframe with number of APF per layer or radial line

vesselID <- function(dataV, dataC, pLab, rrf, make.plot = FALSE){
  outIndex <- NULL
  
  # Calculation of cumul number of cells
  ThisW               <- rrf[, c("idD", "id", "year", "DOY", "date", "layer", "RRF")]
  ThisW               <- na.omit(ThisW)
  ThisW$layer         <- as.numeric(ThisW$layer)
  ThisW               <- ThisW[order(ThisW$layer), ]
  ThisW$cumulCellLay  <- ceiling(cumsum(ThisW$RRF))
  
  ## extraction of the position of the cell within sector
  ThisP           <- pLab[pLab$idD == unique(ThisW$idD), c("idD", "idC", "type", "idSector", "PolySide", "cFile", "layer")]
  names(ThisP)[2] <- "idVessel"
  names(ThisP)[4] <- "idSectorV"
  
  ## vessel indexation
  if(any(!is.na(dataC$idC))){
    names(dataC)[2] <- "idVessel"
    
    # creation of dataframe Vcar + minLay, idSectorMin and nbCellp
    # minLay
    ThisL <- merge(dataC, aggregate(data.frame(minLay = dataV$layer), list(idVessel = dataV$idVessel), min), by = "idVessel")
    # idSectormin
    ThisL$idSectorMin <- paste("c", ThisL$cFile, ".l", ThisL$minLay, sep = "")
    # secside
    ThisL <- merge(ThisL, ThisP, by = c("idD", "idVessel", "cFile", "layer", "idSectorV"))
    # layer side of cells within idSectormin
    cell           <- dataV[which(dataV$idSector %in% ThisL$idSectorMin), c("idSector", "nbcell")]
    names(cell)[1] <- "idSectorMin"
    cell           <- unique(cell)
    cell$CellS     <- NA
    cell$CellI     <- NA
    cell$CellE     <- NA
    
    idSectorMin  <- unique(ThisL$idSectorMin)
    ThisP      <- NULL
    for(i in 1:length(idSectorMin)){
      ThisP       <- pLab[which(pLab$idD == unique(ThisW$idD)), c("idD", "idC", "type", "idSector", "PolySide", "cFile", "layer")]
      ThisP       <- ThisP[which(ThisP$idSector == idSectorMin[i]), ]
      ThisP       <- ThisP[ThisP$type != "Vessel", ]
      # addition of the "layer side" of each cell
      if(nrow(ThisP) == 0){
        cell[which(cell$idSectorMin == idSectorMin[i]), ]$CellS <- 0
        cell[which(cell$idSectorMin == idSectorMin[i]), ]$CellI <- 0
        cell[which(cell$idSectorMin == idSectorMin[i]), ]$CellE <- 0
      }else{
        ThisP$LSide <- NA
        if(length(grep("Sup", ThisP$PolySide)) != 0){
          cell[cell$idSectorMin == idSectorMin[i], ]$CellS <- 1
        }else{
            cell[cell$idSectorMin == idSectorMin[i], ]$CellS <- 0
            }
        if(length(grep("Inf", ThisP$PolySide)) != 0){
          cell[cell$idSectorMin == idSectorMin[i], ]$CellI <- 1
        }else{
            cell[cell$idSectorMin == idSectorMin[i], ]$CellI <- 0
            }
        if(length(grep("Edge", ThisP$PolySide)) != 0){
          cell[cell$idSectorMin == idSectorMin[i], ]$CellE <- 1
        }else{ cell[cell$idSectorMin == idSectorMin[i], ]$CellE <- 0
          }
      }
    }
    rm(i, ThisP)

    # addition of the "layer side" of vessel
    ThisL       <- merge(ThisL, cell, by = "idSectorMin")
    ThisL$LSide <- NA
    if(length(grep("Sup", ThisL$PolySide)) != 0){ ThisL[grep("Sup", ThisL$PolySide), ]$LSide  <- "Sup" }
    if(length(grep("Inf", ThisL$PolySide)) != 0){ ThisL[grep("Inf", ThisL$PolySide), ]$LSide  <- "Inf" }
    if(length(grep("Edge", ThisL$PolySide)) != 0){ ThisL[grep("Edge", ThisL$PolySide), ]$LSide <- "Edge" }
    
    # nbCellp within idSectormin
    ThisL$nbCellp <- 0
    
    for(v in 1:nrow(ThisL)){
      if(ThisL[v, ]$layer != ThisL[v, ]$minLay){
        ThisL[v, ]$nbCellp  <- ThisL[v, ]$CellS + ThisL[v, ]$CellI + ThisL[v, ]$CellE
      }else{
        if(ThisL[v, ]$LSide == "Sup"){
          ThisL[v, ]$nbCellp  <- ThisL[v, ]$CellI + ThisL[v, ]$CellE
        }else if(ThisL[v, ]$LSide == "Edge"){
          ThisL[v, ]$nbCellp <- ThisL[v, ]$CellI
        }
      }
    }
    rm(v)
    
    # cumulated number of cell before idSectorMin: cumulCBef
    ThisL$cumulCBef <- NA
    cell            <- NULL
    cell            <- ThisW[which(ThisW$layer %in% (ThisL$minLay - 1)), c("layer", "cumulCellLay")]
    cell$layer      <- as.numeric(cell$layer) + 1
    names(cell)[1]  <- "minLay"
    ThisL           <- merge(ThisL, cell, by = "minLay", all.x = TRUE)
    
    # rounded number of cell
    ThisL$cumulCBef    <- round(ThisL$cumulCellLay, 0)
    ThisL$cumulCellLay <- NULL
    if(nrow(ThisL[ThisL$minLay == 1,]) != 0){ ThisL[ThisL$minLay == 1,]$cumulCBef <- 0 }
    
    # vessel indexation
    ThisL$vessel.index <- ThisL$cumulCBef + ThisL$nbCellp + 1
    
    # export
    outIndex <- rbind(ThisL, outIndex)
  }else{
    outIndex <- NA
  }
  
  ##### plot for checking #####
  if(make.plot == TRUE){
    #### plot sector ####
    title <- paste(unique(pLab$idD), " vesselID", sep = "")
    layer <- sort(unique(pLab$layer))
    Line  <- sort(unique(pLab$cFile))
    
    pPopo <- plot(Y.cor ~ X.cor, pLab, main =  title, cex.main = 0.8, type = "n",
                  xlim = c(-10, max(pLab$X.cor, na.rm = TRUE) + 25),
                  ylim = c(-14, max(pLab$Y.cor, na.rm = TRUE) * 1.05))
    abline(v = seq(-10, max(pLab$X.cor, na.rm = TRUE) + 30, 10), col = "grey", lty = 3)
    abline(h = seq(0, max(pLab$Y.cor, na.rm = TRUE), 50), col = "grey", lty = 3)
    

    ## Vessels
    if(length(dataV$idSector) != 0){
      # Vessel sector
      for(j in 1:length(dataV$idSector)){
        polygon(c(dataV[j, ]$xLeftInfP,  dataV[j, ]$xLeftwork,  dataV[j, ]$xLeftSupP,
                  dataV[j, ]$xRightSupS, dataV[j, ]$xRightwork, dataV[j, ]$xRightInfS),
                c(dataV[j, ]$yLeftInfS,  dataV[j, ]$yLeftwork,  dataV[j, ]$yLeftSupS, 
                  dataV[j, ]$yRightSupS, dataV[j, ]$yRightwork, dataV[j, ]$yRightInfS),
                border = "wheat2", lwd = 2)
      }
      rm(j)
      
      # vessel points
      points(Y.cor ~ X.cor, pLab[pLab$type == "Vessel", ], col = "darkturquoise", pch = 16)
      
      # index Vessel label
      Vori <- dataV[!is.na(dataV$Vorigin),]
      for(v in 1:nrow(Vori)){
        text((Vori[v, ]$xRightInfS + Vori[v, ]$xLeftInfP + Vori[v, ]$xRightSupS + Vori[v, ]$xLeftSupP) / 4,
             (Vori[v, ]$yRightInfS + Vori[v, ]$yLeftInfS + Vori[v, ]$yRightSupS + Vori[v, ]$yLeftSupS) / 4,
             outIndex[outIndex$idVessel == Vori[v, ]$idVessel, ]$vessel.index, cex = 0.6, font = 2)
      }
      rm(v) 
    } 
    
    ## cells
    points(Y.cor ~ X.cor, pLab[pLab$type == "Fibers", ], col = "brown", pch = 16)  
    points(Y.cor ~ X.cor, pLab[pLab$type == "Axial Parenchyma", ], col = "purple2", pch = 15, cex = 0.8)
    points(Y.cor ~ X.cor, pLab[pLab$type == "APF", ], col = "maroon2", pch = 16)
    points(Y.cor ~ X.cor, pLab[pLab$type == "Cambial", ], col = "seagreen3", pch = 15)
    
    ## layer label and number of cell of the referential radial line
    layer <- sort(as.numeric(unique(pLab$layer)))
    for(j in 1:length(layer)){
      # weighted average number of APF
      text(max(pLab$X.cor, na.rm = TRUE) + 15,
           median(pLab[pLab$layer == layer[j] & pLab$cFile == max(Line), ]$Y.cor),
           round(rrf[rrf$layer == layer[j],]$RRF, 1), cex = 0.9, font = 4, col = "hotpink3")
    }
    rm(j)

    text(max(pLab$X.cor, na.rm = TRUE) + 15,
         max(pLab$Y.cor, na.rm = TRUE) * 1.025,
         round(sum(rrf$RRF, na.rm = TRUE), 0), cex = 0.9,
         font = 4, col = "hotpink3")
  }
  
  
  ### export ####
  return(list(outIndex, ThisW[which.max(ThisW$layer), ]))
}


##### END OF THE FUNCTION vesselID #####



##############################################################################################################