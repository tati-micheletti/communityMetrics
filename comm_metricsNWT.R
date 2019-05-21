# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects, use sim$xxx, and are thus globally available
# to all modules. Functions can be used without sim$ as they are namespaced, like functions
# in R packages. If exact location is required, functions will be: sim$<moduleName>$FunctionName
defineModule(sim, list(
  name = "comm_metricsNWT",
  description = "This module calculates community metrics from the outputs for the bird
                model from Stralberg (unpublished)", #"insert module description here",
  keywords = c("NWT", "birds","diversity","species richness","functional diversity"), # c("insert key words here"),
  authors = c(person(c("Ana", "Raymundo", email = "angeles-ana-paula.raymundo-sanchez.1@ulaval.ca", role = c("aut", "cre"))),
              person("Steve", "Cumming", email = "stevec.boreal@gmail.com", role = c("aut"))),
  childModules = character(0),
  version = list(SpaDES.core = "0.2.3", comm_metricsNWT = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "comm_metricsNWT.Rmd"),
  reqdPkgs = list("data.table", "raster", "rgdal", "sp", "Matrix"),
  parameters = rbind(
    
    defineParameter("frequency", "numeric", 20, NA, NA, "This describes the simulation time step interval"),
    #defineParameter(".plotInitialTime", "numeric", 1, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA, "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", start(sim), NA, NA, "This describes the simulation time at which the first save event should occur"),
   #defineParameter(".saveInterval", "numeric", 10, NA, NA, "This describes the simulation interval at which the save event should occur"),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant")
  ),
  inputObjects = bind_rows(
    #expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
    expectsInput(objectName="birdPrediction", objectClass = "list", 
                 desc = "List per year of the bird species predicted rasters"),
    expectsInput(objectName = "birdsList", objectClass = "character", 
                 desc = "Bird species to be predicted", sourceURL = NA),
    expectsInput(objectName = "caribouArea1", objectClass = "SpatialPolygonsDataFrame",
                 desc = "Study area to predict caribou population to (NWT_Regions_2015_LCs_DC_SS)",
                 sourceURL = "https://drive.google.com/open?id=1Qbt2pOvC8lGg25zhfMWcc3p6q3fZtBtO"),
    expectsInput(objectName = "sp_traitsNWT", objectClass = "data.table",
                 desc = "data.table containing information about life history traits of bird species")
    ),
  outputObjects = bind_rows(
    #createsOutput("objectName", "objectClass", "output object description", ...),
    createsOutput(objectName = "currentdiversityRasters", objectClass = "RasterStack", 
                  desc = "List of Species richness, Shannon-weiner and Simpson diversity indices for the current year"),
    createsOutput(objectName = "diversityStats", objectClass = "data.frame", desc = "mean and standar deviation for each diversity indices")
    # createsOutput(objectName = "raosentropy", objectClass = "list", desc ="List per year of Raos Entropy")
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.comm_metricsNWT = function(sim, eventTime, eventType) {

  switch(eventType,
         init = {

           ## schedule future event(s)
           sim <- distanceMatrix(sim)
           sim <- Init(sim)
          
           sim <- scheduleEvent(sim, eventTime = start(sim), moduleName = "comm_metricsNWT",
                                eventType = "calcDiversity", eventPriority = .last())
           sim  <- scheduleEvent (sim, eventTime = P(sim)$.plotInitialTime, moduleName = "comm_metricsNWT", 
                                  eventType = "plot", eventPriority = .last())
           sim  <- scheduleEvent (sim, eventTime = P(sim)$.saveInitialTime, moduleName = "comm_metricsNWT", 
                                  eventType = "save", eventPriority = .last())
           sim <-scheduleEvent(sim, eventTime = start(sim), moduleName = "comm_metricsNWT", 
                               eventType = "diversityStats", eventPriority = .last())
           
         },
        calcDiversity = {
           sim <- calcDiversityIndices(sim)
           sim <- scheduleEvent(sim, eventTime = time(sim)+ P(sim)$frequency, 
                                moduleName="comm_metricsNWT", eventType = "calcDiversity", eventPriority = .last())
         },
         plot={
           sim <- Plot(sim$currentDiversityRasters)
           ## schedule future event(s)
           sim <- scheduleEvent(sim, eventTime = time(sim) + P(sim)$frequency, 
                                moduleName = "comm_metricsNWT", eventType = "plot", eventPriority = .last())
         },
         diversityStats = {
           sim <- diversityStats(sim)
           sim <- scheduleEvent(sim, eventTime = time(sim) + P (sim)$frequency, 
                                moduleName="comm_metricsNWT", eventType = "diversityStats", eventPriority = .last() )
         },
        save ={
          sim <- Save(sim)
          sim <- scheduleEvent(sim, time(sim) + P(sim)$frequency, 
                                "comm_metricsNWT", "calcDiversity", eventPriority = .last())
          
        },
         warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                       "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

Init<- function(sim){
  sim$currentDiversityRasters <- stack()
  #names(sim$currentdiversityRasters) <- c('shannon', 'simpson', 'richness', 'rao_entropy')
  sim$diversityStats <- data.frame(year = integer(), shannon_mean = double(), simpson_mean = double(), richness_mean = double(), raoentropy= double())
  sim$diversityByPolygon[[as.character(time(sim))]] <- list()
  #sim$diversityByPolygon[[paste0("Year", time(sim))]] <- list()
  return (invisible(sim))
}

diversityPlot = function(sim){
  #plot diversity
  Plot(sim$currentDiversityRasters)
  return (invisible(sim))
}   

calcDiversityIndices <- function(sim){
  birdpredsp <-(sim$birdPrediction[[length(sim$birdPrediction)]])
  bird.abun <- raster::stack(birdpredsp)
  cellSizeHA <- prod(res(bird.abun))/10000
  p <- bird.abun /sum(bird.abun)
  
  shannonRaster <- shannon(p)
  
  simpsonRaster <- simpson(p)
  
  richnessRaster <- richness(bird.abun, cellSizeHA) #second argument should be cell size
  
  SPP <- names(birdpredsp)
  dist.m <- sim$dist.m[SPP,SPP]  #requires that all species in p be in dist.m
  raoRaster <- raoEntropyC(dist.m, p)
  
  sim$currentDiversityRasters <- stack(shannonRaster,simpsonRaster,richnessRaster, raoRaster)
  names(sim$currentDiversityRasters) <- c("shannonRaster","simpsonRaster","richnessRaster", "raoRaster")
   return(invisible(sim))
}

diversityStats <- function(sim){
  sim$diversityStats <- rbind(
    sim$diversityStats,
    data.frame(indiceName = c("shannonRaster","simpsonRaster","richnessRaster", "raoRaster"),
      year = rep(time(sim), length(cellStats(sim$currentDiversityRasters, "mean"))),
      diversity_mean = cellStats(sim$currentDiversityRasters,"mean"))
  )
  
 #browser()
  sim$diversityByPolygon[[paste0("Year", as.character(time(sim)))]] <- extract(sim$currentDiversityRasters, 
                                                                               sim$caribouArea1, fun=mean, na.rm=TRUE)
  return(invisible(sim))
}

Save <- function(sim){
  dir <- file.path(outputPath(sim), "diverstiyIndices")
  dir.create(dir,showWarnings = FALSE)
  fname <- file.path(dir, paste0("diversityIndicesYear", time(sim),".tif"))
  writeRaster(sim$currentDiversityRasters, filename=fname, 
               format ="GTiff", overwrite=TRUE)
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  if (!suppliedElsewhere("traits", sim = sim)){
    sim$traits <- read.csv(file.path(dataPath(sim),"sp_traitsNWT.csv"),header=TRUE)
  }
  if(!suppliedElsewhere("birdPrediction", sim = sim)){
    sim$birdList <- c("REVI", "HETH", "RCKI", "HAFL", "WIWR", "GRCA", "RBNU", "WIWA", 
                      "GRAJ", "RBGR", "WEWP", "GCKI", "PUFI", "WETA", "FOSP", "PISI", 
                      "WCSP", "EVGR", "WBNU", "PIGR", "BTNW", "EAPH", "PHVI", "WAVI", 
                      "BRTH", "EAKI", "BRCR", "PAWA", "VESP", "DEJU", "BRBL", "OVEN", 
                      "VEER", "CSWA", "BOCH", "VATH", "OSFL", "BLPW", "COYE", "TRES", 
                      "BLJA", "OCWA", "TOWA", "TEWA", "BLBW", "CORA", "NOWA", "SWTH", 
                      "BHVI", "CONW", "MOWA", "SWSP", "BHCO", "COGR", "MAWA", "CMWA", 
                      "SOSP", "BCCH", "LISP", "YRWA", "CHSP", "SEWR", "BBWA", "LEFL", 
                      "YBFL", "CEDW", "SAVS", "BAWW", "LCSP", "WWCR", "CCSP", "RWBL", 
                      "BAOR", "HOWR", "WTSP", "CAWA", "RUBL", "AMRO", "HOLA", "AMRE", 
                      "AMGO", "AMCR", "ALFL")  
}

 if (!suppliedElsewhere(object = "studyArea", sim = sim))
  {

    sim[["studyArea"]] <- cloudCache(prepInputs,
                                    url = "https://drive.google.com/open?id=1LUxoY2-pgkCmmNH5goagBp3IMpj6YrdU",
                                    destinationPath = dataPath(sim), 
                                    cloudFolderID = sim$cloudFolderID,
                                    omitArgs = c("destinationPath", "cloudFolderID"))
    
 }
  if (!suppliedElsewhere(object = "rasterToMatch", sim = sim)){
    sim$rasterToMatch <- cloudCache(prepInputs, url = "https://drive.google.com/open?id=1fo08FMACr_aTV03lteQ7KsaoN9xGx1Df", 
                                    studyArea = sim$studyArea,
                                    targetFile = "RTM.tif", destinationPath = dataPath(sim), 
                                    useCloud = getOption("reproducible.useCloud", FALSE),
                                    cloudFolderID = sim$cloudFolderID, overwrite = TRUE, filename2 = NULL,
                                    omitArgs = c("destinationPath", "cloudFolderID", "useCloud", "overwrite", "filename2"))
  }
 if (!suppliedElsewhere("caribouArea1", sim)){
    sim$caribouArea1 <- Cache(prepInputs,
                                   url = extractURL("caribouArea1"), studyArea = sim$studyArea,
                                   destinationPath = dataPath(sim), filename2 = NULL,
                                   rasterToMatch = sim$birdPrediction)
    }      
 
 return(invisible(sim))
}