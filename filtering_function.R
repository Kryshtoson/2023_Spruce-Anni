#### 1. Fuction for distance & simialrity filtering of vegetation plots ####

###ARGUMENTS
#@coord = Coordinates of vegetation plots. A data frame with X, Y or Lon, Lat columns (should be in metres for large datasets). Plot IDs (i.e. PlotObservationID) must be in rownames.
#@spec = Species data in long format. A data frame with the following columns: 'PlotObservationID', 'Taxon_name' and 'cover' (in %). If cover is 1 for all records (representing species presence instead of cover) and sim.method is 'bray', Sorensen similarity between vegetation plots is calculated. 
#@longlat = TRUE if plot coordinates are longitude-latitude decimal degrees, in which case distances between plots are measured in kilometers. Not recommended for large datasets. Deafult is NULL.
#@dist.threshold = Plots closer than this threshold and more similar than sim.threshold will be filtered out. Units depend on coordinate system used in longlat. If longlat is NULL projected coordinate system (usually in metres) is expected. If longlat is TRUE geographical coordinate system (in degrees) is expected and longlat value must be in kilometres. Default is 1000 metres. 
#@sim.threshold = Plots with species composition more similar than this threshold and closer than dist.threshold will be filtered out. Default is 0.8.
#@sim.method = Method used for calculation of pairwise similarity between vegetation plots. Available indices are Bray-Curtis similarity ('bray'), Simpson similarity ('simpson'), Sorensen similarity ('sorensen') and Jaccard similarity ('jaccard'). If sim.method is NOT 'bray', cover values are transformed to presences.   
#@remove = Which plot should be removed, if two plots are closer than sim.threshold and more similar than sim.threshold at the same time? "random" selects and removes one random plot, "less diverse" selects and removes less diverse plot. Default  is 'random'.
#@seed = Seed for repeatability of the result. Default is 1234.
#@parallelize = Should similarity between plots be calculated using parallel computation? If TRUE, 'parallel' and 'pbapply' packages are installed (if not present) and calculations are run on available cores. Parallel computation may speed up filteirng process but it is used only with No. of available cores > 3. Default is FALSE. 

###USAGE
#filtering(coord, spec, longlat=NULL, dist.threshold=1000, sim.threshold=0.8, sim.method=c("bray", "simpson", "sorensen", "jaccard"), remove=c("random", "less diverse"), seed=1234, parallelize=FALSE)

###CODES
similarity <- function(dd, coord, spec, sim.threshold, sim.method)
{
  if(!any(dd == 0))
  {
    spe.sel <- spec[spec$PlotObservationID %in% rownames(coord)[dd], ]  
    tab <- reshape::cast(spe.sel, PlotObservationID ~ Taxon_name, value="cover", fun.aggregate = length, fill = 0)
    rownames(tab) <- tab$PlotObservationID
    tab <- tab[rownames(coord)[dd], ]
    
    if(sim.method[1] == "bray"){s <- c(0, 1-vegan::vegdist(tab[,-1], method="bray")[1:(length(dd)-1)])} ##calculate similarity Bray-Curtis
    if(sim.method[1] == "simpson"){s <- c(0, 1-vegan::betadiver(tab[,-1], method="sim")[1:(length(dd)-1)])} ##calculate similarity Simpson
    if(sim.method[1] == "sorensen"){s <- c(0, vegan::betadiver(tab[,-1], method="sor")[1:(length(dd)-1)])} ##calculate similarity Sorensen
    if(sim.method[1] == "jaccard"){s <- c(0, vegan::betadiver(tab[,-1], method="j")[1:(length(dd)-1)])} ##calculate similarity Jaccard
    
    s[s <= sim.threshold] <- 0 #apply similarity threshold - these plots will NOT be filtered
    
    return(s)
    
  }
  else
  {
    return(c(0, 0)) #for similarity
  }
  
}

filtering <- function(coord, spec, longlat=NULL, dist.threshold=1000, sim.threshold=0.8, sim.method=c("bray", "simpson", "sorensen", "jaccard"), 
                      remove=c("random", "less diverse"), seed=1234, parallelize=FALSE)
{
  #load required packages
  if (!require("spdep")) install.packages("spdep")
  if (!require("vegan")) install.packages("vegan")
  if (!require("reshape")) install.packages("reshape")
  if(parallelize){
    if (!require("parallel")) install.packages("parallel")
    if (!require("pbapply")) install.packages("pbapply")}
  
  #order releves based on No. species or make random order
  if(remove[1] == "random")
  {
    set.seed(seed)
    coord <- coord[sample(nrow(coord)), ]
  }
  if(remove[1] == "less diverse")
  {
    div <- sort(table(spec$PlotObservationID))
    coord <- coord[names(div), ]
  }
  
  ##transform covers to p/a data if desired
  if(sim.method[1] != "bray"){spec$cover <- 1}
  
  #calculate which releves are within dist.threshold
  cat("Searching for neighbours... ...please wait", "\n")
  d <- dnearneigh(as.matrix(coord), d1 = 0, d2 = dist.threshold, bounds = c("GE", "LE"), longlat=longlat)
  
  ##add releve IDs to each element of d
  d <- mapply(append, seq_along(d), d, SIMPLIFY = FALSE)
  
  ##calculate similarity among neighbouring releves************************************************
  if(parallelize)
  {
    no.cores <- detectCores() #check No. cores
    if(no.cores >= 3) {no.cores <- no.cores-1}
    if(no.cores < 3) {no.cores <- 1}
    
    if(no.cores == 1)
    {
      cat("Parallelization unavailable... ...using only one core", "\n")
      cat("Calculating similarity between releves... ...please wait", "\n")
      ds <- lapply(d, similarity, coord, spec, sim.threshold, sim.method)
    }
    if(no.cores > 1)
    {
      cat(paste0("Parallel calculation using ", no.cores, " cores..."), "\n")
      cat("Calculating similarity between releves... ...please wait", "\n")
      cl <- makeCluster(no.cores)
      #clusterExport(cl, list("cast", "vegdist", "betadiver"))
      
      ds <- pblapply(d, FUN = similarity, cl = cl, coord, spec, sim.threshold, sim.method)
      stopCluster(cl)
    }
  }
  else
  {
    cat("Calculating similarity between releves... ...please wait", "\n")
    ds <- lapply(d, similarity, coord, spec, sim.threshold, sim.method)
  }
  
  ###FILTER RELEVES****************************************************************
  
  cat("Removing releves... ...please wait","\n")
  plotID.to.remove <- NULL
  repeat
  {
    id1 <- which.max(unlist(lapply(ds, FUN=max))) #identify 1st releve from the most similar pair
    id2 <- d[[id1]][which(ds[[id1]] > 0)] #identify 2nd releve from the most similar pair
    
    plotID.to.remove <- c(plotID.to.remove, rownames(coord)[id1]) #record releve IDs to remove
    
    ds[[id1]] <- rep(0, length(ds[[id1]])) #assign 0 similarity to this pair
    
    for(i in id2)
    {
      ds[[i]][d[[i]] == id1] <- 0 #assign 0 similarity to this pair
    }
    
    if(all(unlist(lapply(ds, FUN=max)) == 0))
    {
      break
    }
  }
  
  coord.filtered <- coord[setdiff(rownames(coord), plotID.to.remove), ]#selected subset
  
}

###VALUE
#The function returns a data frame with coordinates of preserved plots.

###AUTHOR
#Jan Divíšek (divisekjan@sci.muni.cz)

###EXAMPLE - requires your own data;-)
# 
# res <- filtering(coord=plot_coordinates, spec=species_data, longlat=NULL, dist.threshold=1000, sim.threshold=0.8, sim.method="simpson", 
#                  remove="random", seed=1234, parallelize=TRUE)
# head(res)
# 
# plot(plot_coordinates)
# points(res, pch=16, cex=0.5, col='red')
