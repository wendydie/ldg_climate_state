source("./R/functions/calculate_Info.R")

# return vector of cells that lie within buffer radius of given seed
findPool <- function(seed, dat, siteId, xy, r, nSite, crs = 'epsg:4326'
){
  datSf <- sf::st_as_sf(dat, coords = xy, crs = crs)
  seedRow <- which(dat[, siteId] == seed)[1]
  seedpt <- datSf[seedRow, ]
  # buffer will be more accurate if projected,
  # but wrapping around antimeridian requires lat-long coordinates
  r <- units::set_units(r, 'km')
  buf <- sf::st_buffer(seedpt, dist = r)
  if (crs != 'epsg:4326'){
    buf   <- sf::st_transform(buf,   crs = 'epsg:4326')
    datSf <- sf::st_transform(datSf, crs = 'epsg:4326')
  }
  bufWrap <- sf::st_wrap_dateline(buf, options = c("WRAPDATELINE=YES"))
  
  # find sites within radius of seed site/cell
  poolBool <- sf::st_intersects(datSf, bufWrap, sparse = FALSE)
  pool <- dat[poolBool, siteId]
  return(pool)
}

# function to try all possible starting pts (i.e. all occupied cells)
# save the ID of any cells that contain given pool size within buffer
findSeeds <- function(dat, siteId, xy, r, nSite, crs = 'epsg:4326'
){
  # test whether each occupied site/cell is viable for subsampling
  posSeeds <- dat[,siteId]
  # don't use sapply or object will condense from list to matrix
  # in the special case all pool lengths are equal:
  posPools <- lapply(posSeeds, function(s){
    sPool <- findPool(s, dat, siteId, xy, r, nSite, crs)
    n <- length(sPool)
    if (n >= nSite)
      sPool
  })
  # return pool site/cell IDs for each viable seed point
  # same overall list structure as cookies outputs; names = seed IDs
  names(posPools) <- posSeeds
  Filter(Negate(is.null), posPools)
}


# Buffers function to perform spatial subsampling and optionally return either incidence frequency or full data
#' Perform Spatial Subsampling with Buffers
#'
#' This function generates spatial buffers around specified seed points in a dataset
#' and performs subsampling within these buffers. It can return either an 
#' incidence frequency matrix or the full subsampled data, depending on the 
#' `output` parameter.
#'
#' @param dat A data frame containing input data, including geographic 
#'        coordinates (latitude and longitude) and other associated data (e.g., species or site information).
#'        Example:
#'        dat <- data.frame(
#'          id = c(1, 2, 3),
#'          lat = c(34.05, 36.77, 40.71),
#'          lng = c(-118.25, -119.42, -74.01),
#'          genus = c("Species_A", "Species_B", "Species_C")
#'        )
#'
#' @param xy A character vector specifying the column names in `dat` that 
#'        represent the geographic coordinates (latitude and longitude).
#'        Example:
#'        xy <- c("lat", "lng")
#'
#' @param nSite Integer specifying the number of sites (points) to sample 
#'        within each buffer. This determines how many points are included 
#'        in each buffer subsample.
#'        Example:
#'        nSite <- 10
#'
#' @param r Numeric value specifying the radius (in kilometers) for the 
#'        buffer around each seed point. This defines the spatial extent of 
#'        the buffer.
#'        Example:
#'        r <- 100
#'
#' @param crs A string specifying the coordinate reference system for the 
#'        spatial data. Default is `'epsg:4326'` (WGS84).
#'
#' @param output A string specifying the desired output format. 
#'        Options:
#'        - `'incidence_freq'`: Returns an incidence frequency matrix.
#'        - `'full'`: Returns full data for points within each buffer.
#'        Default: `'locs'`
#'
#' @return A list of length `iter`. Each list element contains either:
#'         - A matrix of incidence frequencies (`output = 'incidence_freq'`), or
#'         - A combined data frame of buffered points (`output = 'full'`).
#'
#' @examples
#' dat <- data.frame(
#'   id = c(1, 2, 3),
#'   lat = c(34.05, 36.77, 40.71),
#'   lng = c(-118.25, -119.42, -74.01),
#'   genus = c("Species_A", "Species_B", "Species_C")
#' )
#' xy <- c("lat", "lng")
#' r <- 200  # Buffer radius in km
#' nSite <- 5
#' iter <- 3
#' buffers(dat, xy, nSite, r, output = 'incidence_freq')
#'
buffers <- function(dat, xy, nSite, r,
                    crs = 'epsg:4326', output = 'incidence_freq') {
  
  # Ensure unique coordinates for processing and add IDs for each site
  coords <- uniqify(dat, xy) |> as.data.frame()
  coords$id <- paste0('loc', 1:nrow(coords))
  
  # Generate spatial buffers for subsampling; this is computationally intensive
  allPools <- findSeeds(coords, 'id', xy, r, nSite, crs)
  if (length(allPools) < 1) {
    # stop('Not enough close sites for any subsample.')
    warning('Not enough close sites for any subsample.')
    return(NULL)  # Skip to the next iteration of the loop
  }
  
  seeds <- names(allPools)
  if (output == 'incidence_freq') {
    # Generate incidence frequency matrix from buffered data
    incfreq_inbuffers <- lapply(seeds, function(seed) {
      # Retrieve buffered points
      pool <- allPools[seed][[1]]
      # samplIds <- sample(sample(pool), nSite, replace = FALSE)
      samplIds <- pool
      
      # Match sampled points to original data
      coordRows <- match(samplIds, coords$id)
      coordLocs <- coords[coordRows, xy]
      x <- xy[1]
      y <- xy[2]
      sampPtStrg <- paste(coordLocs[, x], coordLocs[, y], sep = '/')
      datPtStrg <- paste(dat[, x], dat[, y], sep = '/')
      inSamp <- match(datPtStrg, sampPtStrg)
      dat_inbuffer <- dat[!is.na(inSamp), ]
      
      # Extract seed coordinates
      seed_coord <- coords[coords$id == seed, xy]
      seed_lat <- seed_coord[[xy[1]]]
      seed_lng <- seed_coord[[xy[2]]]
      
      # Generate incidence frequency
      incfreq_inbuffer <- incfreq(dat_inbuffer)
      seed_name <- paste0(seed_lat, "_", seed_lng)
      names(incfreq_inbuffer) <- seed_name
      return(incfreq_inbuffer)
    })
    
    # Combine all incidence frequencies into a single matrix
    out <- do.call(c, incfreq_inbuffers)
    
  } else if (output == 'full') {
    # Return full data for each buffer
    dat_inbuffers <- lapply(seeds, function(seed) {
      # Retrieve buffered points
      pool <- allPools[seed][[1]]
      # samplIds <- sample(sample(pool), nSite, replace = FALSE)
      samplIds <- pool
      
      # Match sampled points to original data
      coordRows <- match(samplIds, coords$id)
      coordLocs <- coords[coordRows, xy]
      x <- xy[1]
      y <- xy[2]
      sampPtStrg <- paste(coordLocs[, x], coordLocs[, y], sep = '/')
      datPtStrg <- paste(dat[, x], dat[, y], sep = '/')
      inSamp <- match(datPtStrg, sampPtStrg)
      dat_inbuffer <- dat[!is.na(inSamp), ]
      
      # Add buffer metadata
      seed_coord <- coords[coords$id == seed, xy]
      seed_lat <- seed_coord[[xy[1]]]
      seed_lng <- seed_coord[[xy[2]]]
      buffer_no <- paste0(seed_lat, "_", seed_lng)
      dat_inbuffer['buffer_no'] <- buffer_no
      
      return(dat_inbuffer)
    })
    
    # Combine all buffered data into a single data frame
    out <- do.call(rbind, dat_inbuffers)
  } else if (output == 'locs'){
    locs_inbuffers <- lapply(seeds, function(seed) {
      # Retrieve buffered points
      pool <- allPools[seed][[1]]
      # samplIds <- sample(sample(pool), nSite, replace = FALSE)
      samplIds <- pool
      
      # Match sampled points to original data
      coordRows <- match(samplIds, coords$id)
      coordLocs <- coords[coordRows, xy]
      x <- xy[1]
      y <- xy[2]
      # Add buffer metadata
      seed_coord <- coords[coords$id == seed, xy]
      seed_lat <- seed_coord[[xy[1]]]
      seed_lng <- seed_coord[[xy[2]]]
      buffer_no <- paste0(seed_lat, "_", seed_lng)
      coordLocs['buffer_no'] <- buffer_no
      return(coordLocs)
    })
    # Combine all buffered data into a single data frame
    out <- do.call(rbind, locs_inbuffers)
  }
  
  return(out)
}
