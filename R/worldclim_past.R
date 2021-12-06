
.worldclim_past <- function(var, res, lon, lat, past ,path,download=TRUE, downloadmethod) {
	if (!res %in% c(2.5, 5, 10)) {
		stop('resolution should be one of: 2.5, 5, 10')
	}
  
  vars <- c('tmin', 'tmax', 'prec', 'bio')
  if(!var %in% vars){
    stop('the variables available are: minimum temp (tn), maximum temp (tx), montly precipitation (pr), bioclimatic variables (bi)')
    }
	var <- c('tn', 'tx', 'pr', 'bi')[match(var, vars)]

	
  if(!past %in% c("lgm", "mid")){
    stop('the past times available are: last glacial maxima (lgm), mid-holocene (mid)')
  }
  
  if (res==2.5) { res <- '2-5' }
  
	path <- paste(path, past,var, '_', res, '/', sep='')
	message("raw data in ",path)

	dir.create(path, showWarnings=FALSE)

  	zip <- paste("cc",past,var, '_', res, 'm.zip', sep='')
  	zipfile <- paste(path, zip, sep='')
		
  	theurl <- paste0("http://biogeo.ucdavis.edu/data/climate/cmip5/",past,"/",zip)
	
		if (var  != 'bio') {
		  tiffiles <- paste0("cc",past,var,1:12,".tif")
		} else {
		  tiffiles <- paste0("cc",past,var,1:19,".tif")
		}
  
	
	files <- paste(path, tiffiles, sep='')
	
	fc <- sum(file.exists(files))
	
  print(files)
  print(fc)
	if ( fc < length(files) ) {
		if (!file.exists(zipfile)) {
			if (download) {
		  	message("trying to download from ",theurl)

				.download(theurl, zipfile,downloadmethod="wget")
				if (!file.exists(zipfile))	{ 
					message("\n Could not download file -- perhaps it does not exist") 
				}
			} else {
				message("File not available locally. Use 'download = TRUE'")
			}
		}	
		utils::unzip(zipfile, exdir=dirname(zipfile))
	}

  st <- raster::stack(paste0(path, tiffiles))
	
	projection(st) <- "+proj=longlat +datum=WGS84"
	return(st)
}
