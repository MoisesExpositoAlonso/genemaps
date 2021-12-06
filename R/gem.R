


# Download recursively geographic data and return as R object
# Author: Moi Exposito-Alonso
# January 2017
#' @export

recursive.getXData <- function(times=c("pres","MP8570"),
                              var="bio" ,res="2.5",path="~/"){
require(raster)
designtimes <- list(
    lgm= "worldclim_past",
    mid = "worldclim_past",
    pres= "worldclim",
    MP2650="CMIP5",
    MP2670="CMIP5",
    MP4550="CMIP5",
    MP4570="CMIP5",
    MP6050="CMIP5",
    MP6070="CMIP5",
    MP8550="CMIP5",
    MP8570="CMIP5")


if( any(times == "all" ) ) {
  times <-  names(designtimes)
  }

if( any(!times %in%  names(designtimes) ) ) {
  message("These are the available dataset for recursie download")
  print(names(designtimes))
}

sapply( times ,function(dataset){
  map.getXData(dataset,thefun=as.character(designtimes[dataset]),var=var,res=res,path=path)
  }
)

}

map.getXData <- function(dataset,thefun, var,res,path){
message("importing dataset ", dataset)
  if(grepl("[A-Z][A-Z][0-9][0-9][0-9][0-9]",dataset)){
    model=substr(dataset, 1, 2)
    rcp=substr(dataset, 3, 4)
    year=substr(dataset, 5, 6)

    tmp<- getXData(name =thefun,model=model,rcp=rcp,year=year,
                  download = T,var=var,res=res,path =path)
  }
  else if(dataset=="pres"){
    tmp<- getXData(name =thefun,
            download = T,var=var,res=res,path =path)
  }
  else{
    tmp<- getXData(name =thefun,past= dataset,
            download = T,var=var,res=res,path =path)
  }

names(tmp) <- paste0(var,c(1:nlayers(tmp) ))
return(tmp)
}

####************************************************************************####

#' getXData exteded Function
#'
#' Download geographic data and return as R object
#  Modified by Moi Exposito-Alonso from Raster package
#
#' @param name
#' @param download
#' @param path
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
getXData <- function(name='GADM', download=TRUE, path='', ...) {
require(raster)

	path <- raster:::.getDataPath(path)
	if (name=='GADM') {
		.GADM(..., download=download, path=path, version=2.8)
	} else if (name=='SRTM') {
		.SRTM(..., download=download, path=path)
	} else if (name=='alt') {
		.raster(..., name=name, download=download, path=path)
	} else if (name=='worldclim') {
		.worldclim(..., download=download, path=path)
	} else if (name=='worldclim_past') {
		.worldclim_past(..., download=download, path=path)
	} else if (name=='CMIP5') {
		.cmip5(..., download=download, path=path)
	} else if (name=='ISO3') {
		ccodes()[,c(2,1)]
	} else if (name=='countries') {
		.countries(download=download, path=path, ...)
	} else {
		stop(name, ' not recognized as a valid name.')
	}
}


.download <- function(aurl, filename, downloadmethod="auto") {

	fn <- paste(tempfile(), '.download', sep='')
	# res <- utils::download.file(url=aurl, destfile=fn, method="auto", quiet = FALSE, mode = "wb", cacheOK = TRUE)
	# res <- utils::download.file(url=aurl, destfile=fn, method="curl", quiet = FALSE, mode = "wb", cacheOK = TRUE) # other people and me problem with method="auto"
	# res <- utils::download.file(url=aurl, destfile=fn, method="wget", quiet = FALSE, mode = "wb", cacheOK = TRUE) # this works in my linux machine
	res <- utils::download.file(url=aurl, destfile=fn, method=downloadmethod, quiet = FALSE, mode = "wb", cacheOK = TRUE) # other people and me problem with method="auto"
	if (res == 0) {
		w <- getOption('warn')
		on.exit(options('warn' = w))
		options('warn'=-1)
		if (! file.rename(fn, filename) ) {
			# rename failed, perhaps because fn and filename refer to different devices
			file.copy(fn, filename)
			file.remove(fn)
		}
	} else {
		stop('could not download the file' )
	}
}

.ISO <- function() {
   ccodes()
}

ccodes <- function() {
	path <- system.file(package="raster")
	#d <- utils::read.csv(paste(path, "/external/countries.csv", sep=""), stringsAsFactors=FALSE, encoding="UTF-8")
	readRDS(file.path(path, "external/countries.rds"))
}


.getCountry <- function(country='') {
	country <- toupper(trim(country[1]))

	cs <- ccodes()
	cs <- sapply(cs, toupper)
	cs <- data.frame(cs, stringsAsFactors=FALSE)
	nc <- nchar(country)

	if (nc == 3) {
		if (country %in% cs$ISO3) {
			return(country)
		} else {
			stop('unknown country')
		}
	} else if (nc == 2) {
		if (country %in% cs$ISO2) {
			i <- which(country==cs$ISO2)
			return( cs$ISO3[i] )
		} else {
			stop('unknown country')
		}
	} else if (country %in% cs[,1]) {
		i <- which(country==cs[,1])
		return( cs$ISO3[i] )
	} else if (country %in% cs[,4]) {
		i <- which(country==cs[,4])
		return( cs$ISO3[i] )
	} else if (country %in% cs[,5]) {
		i <- which(country==cs[,5])
		return( cs$ISO3[i] )
	} else {
		stop('provide a valid name name or 3 letter ISO country code; you can get a list with "ccodes()"')
	}
}


.getDataPath <- function(path) {
	path <- trim(path)
	if (path=='') {
		path <- .dataloc()
	} else {
		if (substr(path, .nchar(path)-1, .nchar(path)) == '//' ) {
			p <- substr(path, 1, .nchar(path)-2)
		} else if (substr(path, .nchar(path), .nchar(path)) == '/'  | substr(path, .nchar(path), .nchar(path)) == '\\') {
			p <- substr(path, 1, .nchar(path)-1)
		} else {
			p <- path
		}
		if (!file.exists(p) & !file.exists(path)) {
			stop('path does not exist: ', path)
		}
	}
	if (substr(path, .nchar(path), .nchar(path)) != '/' & substr(path, .nchar(path), .nchar(path)) != '\\') {
		path <- paste(path, "/", sep="")
	}
	return(path)
}


.GADM <- function(country, level, download, path, version) {
#	if (!file.exists(path)) {  dir.create(path, recursive=T)  }

	country <- .getCountry(country)
	if (missing(level)) {
		stop('provide a "level=" argument; levels can be 0, 1, or 2 for most countries, and higher for some')
	}

	filename <- paste(path, 'GADM_', version, '_', country, '_adm', level, ".rds", sep="")
	if (!file.exists(filename)) {
		if (download) {

			baseurl <- paste0("http://biogeo.ucdavis.edu/data/gadm", version)
			if (version == 2) {
				theurl <- paste(baseurl, '/R/', country, '_adm', level, ".RData", sep="")
			} else {
				theurl <- paste(baseurl, '/rds/', country, '_adm', level, ".rds", sep="")
			}

			.download(theurl, filename)
			if (!file.exists(filename))	{
				message("\nCould not download file -- perhaps it does not exist")
			}
		} else {
			message("File not available locally. Use 'download = TRUE'")
		}
	}
	if (file.exists(filename)) {
		if (version == 2) {
			thisenvir <- new.env()
			data <- get(load(filename, thisenvir), thisenvir)
		} else {
			data <- readRDS(filename)
		}
		return(data)
	} else {
		return(NULL)
	}
}


.countries <- function(download, path, ...) {
#	if (!file.exists(path)) {  dir.create(path, recursive=T)  }
	filename <- paste(path, 'countries.RData', sep="")
	if (!file.exists(filename)) {
		if (download) {
			theurl <- paste("http://biogeo.ucdavis.edu/data/gadm2.6/countries_gadm26.rds", sep="")
			.download(theurl, filename)
			if (!file.exists(filename)) {
				message("\nCould not download file -- perhaps it does not exist")
			}
		} else {
			message("File not available locally. Use 'download = TRUE'")
		}
	}
	if (file.exists(filename)) {
		#thisenvir = new.env()
		#data <- get(load(filename, thisenvir), thisenvir)
		data <- readRDS(filename)
		return(data)
	}
}


.cmip5 <- function(var, model, rcp, year, res, lon, lat, path, download=TRUE) {
	if (!res %in% c(0.5, 2.5, 5, 10)) {
		stop('resolution should be one of: 2.5, 5, 10')
	}
	if (res==2.5) {
		res <- '2_5m'
    } else if (res == 0.5) {
        res <- "30s"
    } else {
		res <- paste(res, 'm', sep='')
	}

	var <- tolower(var[1])
	vars <- c('tmin', 'tmax', 'prec', 'bio')
	stopifnot(var %in% vars)
	var <- c('tn', 'tx', 'pr', 'bi')[match(var, vars)]

	model <- toupper(model)
	models <- c('AC', 'BC', 'CC', 'CE', 'CN', 'GF', 'GD', 'GS', 'HD', 'HG', 'HE', 'IN', 'IP', 'MI', 'MR', 'MC', 'MP', 'MG', 'NO')
	stopifnot(model %in% models)

	rcps <- c(26, 45, 60, 85)
	stopifnot(rcp %in% rcps)
	stopifnot(year %in% c(50, 70))

	m <- matrix(c(0,1,1,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,1,0,1,1,1,0,0,1,1,1,1,0,1,1,1,1,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1), ncol=4)
	i <- m[which(model==models), which(rcp==rcps)]
	if (!i) {
		warning('this combination of rcp and model is not available')
		return(invisible(NULL))
	}

	path <- paste(path, '/cmip5/', res, '/', sep='')
	dir.create(path, recursive=TRUE, showWarnings=FALSE)

	zip <- tolower(paste(model, rcp, var, year, '.zip', sep=''))
	theurl <- paste('http://biogeo.ucdavis.edu/data/climate/cmip5/', res, '/', zip, sep='')

	zipfile <- paste(path, zip, sep='')
	if (var == 'bi') {
		n <- 19
	} else {
		n <- 12
	}
	tifs <- paste(extension(zip, ''), 1:n, '.tif', sep='')
	files <- paste(path, tifs, sep='')
	fc <- sum(file.exists(files))
	if (fc < n) {
		if (!file.exists(zipfile)) {
			if (download) {
				.download(theurl, zipfile)
				if (!file.exists(zipfile))	{
					message("\n Could not download file -- perhaps it does not exist")
				}
			} else {
				message("File not available locally. Use 'download = TRUE'")
			}
		}
		utils::unzip(zipfile, exdir=dirname(zipfile))
	}
	stack(paste(path, tifs, sep=''))
}

#.cmip5(var='prec', model='BC', rcp=26, year=50, res=10, path=getwd())

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

	st <- stack(paste0(path, tiffiles))

	projection(st) <- "+proj=longlat +datum=WGS84"
	return(st)
}


.worldclim <- function(var, res, lon, lat, path, download=TRUE) {
	if (!res %in% c(0.5, 2.5, 5, 10)) {
		stop('resolution should be one of: 0.5, 2.5, 5, 10')
	}
	if (res==2.5) { res <- '2-5' }

	stopifnot(var %in% c('tmean', 'tmin', 'tmax', 'prec', 'bio', 'alt'))
	path <- paste(path, 'wc', res, '/', sep='')
	dir.create(path, showWarnings=FALSE)

	if (res==0.5) {
		lon <- min(180, max(-180, lon))
		lat <- min(90, max(-60, lat))
		rs <- raster(nrows=5, ncols=12, xmn=-180, xmx=180, ymn=-60, ymx=90 )
		row <- rowFromY(rs, lat) - 1
		col <- colFromX(rs, lon) - 1
		rc <- paste(row, col, sep='')
		zip <- paste(var, '_', rc, '.zip', sep='')
		zipfile <- paste(path, zip, sep='')
		if (var  == 'alt') {
			bilfiles <- paste(var, '_', rc, '.bil', sep='')
			hdrfiles <- paste(var, '_', rc, '.hdr', sep='')
		} else if (var  != 'bio') {
			bilfiles <- paste(var, 1:12, '_', rc, '.bil', sep='')
			hdrfiles <- paste(var, 1:12, '_', rc, '.hdr', sep='')
		} else {
			bilfiles <- paste(var, 1:19, '_', rc, '.bil', sep='')
			hdrfiles <- paste(var, 1:19, '_', rc, '.hdr', sep='')
		}
		theurl <- paste('http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/tiles/cur/', zip, sep='')
	} else {
		zip <- paste(var, '_', res, 'm_bil.zip', sep='')
		zipfile <- paste(path, zip, sep='')
		if (var  == 'alt') {
			bilfiles <- paste(var, '.bil', sep='')
			hdrfiles <- paste(var, '.hdr', sep='')
		} else if (var  != 'bio') {
			bilfiles <- paste(var, 1:12, '.bil', sep='')
			hdrfiles <- paste(var, 1:12, '.hdr', sep='')
		} else {
			bilfiles <- paste(var, 1:19, '.bil', sep='')
			hdrfiles <- paste(var, 1:19, '.hdr', sep='')
		}
		theurl <- paste('http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/', zip, sep='')
	}
	files <- c(paste(path, bilfiles, sep=''), paste(path, hdrfiles, sep=''))
	fc <- sum(file.exists(files))


	if ( fc < length(files) ) {
		if (!file.exists(zipfile)) {
			if (download) {
				.download(theurl, zipfile)
				if (!file.exists(zipfile))	{
					message("\n Could not download file -- perhaps it does not exist")
				}
			} else {
				message("File not available locally. Use 'download = TRUE'")
			}
		}
		utils::unzip(zipfile, exdir=dirname(zipfile))
		for (h in paste(path, hdrfiles, sep='')) {
			x <- readLines(h)
			x <- c(x[1:14], 'PIXELTYPE     SIGNEDINT', x[15:length(x)])
			writeLines(x, h)
		}
	}
	if (var  == 'alt') {
		st <- raster(paste(path, bilfiles, sep=''))
	} else {
		st <- stack(paste(path, bilfiles, sep=''))
	}
	projection(st) <- "+proj=longlat +datum=WGS84"
	return(st)
}


.raster <- function(country, name, mask=TRUE, path, download, keepzip=FALSE, ...) {

	country <- .getCountry(country)
	path <- .getDataPath(path)
	if (mask) {
		mskname <- '_msk_'
		mskpath <- 'msk_'
	} else {
		mskname<-'_'
		mskpath <- ''
	}
	filename <- paste(path, country, mskname, name, ".grd", sep="")
	if (!file.exists(filename)) {
		zipfilename <- filename
		extension(zipfilename) <- '.zip'
		if (!file.exists(zipfilename)) {
			if (download) {
				theurl <- paste("http://biogeo.ucdavis.edu/data/diva/", mskpath, name, "/", country, mskname, name, ".zip", sep="")
				.download(theurl, zipfilename)
				if (!file.exists(zipfilename))	{
					message("\nCould not download file -- perhaps it does not exist")
				}
			} else {
				message("File not available locally. Use 'download = TRUE'")
			}
		}
		ff <- utils::unzip(zipfilename, exdir=dirname(zipfilename))
		if (!keepzip) {
			file.remove(zipfilename)
		}
	}
	if (file.exists(filename)) {
		rs <- raster(filename)
	} else {
		#patrn <- paste(country, '.', mskname, name, ".grd", sep="")
		#f <- list.files(path, pattern=patrn)
		f <- ff[substr(ff, .nchar(ff)-3, .nchar(ff)) == '.grd']
		if (length(f)==0) {
			warning('something went wrong')
			return(NULL)
		} else if (length(f)==1) {
			rs <- raster(f)
		} else {
			rs <- sapply(f, raster)
			message('returning a list of RasterLayer objects')
			return(rs)
		}
	}
	projection(rs) <- "+proj=longlat +datum=WGS84"
	return(rs)
}



.SRTM <- function(lon, lat, download, path) {
	stopifnot(lon >= -180 & lon <= 180)
	stopifnot(lat >= -60 & lat <= 60)

	rs <- raster(nrows=24, ncols=72, xmn=-180, xmx=180, ymn=-60, ymx=60 )
	rowTile <- rowFromY(rs, lat)
	colTile <- colFromX(rs, lon)
	if (rowTile < 10) { rowTile <- paste('0', rowTile, sep='') }
	if (colTile < 10) { colTile <- paste('0', colTile, sep='') }

	f <- paste('srtm_', colTile, '_', rowTile, sep="")
	zipfilename <- paste(path, "/", f, ".ZIP", sep="")
	tiffilename <- paste(path, "/", f, ".TIF", sep="")

	if (!file.exists(tiffilename)) {
		if (!file.exists(zipfilename)) {
			if (download) {
				theurl <- paste("ftp://xftp.jrc.it/pub/srtmV4/tiff/", f, ".zip", sep="")
				test <- try (.download(theurl, zipfilename) , silent=TRUE)
				if (class(test) == 'try-error') {
					theurl <- paste("http://hypersphere.telascience.org/elevation/cgiar_srtm_v4/tiff/zip/", f, ".ZIP", sep="")
					test <- try (.download(theurl, zipfilename) , silent=TRUE)
					if (class(test) == 'try-error') {
						theurl <- paste("http://srtm.csi.cgiar.org/SRT-ZIP/SRTM_V41/SRTM_Data_GeoTiff/", f, ".ZIP", sep="")
						.download(theurl, zipfilename)
					}
				}
			} else {message('file not available locally, use download=TRUE') }
		}
		if (file.exists(zipfilename)) {
			utils::unzip(zipfilename, exdir=dirname(zipfilename))
			file.remove(zipfilename)
		}
	}
	if (file.exists(tiffilename)) {
		rs <- raster(tiffilename)
		projection(rs) <- "+proj=longlat +datum=WGS84"
		return(rs)
	} else {
		stop('file not found')
	}
}



####************************************************************************####


# make_Euroclim<-function(code="pres",refresh=FALSE){
# require(raster)
#   stopifnot(code %in% c('pres',
#                            'MP2650',
#                            'MP2670',
#                            'MP4550',
#                            'MP4570',
#                            'MP6050',
#                            'MP6070',
#                            'MP8550',
#                            'MP8570'
#                            ) )
#
#   require(raster)
#
#   if(code=='pres'){
#     ras=raster::getData(name ="worldclim" ,res="2.5",path="~/", var="bio",download = refresh)
#   }else{
#     model=substr(code, 1, 2)
#     rcp=substr(code, 3, 4)
#     year=substr(code, 5, 6)
#     ras=raster::getData(name ="CMIP5" ,res="2.5",path="~/",var="bio",model=model,rcp=rcp,year=year,download = refresh) # MP for Max Planck
#   }
#
#   Euroclim<-stack(sapply(1:19,function(i) cropenvironment(ras[[i]]) ))
#
#   return(Euroclim)
# }
#

make_Euroclim<-function(clim="pres"){
  if(clim=='pres'){clim<-''}
  clim<-stack(paste0('dataint/',clim,'newclim2.grd'))

return(clim)
}

load_bio<-function(refresh=FALSE){
    pres=raster::getData(name ="worldclim" ,res="2.5",path="~/", var="bio",download = refresh)
}



#
# overlapbio<-function(lon, lat, refresh=FALSE){
#
#   if( ! "pres" %in% ls()){
#     pres=load_bio(refresh)
#   }
# lon<-fn(lon)
# lat<-fn(lat)
# locs<-data.frame(lon,lat)
#
# res<- sapply(1:19,function(i) raster::extract( pres[[i]], locs))
#
# colnames(res) <- paste0('bio',1:19)
#
# return(res)
#
# }

overlapbio<-function(lon, lat, rast=NULL, refresh=FALSE){
  if(is.null(rast)){
    if( ! "pres" %in% ls()){
      pres=load_bio(refresh)
    }
  }else{
    pres=rast
  }

  lon<-fn(lon)
  lat<-fn(lat)
  locs<-data.frame(lon,lat)

  res<- sapply(1:length(names(pres)),function(i) raster::extract( pres[[i]], locs))
  if(class(res) == "numeric"){
    names(res) <- names(pres)
  }else{
    colnames(res) <- names(pres)
  }
  return(res)
}

make_fieldenvs<-function(){

  if( ! "Euroclim" %in% ls()){
  Euroclim<-make_Euroclim()
  }

  locs<-fieldlocations() %>% t() %>% as.data.frame %>% dplyr::select(lon, lat)

  fieldenvs<- sapply(1:length(names(Euroclim)),function(i) raster::extract( Euroclim[[i]], locs))
  colnames(fieldenvs) <- names(Euroclim)
  rownames(fieldenvs)<-rownames(locs)
  fieldenvs
  return(fieldenvs)

}


#####**********************************************************************#####
#### Selection load ####


make_metrics_dummyraster<-function(ras=refraster,metrics=allmetrics[onesnp,allmetcols]){
  mycols<-colnames(metrics)
  metricsmap<-lapply(allmetcols, function(i){
    temp<-ras
    names(temp)<-i
    temp<-temp*metrics[,i]
    return(temp)
  }) %>% raster::stack()
  return(metricsmap)
}

#####**********************************************************************#####
#### Utilities ####

makebaseraster<-function(refraster=NULL){
  baseraster<- refraster
  if( is.null(baseraster) | class(baseraster) !='RasterLayer'  ){
     baseraster<-readRDS("dataint/baseraster.rda")
  }
  baseraster[!is.na(baseraster)] <-1
  return(baseraster)
}

plotbaseraster<-function(baseraster=NULL, reference=TRUE, color='gray80'){
  if( is.null(baseraster) | class(baseraster) !='RasterLayer'  ){
     baseraster<-readRDS("dataint/baseraster.rda")
  }
  moiR::envirplot(baseraster,vecol=color,discrete=T,plotlegend=F)
  addreferencebar()
}

addreferencebar<-function(startx= 15, starty=35, size=500 ){
  kmbydegree<-88.42545 # calculated at 37 degrees latitude
  endx=startx +(size/kmbydegree)
  lines(x = c(startx,endx ) ,y= c(starty,starty) ,col='darkgrey')
  # text(x=startx, y=starty-1,'0',col='darkgrey')
  # text(x=endx, y=starty-1,paste("    ",size, 'Km'),col='darkgrey')
  text(x=endx-(endx-startx)/2, y=starty-1,paste(size, 'km'),col='darkgrey')
}

make_georaster<-function(motherraster){
  d.lat <- rasterToPoints(motherraster)
  d.lat[,3] =d.lat[,2]
  colnames(d.lat)[3] ="latitude"

  d.lon <- rasterToPoints(motherraster)
  d.lon[,3] =d.lon[,2]
  colnames(d.lon)[3] ="latitude"

  r.lat=rasterFromXYZ(d.lat)
  r.lon=rasterFromXYZ(d.lon)

  projection(r.lat)=projection(motherraster)
  projection(r.lon)=projection(motherraster)

  return(stack(list(latitude=r.lat,longitude=r.lon)))
}




#####**********************************************************************#####
#' generate a raster of densities
#'
#' @param longitude
#' @param latitude
#' @param refraster
#' @param dilute
#'
#' @return
#' @export
#'
#' @examples
getdensityraster <- function(longitude,latitude,refraster=NULL,dilute=5,method='bilinear'){

  require(raster)
  original=refraster
  if(is.null(refraster)){
    refraster<- rbioclim::getData(name = "worldclim",var='bio', res=2.5)
    refraster[!is.na(refraster)]=1
  }
  if(!is.null(dilute)){
  refraster<-aggregate(x = refraster,fact=dilute)
  }

coords = cbind(longitude,latitude)
sp = SpatialPoints(coords)

x <- rasterize(sp, refraster, fun='count')
x=x+5
x[is.na(x)]<-0

x <- resample(x, original, method)

return(x)
}

#####**********************************************************************#####
#' remove predictions where it does not exit
#'
#' @param myraster
#' @param maks
#' @param densitylimit
#' @param code
#'
#' @return
#' @export
#'
#' @examples
maskpresence <- function(myraster, mymask,densitylimit=0,code=NA){
  # maks[maks <= densitylimit]<-NA
  # myraster[is.na(maks)] <-  code
  updatemask=mymask
  updatemask[updatemask < densitylimit] <- NA
  myrasternew=myraster
  myrasternew[is.na(updatemask)]<-NA
  myrasternew<-raster::mask(myraster,updatemask)
return(myrasternew)
}


maskmap<-function(myraster, dens=T, envrange=F){
  rtmp<-myraster
  if(dens==T){
    if(!'aramask' %in% ls()){  mymask<-readRDS('maps/aramask.rda')  }
    mymask[is.na(mymask)] <-NA
    rtmp<-raster::mask(rtmp,mymask )
  }
  if(envrange==T){
    if(!'envmask' %in% ls()){  mymask<-readRDS('maps/envmask.rda')  }
    mymask[is.na(mymask)] <-NA
    rtmp<-raster::mask(rtmp,mymask )
  }
  return(rtmp)
}

#####**********************************************************************#####

#' gem
#'
#' @param X
#' @param env
#'
#' @return
#' @export
#'
#' @examples
gem<-function(X,env,envmap, cores=NULL,parallel=TRUE, n.trees=50, interaction.depth=2, shrinkage=0.1, n.minobsinnode=10, whatoutput='prob'){

  #----------------------------------------------------------------------------#
  stopifnot(whatoutput %in% c("raw",'prob'))
  stopifnot(class(X) %in% c("matrix",'data.frame','numeric'))

  if(class(X) == 'numeric'){
    nSNPs<-1
    SNPnames<- names(X)
    X<-data.frame(X=X, dummy=0)
  }else if(class(X) == 'matrix' |class(X) == 'data.frame' ){
    nSNPs<-ncol(X)
    SNPnames<- colnames(X)
  }else{
    stop('Dont know what to do')
  }

  message('Chose output from predictions: ', whatoutput)
  #----------------------------------------------------------------------------#
  require(caret)
  require(raster)

  # message('Using Cross-Validation to get the best parameter set in a random SNP')
  # i=sample(size=1, 1:ncol(X))
  # fitcontrol<-caret::trainControl(method='repeatedcv', number=10, repeats=3,verboseIter = FALSE)
  # gdmtrain<-caret::train( y =factor(X[,i]) , x= env , method = 'gbm', trControl = fitcontrol,verbose=FALSE)
  #
  # print(gdmtrain)
  # print(gdmtrain$bestTune)

  parlist<-data.frame(n.trees=n.trees, interaction.depth=interaction.depth, shrinkage=shrinkage, n.minobsinnode=n.minobsinnode)
  rownames(parlist)<-length(parlist)

  message('Model and prediction for every SNP...')

  if(parallel==TRUE){

  cl<-startcluster(cores)#<- start cluster
  clusterExport(cl=cl, list("envmap","X","parlist","env"),envir=environment())
  res<-(
    parLapply(
      cl,1:nSNPs,function(i){
                      message('modeling SNP ',i)
                      gdmfit<-caret::train( y =factor(X[,i]) , x= env , method = 'gbm', tuneGrid =  parlist, verbose = FALSE)
                      allelemap<-raster::predict(envmap,gdmfit ,type=whatoutput)
                      accu<-median(gdmfit$resample$Accuracy)
                      varimp<-as.character.factor(summary(gdmfit)[1,"var"])
                    return( list(map=allelemap,accu=accu, varimp=varimp) )
                      })
                    )
  stopCluster(cl) #<- stop cluster

  }else{
  res<-( lapply(1:nSNPs,function(i){
                  message('modeling SNP ',i)
                  gdmfit<-caret::train( y =factor(X[,i]) , x= env , method = 'gbm', tuneGrid =  parlist, verbose = FALSE)
                  allelemap<-predict(envmap,gdmfit ,type=whatoutput )
                  accu<-median(gdmfit$resample$Accuracy)
                  varimp<-as.character.factor(summary(gdmfit)[1,"var"])
                  return( list(map=allelemap,accu=accu, varimp=varimp) )
                  })
                )
  }


  message('...done')

  mapstack<- stack(sapply(res,function(i) i[['map']]))
  accuracylist<-unlist(sapply(res,function(i) i[['accu']]))
  importancelist<-unlist(sapply(res,function(i) i[['varimp']]))

  names(mapstack)<-SNPnames

  return(
          list(maps=mapstack , accuracy=accuracylist, importance=importancelist)
         )
}

  # else if(method=='rf)
  # require(randomForest)
  # rf<-randomForest( y =factor(X[,i]) , x= env )
  # adaptalleles<- predict(Euroclim,rf )

  # gemstack<-stack(sapply(1:ncol(X), FUN = function(a){
  #   rfm<-randomForest(
  #                 # x= dplyr::select(df, contains('bio'),contains('lon'), contains('lat')) ,
  #                 x= dplyr::select(df, contains('bio')) ,
  #                 y=factor(df$a),
  #                 ntree = 50
  #                 )
  # rfm
  # })




#' Title
#' #'
#' #' @param mydata
#' #' @param var
#' #' @param mypredictors
#' #' @param cross
#' #' @param seed
#' #' @param myfilename
#' #' @param type
#' #' @param sink
#' #' @param verbose
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#'
#' gem<-function(mydata,var,mypredictors,cross=5,seed=1, type="classification",sink=T, out='tmptables' ,filename="",verbose=T){
#' seed(seed)
#'
#' # checking arguments
#' if(!type %in%  c('regression','classification')){message('provide an appropriate random forest type')}
#'
#' # load llibraries
#' require(caret)
#' # require(randomForest)
#'
#' # generate the cross validation partitions
#' if(verbose) message(paste('creating',cross, 'partitions of the data ...'))
#' group<-createFolds(mydata[,'beta'] , k=cross,list=FALSE)
#'
#' # set up several parameters
#' allpredicted=c()
#' varcol=which(colanmes(mydata)==var)
#'
#' if(verbose) message(paste('running Random Forest...'))
#' rfmodels<-lapply(1:cross,FUN =
#'           function(l){
#'             message(paste('# ',l))
#'
#'             train <- mydata[group != l,]
#'             test <- mydata[group == l,-varcol]
#'
#'             if(type=="classification"){ mydata[,var] = factor(mydata[,var])
#'             }else{mydata[,var] = as.numeric(mydata[,var])}
#'
#'             myformula<-formula(paste(var," ~ ", paste(mypredictors,collapse = '+')))
#'
#'             fit <- train(myformula, data=mydata,method='rf', prox=TRUE )
#'             # fit <- randomForest(myformula, data=mydata,ntree=500,na.action=na.omit,keep.inbag = TRUE,importance=TRUE) # removed proximity
#'
#'             mypred <-predict(fit,test)
#'             allpredicted<-c(allpredicted,mypred)
#'
#'           return(fit)
#' } )
#'
#' # combine all random forest
#' if(verbose) message("merging the random forest")
#' rfall<-randomForest::combine(rfmodels[[1]],rfmodels[[2]],rfmodels[[3]],rfmodels[[4]],rfmodels[[5]])
#'
#' # Generate accuracy estimation
#' if(type=="classification"){
#'     myr2= table(ypredict == mydata[,phenoname])["TRUE"] / sum(table(ypredict == mydata[,phenoname]))
#'     if(verbose) message(paste("predictive power, proportion of right predictions=",myr2))
#' }else{
#'     myr2<-summary(lm(allpredicted~mydata[,phenoname]))$r.squared
#'     if(verbose) message(paste("predictive power, R2=",myr2))
#' }
#'
#' # sink info
#' if(sink==T){
#'     sortedimportance<-data.frame(importance(rfall) )
#'     sortedimportance$r2<-myr2
#'     sortedimportance$var<-row.names(sortedimportance)
#'     write.tsv(sortedimportance,file = paste0(out, '/',var,"_",filename,"randomforest.tsv"))
#' }
#'
#' return(rfall)
#' }
#'
#'
#'
#' #' Title
#' #'
#' @param mybioclim
#' @param xlim
#' @param ylim
#' @param replace
#' @param addPCA
#'
#' @return
#' @export
#'
#' @examples
cropenvironment <- function(mybioclim,xlim=c(-10.5,+ 53),ylim=c(32,65))  {
  require(raster)

  # xlim=c(-10.5,+ 35),ylim=c(34,65)

  Range=extent(c(xlim,ylim))
  EuropeClim = crop(mybioclim, Range)

  return(EuropeClim)
  }


#' #####**********************************************************************#####
#'
#' #' Title
#' #'
#' #' @param motherraster
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' make_georaster<-function(motherraster){
#' d.lat <- rasterToPoints(motherraster)
#' d.lat[,3] =d.lat[,2]
#' colnames(d.lat)[3] ="latitude"
#'
#' d.lon <- rasterToPoints(motherraster)
#' d.lon[,3] =d.lon[,2]
#' colnames(d.lon)[3] ="latitude"
#'
#' r.lat=rasterFromXYZ(d.lat)
#' r.lon=rasterFromXYZ(d.lon)
#'
#' projection(r.lat)=projection(motherraster)
#' projection(r.lon)=projection(motherraster)
#'
#' return(stack(list(latitude=r.lat,longitude=r.lon)))
#' }
#'
#'
#' #' Title
#' #'
#' #' @param bioclim
#' #' @param newlayerslist
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' morelayers_2_bioclim<-function(bioclim,newlayerslist){
#'
#' CLIMPC<-stack(bioclim,newlayerslist)
#' names(CLIMPC)<-c(names(bioclim),names(newlayerslist) )
#'
#' return(CLIMPC)
#' }
