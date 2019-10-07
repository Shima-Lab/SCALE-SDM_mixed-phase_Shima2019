## This is a R script to plot 1D number density distribution 
## using kernel density estimation
# Load required libraries
library(ncdf4)

# check arguments
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0) {
	cat(sprintf("No argument specified. Analyse all ice particles\n"))
	ice_category <- "all"
}else if( (args[1]=="graupel") || (args[1]=="snow") || (args[1]=="ice") || (args[1]=="all")) {
	cat(sprintf("Analyse %s\n",args[1]))
	ice_category <- args[1]
}else{
	cat(sprintf("Unsupported argument \"%s\" specified \n",args[1]))
	cat(sprintf("Use one of the following: graupel, snow, ice, all\n"))
	quit()
}

# Make the list of files
tmp_files = dir("../",pattern="^SD_output_NetCDF_")
tmp = strsplit(tmp_files,"\\SD_output_NetCDF_|\\.000.pe")
alltimes = unique(matrix(unlist(tmp),nrow=3)[2,])
allmpiranks = unique(matrix(unlist(tmp),nrow=3)[3,])

allfiles = matrix(tmp_files,ncol=length(alltimes))
rownames(allfiles) = allmpiranks
colnames(allfiles) = alltimes

# set parameters
VALID2INVALID <- -999.0 

# read parameters
tmp <- read.table(textConnection(gsub("="," ",gsub(",", " ", readLines("../run.conf")))),header=FALSE,fill=TRUE)
DX <- as.numeric(gsub(".[dD]",".e",tmp[tmp$V1=="DX",2]))
DY <- as.numeric(gsub(".[dD]",".e",tmp[tmp$V1=="DY",2]))
DZ <- as.numeric(gsub(".[dD]",".e",tmp[tmp$V1=="DZ",2]))
IMAX <- as.numeric(gsub(".[dD]",".e",tmp[tmp$V1=="IMAX",2]))
JMAX <- as.numeric(gsub(".[dD]",".e",tmp[tmp$V1=="JMAX",2]))
KMAX <- as.numeric(gsub(".[dD]",".e",tmp[tmp$V1=="KMAX",2]))
sdm_zlower <- as.numeric(gsub(".[dD]",".e",tmp[tmp$V1=="sdm_zlower",2]))
sdm_zupper <- as.numeric(gsub(".[dD]",".e",tmp[tmp$V1=="sdm_zupper",2]))
KMAX_sdm <- round((sdm_zupper - sdm_zlower)/DZ)
PRC_NUM_X <- as.numeric(gsub(".[dD]",".e",tmp[tmp$V1=="PRC_NUM_X",2])) 
PRC_NUM_Y <- as.numeric(gsub(".[dD]",".e",tmp[tmp$V1=="PRC_NUM_Y",2])) 

# check the number of MPI processes
if(length(allmpiranks)!=(PRC_NUM_X*PRC_NUM_Y)){
	 cat(sprintf("Number of MPI processes is not consistent\n"))
	 cat(sprintf("%s != %s * %s\n",length(allmpiranks),PRC_NUM_X,PRC_NUM_Y))
	 quit()
}

# loop of time
for(time in alltimes){
#for(time in alltimes[which(as.numeric(alltimes)==2400)]){
	 cat(sprintf("processing the time = %s [sec]\n",time))

	 data <- NULL
	 # loop of MPI rank
	 for(mpirank in allmpiranks){
		# Read data
	 	file <- allfiles[mpirank,time]
		ncin <- nc_open(paste("../",file,sep=""))	

		allvarnames <- attributes(ncin$var)$names
		data_tmp <- NULL
		for(varname in allvarnames){
			data_tmp <- cbind(data_tmp,ncvar_get(ncin,varname))
		}
		colnames(data_tmp) <- allvarnames
		data <- rbind(data,data_tmp)
		nc_close(ncin)
	}
	data <- data.frame(data)
	valid_data <- data[data$sd_z>VALID2INVALID,]

	# terminal velocity of ice
	all_ice_data <- valid_data[valid_data$sd_liqice==10,]

	# categorize ice particles
	mass <- all_ice_data$sdi_rho*(4.0/3.0)*pi*all_ice_data$sdi_re**2*all_ice_data$sdi_rp
	if( ice_category=="graupel" ){
	    typed_ice_id <- which( (mass*0.3)<all_ice_data$sdi_mrime )
	}else if( ice_category=="snow" ){
	    typed_ice_id <- which( ((mass*0.3)>=all_ice_data$sdi_mrime) & (all_ice_data$sdi_nmono>10) )
	}else if( ice_category=="ice" ){
	    typed_ice_id <- which( ((mass*0.3)>=all_ice_data$sdi_mrime) & (all_ice_data$sdi_nmono<=10) )
	}else if( ice_category=="all" ){
	    typed_ice_id <- which( mass>0.0 )
	}else{
	    cat(sprintf("Wrong ice category \"%s\" specified \n",typed_ice_id))
	    quit()
	}
	ice_data <- all_ice_data[typed_ice_id,]

	if(nrow(ice_data)>0){
		ice_radi <- mapply(max,ice_data$sdi_re,ice_data$sdi_rp)
		ice_termv<- ice_data$sd_vz # terminal velocity of ice [m/s]
		n_data_skip <- max(round(length(ice_termv)/10000),1)

		pdf(paste("term_vel_ice.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
        	plot(2.0*ice_radi[seq(1, length(ice_radi), n_data_skip)],ice_termv[seq(1, length(ice_termv), n_data_skip)],
		   log="x",pch='.',
		   xlim=c(1e-9,1e-1),
		   # ylim=c(0.0,max(ice_termv[seq(1, length(ice_termv), n_data_skip)])),
		   ylim=c(0.0,20.0),
		   main=paste("Terminal velocity of ice particles (T=",
		              sprintf("%05d",as.numeric(time)),
		              "s)"),
                   xlab="Dimension of ice particles [m]",
                   ylab="Terminal velocity of ice particles [m/s]")
	}else{
		pdf(paste("term_vel_ice.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
        	plot(1e-100,-1,type="n",log="x",
		   xlim=c(1e-9,1e-1),
		   ylim=c(0.0,20.0),
		   main=paste("Terminal velocity of ice particles (T=",
		              sprintf("%05d",as.numeric(time)),
		              "s)"),
                   xlab="Dimension of ice particles [m]",
                   ylab="Terminal velocity of ice particles [m/s]")
	}
	# plot fall velocities of typical ice particle types 
	linetype <- rep(1,10)  # all solid lines
	colortype <- seq(2,11) # skip 1 (black)
	n <- 0 

	n <- n+1
	legendtxt <- c("hail (M96)")
        plot(function(x){(1184.0/100.0)*(x*100.0)**0.5},5e-3,25e-3,lty=linetype[n],col=colortype[n],lwd=3,add=T)

	n <- n+1
	legendtxt <- append(legendtxt, c("lump graupel (LH74)"))
       	plot(function(x){1.3*(x*1000.0)**0.66},0.5e-3,3e-3,lty=linetype[n],col=colortype[n],lwd=3,add=T)

	n <- n+1
	legendtxt <- append(legendtxt, c("side planes (LH74)"))
       	plot(function(x){(0.81)*(x*1000.0)**0.99},0.4e-3,1.2e-3,lty=linetype[n],col=colortype[n],lwd=3,add=T)

	n <- n+1
	legendtxt <- append(legendtxt, c("plate (H72)"))
       	plot(function(x){(-0.55+90.69*(x*1000.0)-23.44*(x*1000.0)**2-5.26*(x*1000.0)**3)/100.0},1.0e-6,1.0e-3,lty=linetype[n],col=colortype[n],lwd=3,add=T)

	n <- n+1
	legendtxt <- append(legendtxt, c("aggregates of planar crystals (H02)"))
       	plot(function(x){(136.0/100.0)*(x*100.0)**0.4},6.0e-3,36.0e-3,lty=linetype[n],col=colortype[n],lwd=3,add=T)

	n <- n+1
	legendtxt <- append(legendtxt, c("dendrite (H72)"))
       	plot(function(x){(25.9/100.0)*(x*1000.0)**0.339 },1e-3,6e-3,lty=linetype[n],col=colortype[n],lwd=3,add=T)

	legend("topleft",
	lty=linetype[1:n],col=colortype[1:n],legend=legendtxt[1:n])

       	dev.off()
}
