## This is a R script to plot 1D number density distribution 
## using kernel density estimation

# Make the list of files
tmp_files = dir("../",pattern="^SD_output_ASCII_")
tmp = strsplit(tmp_files,"\\SD_output_ASCII_|\\.000.pe")
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

sink("numbers.txt")
sink()

# check the number of MPI processes
if(length(allmpiranks)!=(PRC_NUM_X*PRC_NUM_Y)){
	 cat(sprintf("Number of MPI processes is not consistent\n"))
	 cat(sprintf("%s != %s * %s\n",length(allmpiranks),PRC_NUM_X,PRC_NUM_Y))
	 quit()
}

# loop of time
for(time in alltimes){
#for(time in alltimes[which(as.numeric(alltimes)==1800)]){
	 cat(sprintf("processing the time = %s [sec]\n",time))

	 data <- NULL
	 # loop of MPI rank
	 for(mpirank in allmpiranks){
	 	file <- allfiles[mpirank,time]

		# Read data
		tmp <- read.table(paste("../",file,sep=""))
		data <- rbind(data, tmp)
	
	}
	valid_data <- data[data$V3>VALID2INVALID,]

	# number of SDs
	SD_num <- length(valid_data$V1)
	SD_num_dens <- SD_num/IMAX/JMAX/KMAX_sdm/length(allmpiranks)

	# number density of non-IN aerosols (soluble, (soluble+IN inactive mineral dust))
	nonIN_data <- valid_data[valid_data$V10<(-37.99999),]
	nonIN_num <- sum(nonIN_data$V11)
	nonIN_num_dens <- nonIN_num/IMAX/JMAX/KMAX_sdm/length(allmpiranks)/DX/DY/DZ
	
	# number density of IN active aerosols
	IN_data <- valid_data[valid_data$V10>(-37.99999),]
	IN_num <- sum(IN_data$V11)
	IN_num_dens <- IN_num/IMAX/JMAX/KMAX_sdm/length(allmpiranks)/DX/DY/DZ

	sink("numbers.txt", append=TRUE)
	cat(sprintf("###### SD/RD numbers at time = %s [s] ######\n",time))
	cat(sprintf("SD number = %d\n",SD_num))
	cat(sprintf("SD number density = %g [/grid]\n",SD_num_dens))
	cat(sprintf("total aerosol number = %g\n",nonIN_num+IN_num))
	cat(sprintf("total aerosol number density = %g [/m^3]\n",nonIN_num_dens+IN_num_dens))
	cat(sprintf("non IN aerosol number = %g\n",nonIN_num))
	cat(sprintf("non IN aerosol number density = %g [/m^3]\n",nonIN_num_dens))
	cat(sprintf("IN aerosol number = %g\n",IN_num))
	cat(sprintf("IN aerosol number density = %g [/m^3]\n\n",IN_num_dens))
	sink()

	# Probability Density of Freezing Temperature
	prob_dens_tf <- density(IN_data$V10,bw="nrd",weight=IN_data$V11/IN_num)
	#jpeg(paste("prob_dens_freezing_temp.",sprintf("%05d",as.numeric(time)),".jpg",sep=""))
	pdf(paste("prob_dens_freezing_temp.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
        plot(prob_dens_tf,
                   xlim=c(-38.0,-10.0),
		   main="Probability Density of Freezing Temperature ([/degC])",
                   xlab="Freezing temperature [degC])",
                   ylab="Probability Density [/degC])")
        dev.off()

	# Number density distribution of ammonium bisulfate\n (sum of nonIN and IN aerosols))"
	rho_amsul <- 1.78E+6 # density of ammonium bisulfate [g/m3]
	dry_radi_amsul <- (valid_data$V6/(4.0/3.0)/pi/rho_amsul)**(1.0/3.0)
	num_dens_amsul <- density(log(dry_radi_amsul),bw="nrd",weight=valid_data$V11/(nonIN_num+IN_num))
	num_dens_amsul$y <- (nonIN_num_dens+IN_num_dens)*num_dens_amsul$y
	num_dens_amsul$x <- exp(num_dens_amsul$x)
	#jpeg(paste("num_dens_amsul.",sprintf("%05d",as.numeric(time)),".jpg",sep=""))
	pdf(paste("num_dens_amsul.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
        plot(num_dens_amsul,
		   log="x",
		   xlim=c(1e-8,1e-6),
		   main="Number density distribution of ammonium bisulfate\n (sum of nonIN and IN aerosols))",
                   xlab="Dry radius of ammonium bisulfate [m])",
                   ylab="Number density distribution dN/dlogr ([/unit logr/m^3]")
        dev.off()

	# Mass density distribution of droplets"
	drop_data <- valid_data[valid_data$V12=="LIQ",]
	drop_radi <- drop_data$V5
	drop_num  <- sum(drop_data$V11)
	drop_num_dens <- drop_num/IMAX/JMAX/KMAX_sdm/length(allmpiranks)/DX/DY/DZ

	rho_liqw <- 1.0E+6 # density of liquid water [g/m3]
	mass_dens_drop <- density(log(drop_radi),bw="nrd",weight=drop_data$V11/drop_num)
	mass_dens_drop$x <- exp(mass_dens_drop$x)
	mass_dens_drop$y <- drop_num_dens*rho_liqw*(4.0/3.0)*pi*mass_dens_drop$x**3*mass_dens_drop$y
	pdf(paste("mass_dens_drop.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
        plot(mass_dens_drop,
		   log="x",
		   xlim=c(1e-9,1e-2),
		   main="Mass density distribution of droplets",
                   xlab="Radius of droplet [m])",
                   ylab="Mass density distribution of m(r)dN/dlogr ([g/unit logr/m^3]")
        dev.off()

	# terminal velocity of droplets
#	drop_data <- valid_data[valid_data$V12=="LIQ",]
#	drop_radi <- drop_data$V5
	drop_termv<- drop_data$V4 # terminal velocity of droplets [m/s]
	n_data_skip <- max(round(length(drop_termv)/10000),1)

	pdf(paste("term_vel_drop.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
        plot(2.0*drop_radi[seq(1, length(drop_radi), n_data_skip)],drop_termv[seq(1, length(drop_termv), n_data_skip)],
		   log="x",pch='.',
		   xlim=c(1e-9,1e-1),
		   main="Terminal velocity of droplets",
                   xlab="Diameter of droplet [m]",
                   ylab="Terminal velocity of droplets [m/s]")
        dev.off()

	# terminal velocity of ice
	ice_data <- valid_data[valid_data$V12=="ICE",]
	if(nrow(ice_data)>0){
		ice_radi <- mapply(max,ice_data$V7,ice_data$V8)
		ice_termv<- ice_data$V4 # terminal velocity of ice [m/s]
		n_data_skip <- max(round(length(ice_termv)/10000),1)

		pdf(paste("term_vel_ice.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
        	plot(2.0*ice_radi[seq(1, length(ice_radi), n_data_skip)],ice_termv[seq(1, length(ice_termv), n_data_skip)],
		   log="x",pch='.',
		   xlim=c(1e-9,1e-1),
		   ylim=c(0.0,max(ice_termv[seq(1, length(ice_termv), n_data_skip)])),
		   main="Terminal velocity of ice particles",
                   xlab="Dimension of ice particles [m]",
                   ylab="Terminal velocity of ice particles [m/s]")
	}else{
		pdf(paste("term_vel_ice.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
        	plot(1e-100,-1,type="n",log="x",
		   xlim=c(1e-9,1e-1),
		   ylim=c(0.0,20.0),
		   main="Terminal velocity of ice particles",
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

	# x-z position of SDs
	if(SD_num>0){
		valid_x <- valid_data$V1 # x coordinate of valid SD [m]
		valid_z <- valid_data$V3 # z coordinate of valid SD [m]
		n_data_skip <- max(round(length(valild_z)/10000),1)

		pdf(paste("xz_SD.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
        	plot(valid_x[seq(1, length(valid_x), n_data_skip)],valid_z[seq(1, length(valid_z), n_data_skip)],
		   pch='.',
		   xlim=c(0.0,PRC_NUM_X*DX*IMAX),
		   ylim=c(0.0,DZ*KMAX),
		   main="Positions of SDs on x-z plane",
                   xlab="x [m]",
                   ylab="z [m]")
        	dev.off()
	}
}
