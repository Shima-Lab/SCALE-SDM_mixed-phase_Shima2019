## This is a R script to plot 2D histogram of ice particles
# Load required libraries
library(fields)
library(ncdf4)
library(scales)

# check arguments
args = commandArgs(trailingOnly=TRUE)
if(length(args)!=0) {
	cat(sprintf("Use without any arguments\n"))
	quit()
}

# set parameters
VALID2INVALID <- -999.0

# Setting of bins
nbins = 100
min_radius = 5.0e-7 # [m]
max_radius = 1.0e-1 # [m]
min_maxdim = 2.0*min_radius
max_maxdim = 2.0*max_radius
min_density= 5.0e0 # [kg/m^3]
max_density= 1.1e+3 # [kg/m^3]
min_mass = 1.0e-16  # [kg]
max_mass = 1.0e-1   # [kg]
min_vel  = 1.0e-3   # [m/s]
max_vel  = 50.0e0   # [m/s]
min_aspct= 1.0e-2   # []
max_aspct= 1.0e+2   # []
min_massratio = 1.0e-4 # []
max_massratio = 1.2    # []

equr.bin = seq(log10(min_radius), log10(max_radius), length=nbins)
polr.bin = seq(log10(min_radius), log10(max_radius), length=nbins)
maxr.bin = seq(log10(min_radius), log10(max_radius), length=nbins)
maxdim.bin  = seq(log10(min_maxdim), log10(max_maxdim), length=nbins)
dens.bin = seq(log10(min_density),log10(max_density),length=nbins)
mass.bin = seq(log10(min_mass),log10(max_mass),length=nbins)
vel.bin = seq(log10(min_vel),log10(max_vel),length=nbins)
aspct.bin = seq(log10(min_aspct),log10(max_aspct),length=nbins)
massratio.bin = seq(log10(min_massratio),log10(max_massratio),length=nbins)

d_equr.bin = equr.bin[2] - equr.bin[1]
d_polr.bin = polr.bin[2] - polr.bin[1]
d_maxr.bin = maxr.bin[2] - maxr.bin[1]
d_maxdim.bin  = maxdim.bin[2]  - maxdim.bin[1]
d_dens.bin = dens.bin[2] - dens.bin[1]
d_mass.bin = mass.bin[2] - mass.bin[1]
d_vel.bin = vel.bin[2] - vel.bin[1]
d_aspct.bin = aspct.bin[2] - aspct.bin[1]
d_massratio.bin = massratio.bin[2] - massratio.bin[1]

# Make the list of files
tmp_files = dir("../",pattern="^SD_output_NetCDF_")
tmp = strsplit(tmp_files,"\\SD_output_NetCDF_|\\.000.pe")
alltimes = unique(matrix(unlist(tmp),nrow=3)[2,])
allmpiranks = unique(matrix(unlist(tmp),nrow=3)[3,])

allfiles = matrix(tmp_files,ncol=length(alltimes))
rownames(allfiles) = allmpiranks
colnames(allfiles) = alltimes

## set color pallets

pal <- c(colorRampPalette(c(rgb(0,0,0.7,0.2), rgb(0,0,1,0.9)), alpha = TRUE),
         colorRampPalette(c(rgb(0,0.7,0,0.2), rgb(0,1,0,0.9)), alpha = TRUE),
         colorRampPalette(c(rgb(0.7,0,0,0.2), rgb(1,0,0,0.9)), alpha = TRUE))
names(pal)=c("ice","snow","graupel")

# loop of time
for(time in alltimes){
#for(time in alltimes[which(as.numeric(alltimes)==3000)]){
    cat(sprintf("processing the time = %s [sec]\n",time))

    # open pdf files and draw the axes
    pdf(paste("aspect.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
    dev_aspect <- dev.cur()
    image(maxdim.bin+d_maxdim.bin/2.0, aspct.bin+d_aspct.bin/2.0, diag(nbins)*0.0e0,
		   zlim=c(4.0,7.0),
    		   col = rgb(1,1,1,0),
		   main=paste("Aspect Ratio Distribution (T=",
		              sprintf("%05d",as.numeric(time)),
		              "s)\n(Mass Density log([kg/unit_log10(max_D)/unit_log10(phi)]))"),
		   xlab="Size (maximum dimension) log10([m])",
		   ylab="Aspect ratio (c/a) log10([m])")

    pdf(paste("density.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
    dev_density <- dev.cur()
    image(maxdim.bin+d_maxdim.bin/2.0, dens.bin+d_dens.bin/2.0, diag(nbins)*0.0e0,
		   zlim=c(4.0,7.0),
    		   col = rgb(1,1,1,0),
		   main=paste("Size-Density Distribution (T=",
		              sprintf("%05d",as.numeric(time)),
		              "s)\n(Mass Density log([kg/unit_log10(max_D)/unit_log10(density)]))"),
		   xlab="Size (maximum dimension) log10([m])",
		   ylab="Density of ice particle log10([kg/m^3])")

#    pdf(paste("mass.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
#    dev_mass <- dev.cur()
#    image(maxdim.bin+d_maxdim.bin/2.0, mass.bin+d_mass.bin/2.0, diag(nbins)*0.0e0,
#		   zlim=c(4.0,7.0), #nlevel=100,
#    		   col = rgb(1,1,1,0),
#		   main=paste("Mass-Dimension Distribution (T=",
#		              sprintf("%05d",as.numeric(time)),
#                              "s)\n(Mass Density log([kg/unit_log10(max_D)/unit_log10(mass)]))"),
#		   xlab="Size (maximum dimension) log10([m])",
#		   ylab="Mass of ice particle log10([kg])")

    pdf(paste("massratio.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
    dev_massratio <- dev.cur()
    image(maxdim.bin+d_maxdim.bin/2.0, massratio.bin+d_massratio.bin/2.0, diag(nbins)*0.0e0,
		   zlim=c(4.0,7.0), #nlevel=100,
    		   col = rgb(1,1,1,0),
		   main=paste("Mass-Dimension Distribution (T=",
		              sprintf("%05d",as.numeric(time)),
                              "s)\n(Mass Density log([kg/unit_log10(max_D)/unit_log10(massratio)]))"),
		   xlab="Size (maximum dimension) log10([m])",
		   ylab="Normalized mass of ice particle log10([])")

    pdf(paste("term_vel_ice.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
    dev_vel <- dev.cur()
    image(maxdim.bin+d_maxdim.bin/2.0, vel.bin+d_vel.bin/2.0, diag(nbins)*0.0e0,
		   zlim=c(4.0,7.0), #nlevel=100,
    		   col = rgb(1,1,1,0),
		   main=paste("Terminal velocity of ice particles (T=",
				     sprintf("%05d",as.numeric(time)),
				     "s)\n(Mass Density log([kg/unit_log10(max_D)/unit_log10(velocity)]))"),
		   xlab="Dimension of ice particles log10([m])",
		   ylab="Terminal velocity of ice particles log10([m/s])")

#    for(ice_category in c("snow")){
    for(ice_category in c("graupel", "ice", "snow")){

    	 # initialize histogram
    	 aspect_dist <- diag(nbins)*0.0 + 1.0e-100
    	 density_dist <- diag(nbins)*0.0 + 1.0e-100
    	 mass_dist <- diag(nbins)*0.0 + 1.0e-100
    	 vel_dist <- diag(nbins)*0.0 + 1.0e-100
    	 massratio_dist <- diag(nbins)*0.0 + 1.0e-100

	 # loop of MPI rank
	 for(mpirank in allmpiranks){
#	 for(mpirank in allmpiranks[80]){
		# open file
	 	file <- allfiles[mpirank,time]
		ncin <- nc_open(paste("../",file,sep=""))

		# extract the id of ice particles
		all_sd_z      <- ncvar_get(ncin,"sd_z")
		all_sd_liqice <- ncvar_get(ncin,"sd_liqice")
		ice_id <- which(all_sd_liqice==10 & all_sd_z>VALID2INVALID)

		## equatorial radius
		all_sdi_re     <- ncvar_get(ncin,"sdi_re")
		equr <- all_sdi_re[ice_id]
		## polar radius
		all_sdi_rp     <- ncvar_get(ncin,"sdi_rp")
		polr <- all_sdi_rp[ice_id]
		## density of particle
		all_sdi_rho     <- ncvar_get(ncin,"sdi_rho")
		dens <- all_sdi_rho[ice_id]
		## terminal velocity of ice particles
		all_sdi_tvel     <- ncvar_get(ncin,"sd_vz")
		icetvel <- all_sdi_tvel[ice_id]
		## multiplicity of super-droplet
		all_sd_n     <- ncvar_get(ncin,"sd_n")
		mult <- all_sd_n[ice_id]
		## rime mass
		all_sdi_mrime     <- ncvar_get(ncin,"sdi_mrime")
		mrime <- all_sdi_mrime[ice_id]
		## number of monomers
		all_sdi_nmono     <- ncvar_get(ncin,"sdi_nmono")
		nmono <- all_sdi_nmono[ice_id]

		# close file
		nc_close(ncin)

		# categorize ice particles
		mass <- dens*(4.0/3.0)*pi*equr**2*polr
		if( ice_category=="graupel" ){
		    typed_ice_id <- which( (mass*0.3)<mrime )
		}else if( ice_category=="snow" ){
		    typed_ice_id <- which( ((mass*0.3)>=mrime) & (nmono>10) )
		}else if( ice_category=="ice" ){
		    typed_ice_id <- which( ((mass*0.3)>=mrime) & (nmono<=10) )
		}else if( ice_category=="all" ){
		    typed_ice_id <- which( mass>0.0 )
		}else{
		    cat(sprintf("Wrong ice category \"%s\" specified \n",typed_ice_id))
	    	    quit()
		}

		# if no ice exists, skip to the next MPI rank
		if(length(typed_ice_id)==0){ next }

		equr <- equr[typed_ice_id]
		polr <- polr[typed_ice_id]
		dens <- dens[typed_ice_id]
		mult <- mult[typed_ice_id]
		icetvel <- icetvel[typed_ice_id]

		# calculate the mass density distribution of aspect ratio
		ptl_num <- length(mult)
		maxdim <- numeric(ptl_num)
		imaxdim <- numeric(ptl_num)
		iaspct <- numeric(ptl_num)
		maxdim[1:ptl_num] <- 2.0*pmax(equr[1:ptl_num],polr[1:ptl_num])
		imaxdim[1:ptl_num] <- findInterval(log10(maxdim[1:ptl_num]), maxdim.bin, all.inside=TRUE)
		iaspct[1:ptl_num] <- findInterval(log10(polr[1:ptl_num]/equr[1:ptl_num]), aspct.bin, all.inside=TRUE)
		for (n in 1:ptl_num){
    		    aspect_dist[imaxdim[n],iaspct[n]] <- aspect_dist[imaxdim[n],iaspct[n]] +
		    			        mult[n]*dens[n]*(4.0/3.0)*pi*equr[n]**2*polr[n]
		}
		
		# calculate the mass density distribution of desity-dimension
		ptl_num <- length(mult)
		maxdim <- numeric(ptl_num)
		imaxdim <- numeric(ptl_num)
		idens <- numeric(ptl_num)
		maxdim[1:ptl_num] <- 2.0*pmax(equr[1:ptl_num],polr[1:ptl_num])
		imaxdim[1:ptl_num] <- findInterval(log10(maxdim[1:ptl_num]), maxdim.bin, all.inside=TRUE)
		idens[1:ptl_num] <- findInterval(log10(dens[1:ptl_num]), dens.bin, all.inside=TRUE)
		for (n in 1:ptl_num){
    		    density_dist[imaxdim[n],idens[n]] <- density_dist[imaxdim[n],idens[n]] +
		    			        mult[n]*dens[n]*(4.0/3.0)*pi*equr[n]**2*polr[n]
		}

#		# calculate the mass density distribution of mass-dimension
#		ptl_num <- length(mult)
#		maxdim <- numeric(ptl_num)
#		mass <- numeric(ptl_num)
#		imaxdim <- numeric(ptl_num)
#		imass <- numeric(ptl_num)
#		maxdim[1:ptl_num] <- 2.0*pmax(equr[1:ptl_num],polr[1:ptl_num])
#		mass[1:ptl_num] <- dens[1:ptl_num]*(4.0/3.0)*pi*equr[1:ptl_num]**2*polr[1:ptl_num]
#		imaxdim[1:ptl_num] <- findInterval(log10(maxdim[1:ptl_num]), maxdim.bin, all.inside=TRUE)
#		imass[1:ptl_num] <- findInterval(log10(mass[1:ptl_num]), mass.bin, all.inside=TRUE)
#		for (n in 1:ptl_num){
#    		    mass_dist[imaxdim[n],imass[n]] <- mass_dist[imaxdim[n],imass[n]] +
#		    			        mult[n]*dens[n]*(4.0/3.0)*pi*equr[n]**2*polr[n]
#		}

		# calculate the mass density distribution of mass-dimension (massratio)
		ptl_num <- length(mult)
		maxdim <- numeric(ptl_num)
		massratio <- numeric(ptl_num)
		imaxdim <- numeric(ptl_num)
		imassratio <- numeric(ptl_num)
		maxdim[1:ptl_num] <- 2.0*pmax(equr[1:ptl_num],polr[1:ptl_num])
		massratio[1:ptl_num] <- (dens[1:ptl_num]*(4.0/3.0)*pi*equr[1:ptl_num]**2*polr[1:ptl_num])/(916.8*(pi/6.0)*maxdim[1:ptl_num]**3)
		imaxdim[1:ptl_num] <- findInterval(log10(maxdim[1:ptl_num]), maxdim.bin, all.inside=TRUE)
		imassratio[1:ptl_num] <- findInterval(log10(massratio[1:ptl_num]), massratio.bin, all.inside=TRUE)
		for (n in 1:ptl_num){
    		    massratio_dist[imaxdim[n],imassratio[n]] <- massratio_dist[imaxdim[n],imassratio[n]] +
		    			        mult[n]*dens[n]*(4.0/3.0)*pi*equr[n]**2*polr[n]
		}

		# calculate the mass density distribution of terminal velocity-dimension
		ptl_num <- length(mult)
		maxdim <- numeric(ptl_num)
		imaxdim <- numeric(ptl_num)
		itvel <- numeric(ptl_num)
		maxdim[1:ptl_num] <- 2.0*pmax(equr[1:ptl_num],polr[1:ptl_num])
		imaxdim[1:ptl_num] <- findInterval(log10(maxdim[1:ptl_num]), maxdim.bin, all.inside=TRUE)
		itvel[1:ptl_num] <- findInterval(log10(icetvel[1:ptl_num]), vel.bin, all.inside=TRUE)
		for (n in 1:ptl_num){
    		    vel_dist[imaxdim[n],itvel[n]] <- vel_dist[imaxdim[n],itvel[n]] +
		    			        mult[n]*dens[n]*(4.0/3.0)*pi*equr[n]**2*polr[n]
		}

	}
	aspect_dist[,] <- aspect_dist[,]/d_maxdim.bin/d_aspct.bin
	density_dist[,] <- density_dist[,]/d_maxdim.bin/d_dens.bin
#	mass_dist[,] <- mass_dist[,]/d_maxdim.bin/d_mass.bin
	vel_dist[,] <- vel_dist[,]/d_maxdim.bin/d_vel.bin
	massratio_dist[,] <- massratio_dist[,]/d_maxdim.bin/d_massratio.bin
	
	# plot the equatorial and polar radius distribution
	dev.set(dev_aspect)
	par(new=T)
	image(maxdim.bin+d_maxdim.bin/2.0, aspct.bin+d_aspct.bin/2.0, pmin(
			  	    apply(log10(aspect_dist),c(1,2),function(x){ifelse(x>(-90.0),max(x,4.0),x)}), 
				    7.0),
		   col=pal[[ice_category]](100),
		   zlim=c(4.0,7.0), # nlevel=100,
		   main=paste("Aspect Ratio Distribution (T=",
		              sprintf("%05d",as.numeric(time)),
		              "s)\n(Mass Density log([kg/unit_log10(max_D)/unit_log10(phi)]))"),
		   xlab="Size (maximum dimension) log10([m])",
		   ylab="Aspect ratio (c/a) log10([m])",
		   ann=F,useRaster=TRUE)

	# plot the mass density distribution of size-density relation
	dev.set(dev_density)
	par(new=T)
	image(maxdim.bin+d_maxdim.bin/2.0, dens.bin+d_dens.bin/2.0, pmin(
			  	    apply(log10(density_dist),c(1,2),function(x){ifelse(x>(-90.0),max(x,4.0),x)}), 
				    7.0),
		   col=pal[[ice_category]](100),
		   zlim=c(4.0,7.0), # nlevel=100,
		   main=paste("Size-Density Distribution (T=",
		              sprintf("%05d",as.numeric(time)),
		              "s)\n(Mass Density log([kg/unit_log10(max_D)/unit_log10(density)]))"),
		   xlab="Size (maximum dimension) log10([m])",
		   ylab="Density of ice particle log10([kg/m^3])",
		   ann=F,useRaster=TRUE)

#	# plot the mass density distribution of mass-dimension relation
#	dev.set(dev_mass)
#	par(new=T)
#	image(maxdim.bin+d_maxdim.bin/2.0, mass.bin+d_mass.bin/2.0, pmin(
#			  	    apply(log10(mass_dist),c(1,2),function(x){ifelse(x>(-90.0),max(x,4.0),x)}), 
#				    7.0)								  ,
#		   col=pal[[ice_category]](100),
#		   zlim=c(4.0,7.0), #nlevel=100,
#		   main=paste("Mass-Dimension Distribution (T=",
#		              sprintf("%05d",as.numeric(time)),
#                              "s)\n(Mass Density log([kg/unit_log10(max_D)/unit_log10(mass)]))"),
#		   xlab="Size (maximum dimension) log10([m])",
#		   ylab="Mass of ice particle log10([kg])",
#		   ann=F,useRaster=TRUE)

	# plot the mass density distribution of mass-dimension relation (massratio)
	dev.set(dev_massratio)
	par(new=T)
	image(maxdim.bin+d_maxdim.bin/2.0, massratio.bin+d_massratio.bin/2.0, pmin(
			  	    apply(log10(massratio_dist),c(1,2),function(x){ifelse(x>(-90.0),max(x,4.0),x)}), 
				    7.0)								  ,
		   col=pal[[ice_category]](100),
		   zlim=c(4.0,7.0), #nlevel=100,
		   main=paste("Mass-Dimension Distribution (T=",
		              sprintf("%05d",as.numeric(time)),
                              "s)\n(Mass Density log([kg/unit_log10(max_D)/unit_log10(massratio)]))"),
		   xlab="Size (maximum dimension) log10([m])",
		   ylab="Normalized mass of ice particle log10([])",
		   ann=F,useRaster=TRUE)

	# plot the mass density distribution of the "terminal velocity"-dimension relation
	dev.set(dev_vel)
	par(new=T)
	image(maxdim.bin+d_maxdim.bin/2.0, vel.bin+d_vel.bin/2.0, pmin(
			  	    apply(log10(vel_dist),c(1,2),function(x){ifelse(x>(-90.0),max(x,4.0),x)}), 
				    7.0),
		   col=pal[[ice_category]](100),
		   zlim=c(4.0,7.0), #nlevel=100,
		   main=paste("Terminal velocity of ice particles (T=",
				     sprintf("%05d",as.numeric(time)),
				     "s)\n(Mass Density log([kg/unit_log10(max_D)/unit_log10(velocity)]))"),
		   xlab="Dimension of ice particles log10([m])",
		   ylab="Terminal velocity of ice particles log10([m/s])",
		   ann=F,useRaster=TRUE)

    }

#    # plot the mass density distribution of mass-dimension relation
#    dev.set(dev_mass)
#    # plot various typical mass-dimension relations
#    linetype <- rep(1,20)  # all solid lines
#    colortype <- seq(2,21) # skip 1 (black)
#    n <- 0 
#    n <- n+1
#    legendtxt <- c("ice spheres")
#    plot(function(x){log10(0.480*1000.0*(10**x)**3)},log10(min_maxdim),log10(max_maxdim),lty=linetype[n],col=colortype[n],lwd=2,add=T)
#
#    n <- n+1
#    legendtxt <- append(legendtxt,c("hail (M96)"))
#    plot(function(x){log10(0.466*1000.0*(10**x)**3)},log10(0.5e-2),log10(2.5e-2),lty=linetype[n],col=colortype[n],lwd=2,add=T)
#
#    n <- n+1
#    legendtxt <- append(legendtxt,c("lump graupel (M96)"))
#    plot(function(x){log10(0.049e-3*((10**x)*100.0)**2.8)},log10(0.05e-2),log10(0.3e-2),lty=linetype[n],col=colortype[n],lwd=2,add=T)
#
#    n <- n+1
#    legendtxt <- append(legendtxt,c("densely rimed dendrites (M96)"))
#    plot(function(x){log10(0.003e-3*((10**x)*100.0)**2.3)},log10(0.18e-2),log10(0.4e-2),lty=linetype[n],col=colortype[n],lwd=2,add=T)
#
##	n <- n+1
##	legendtxt <- append(legendtxt,c("aggregates (SK10)"))
##        plot(function(x){log10(0.00183e-3*((10**x)*100.0)**2.04)},log10(min_maxdim),log10(max_maxdim),lty=linetype[n],col=colortype[n],lwd=2,add=T)
#
##	n <- n+1
##	legendtxt <- append(legendtxt,c("aggregates of unrimed side planes (LH74)"))
##        plot(function(x){log10(0.04e-6*((10**x)*1000.0)**1.4)},log10(0.05e-2),log10(0.4e-2),lty=linetype[n],col=colortype[n],lwd=2,add=T)
#
#    n <- n+1
#    legendtxt <- append(legendtxt,c("aggregates (M96)"))
#    plot(function(x){log10(0.0028e-3*((10**x)*100.0)**2.1)},log10(0.08e-2),log10(0.45e-2),lty=linetype[n],col=colortype[n],lwd=2,add=T)
#
#    n <- n+1
#    legendtxt <- append(legendtxt,c("hexagonal columns (small) (M96)"))
#    plot(function(x){log10(0.1677e-3*((10**x)*100.0)**2.91)},log10(0.003e-2),log10(0.01e-2),lty=linetype[n],col=colortype[n],lwd=2,add=T)
#
#    n <- n+1
#    legendtxt <- append(legendtxt,c("hexagonal columns (middle) (M96)"))
#    plot(function(x){log10(0.00166e-3*((10**x)*100.0)**1.91)},log10(0.01e-2),log10(0.03e-2),lty=linetype[n],col=colortype[n],lwd=2,add=T)
#
#    n <- n+1
#    legendtxt <- append(legendtxt,c("hexagonal columns (big) (M96)"))
#    plot(function(x){log10(0.000907e-3*((10**x)*100.0)**1.74)},log10(0.03e-2),log10(0.3e-2),lty=linetype[n],col=colortype[n],lwd=2,add=T)
#
#    n <- n+1
#    legendtxt <- append(legendtxt,c("hexagonal plates (M96)"))
#    plot(function(x){log10(0.00739e-3*((10**x)*100.0)**2.45)},log10(0.0015e-2),log10(0.3e-2),lty=linetype[n],col=colortype[n],lwd=2,add=T)
#
#    n <- n+1
#    legendtxt <- append(legendtxt,c("broad-branched crystals (small) (M96)"))
#    plot(function(x){log10(0.00583e-3*((10**x)*100.0)**2.42)},log10(0.001e-2),log10(0.01e-2),lty=linetype[n],col=colortype[n],lwd=2,add=T)
#
#    n <- n+1
#    legendtxt <- append(legendtxt,c("broad-branched crystals (big) (M96)"))
#    plot(function(x){log10(0.000516e-3*((10**x)*100.0)**1.8)},log10(0.01e-2),log10(0.1e-2),lty=linetype[n],col=colortype[n],lwd=2,add=T)
#
#    n <- n+1
#    legendtxt <- append(legendtxt,c("vapor spheres")) # 1e-3 kg/m^3 is the saturated vapor density at T=-20C
#    plot(function(x){log10(1e-3*(pi/6.0)*(10**x)**3)},log10(min_maxdim),log10(max_maxdim),lty=linetype[n],col=colortype[n],lwd=2,add=T)
#
#    op <- par(cex = 0.8)
#    legend("topleft",
#    lty=linetype[1:n],col=colortype[1:n],legend=legendtxt[1:n],
#    bty="n")
#    op <- par(cex = 1.0)

    # plot the mass density distribution of mass-dimension relation (massratio)
    dev.set(dev_massratio)
    # plot various typical mass-dimension relations (massratoi)
    linetype <- seq(1,50)
    line_alpha <- 1.0
    rainbowalpha <- alpha(rainbow(50),line_alpha)
    colortype <- rep(1,50)

    n <- 0
    legendtxt <- character()
 
#    n <- n+1
#    legendtxt <- c("ice spheres")
#    plot(function(x){log10((0.480*1000.0*(10**x)**3)/(0.480*1000.0*(10**x)**3))},log10(min_maxdim),log10(max_maxdim),lty=linetype[n],col=colortype[n],lwd=2,add=T)

    n <- n+1
    legendtxt <- append(legendtxt,c("hexagonal column (M96)"))
    colortype[n] <- rainbowalpha[n+28]
    plot(function(x){log10((0.1677e-3*((10**x)*100.0)**2.91)/(0.480*1000.0*(10**x)**3))},log10(0.003e-2),log10(0.01e-2),lty=1,col=colortype[n],lwd=3,add=T)

    plot(function(x){log10((0.00166e-3*((10**x)*100.0)**1.91)/(0.480*1000.0*(10**x)**3))},log10(0.01e-2),log10(0.03e-2),lty=1,col=colortype[n],lwd=3,add=T)

    plot(function(x){log10((0.000907e-3*((10**x)*100.0)**1.74)/(0.480*1000.0*(10**x)**3))},log10(0.03e-2),log10(0.3e-2),lty=1,col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt,c("hexagonal plate (M96)"))
    colortype[n] <- rainbowalpha[n+28]
    plot(function(x){log10((0.00739e-3*((10**x)*100.0)**2.45)/(0.480*1000.0*(10**x)**3))},log10(0.0015e-2),log10(0.3e-2),lty=1,col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt,c("broad-branched crystal (M96)"))
    colortype[n] <- rainbowalpha[n+28]
    plot(function(x){log10((0.00583e-3*((10**x)*100.0)**2.42)/(0.480*1000.0*(10**x)**3))},log10(0.001e-2),log10(0.01e-2),lty=1,col=colortype[n],lwd=3,add=T)

    plot(function(x){log10((0.000516e-3*((10**x)*100.0)**1.8)/(0.480*1000.0*(10**x)**3))},log10(0.01e-2),log10(0.1e-2),lty=1,col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt,c("crystal with sector like branches (M96)"))
    colortype[n] <- rainbowalpha[n+28]
    plot(function(x){log10((0.00614e-3*((10**x)*100.0)**2.42)/(0.480*1000.0*(10**x)**3))},log10(10.0e-6),log10(40.0e-6),lty=1,col=colortype[n],lwd=3,add=T)

    plot(function(x){log10((0.00142e-3*((10**x)*100.0)**2.02)/(0.480*1000.0*(10**x)**3))},log10(40.0e-6),log10(2000.0e-6),lty=1,col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt,c("stellar crystal with broad arms (M96)"))
    colortype[n] <- rainbowalpha[n+28]
    plot(function(x){log10((0.00583e-3*((10**x)*100.0)**2.42)/(0.480*1000.0*(10**x)**3))},log10(10.0e-6),log10(90.0e-6),lty=1,col=colortype[n],lwd=3,add=T)

    plot(function(x){log10((0.000270e-3*((10**x)*100.0)**1.67)/(0.480*1000.0*(10**x)**3))},log10(90.0e-6),log10(1500.0e-6),lty=1,col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt,c("side plane (M96)"))
    colortype[n] <- rainbowalpha[n+28]
    plot(function(x){log10((0.00419e-3*((10**x)*100.0)**2.3)/(0.480*1000.0*(10**x)**3))},log10(300.0e-6),log10(2500.0e-6),lty=1,col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt,c("dendritic crystal (HK87)"))
    colortype[n] <- rainbowalpha[n+28]
    plot(function(x){log10((6.12e-7*((10**x)*100.0)**2.29)/(0.480*1000.0*(10**x)**3))},log10(0.6e-3),log10(5.3e-3),lty=1,col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt,c("dendritic crystal (K89)"))
    colortype[n] <- rainbowalpha[n+28]
    plot(function(x){log10((4.82e-7*((10**x)*100.0)**1.97)/(0.480*1000.0*(10**x)**3))},log10(0.8e-3),log10(6.8e-3),lty=1,col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt,c("elementary needle (M90)"))
    colortype[n] <- rainbowalpha[n+28]
    plot(function(x){log10((0.0049e-6*((10**x)*1000.0)**1.8)/(0.480*1000.0*(10**x)**3))},log10(0.6e-3),log10(2.7e-3),lty=1,col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt,c("hail (M96)"))
    colortype[n] <- rainbowalpha[n-8]
    plot(function(x){log10((0.466*1000.0*(10**x)**3)/(0.480*1000.0*(10**x)**3))},log10(0.5e-2),log10(2.5e-2),lty=1,col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt,c("lump graupel (M96)"))
    colortype[n] <- rainbowalpha[n-8]
    plot(function(x){log10((0.049e-3*((10**x)*100.0)**2.8)/(0.480*1000.0*(10**x)**3))},log10(0.05e-2),log10(0.3e-2),lty=1,col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt,c("lump graupel (HK87)"))
    colortype[n] <- rainbowalpha[n-8]
    plot(function(x){log10((1.07e-4*((10**x)*100.0)**3.10)/(0.480*1000.0*(10**x)**3))},log10(0.4e-3),log10(9.0e-3),lty=1,col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt,c("densely rimed plate (HK87)"))
    colortype[n] <- rainbowalpha[n-8]
    plot(function(x){log10((9.53e-5*((10**x)*100.0)**3.8)/(0.480*1000.0*(10**x)**3))},log10(0.7e-3),log10(2.2e-3),lty=1,col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt,c("densely rimed dendrite (M96)"))
    colortype[n] <- rainbowalpha[n-8]
    plot(function(x){log10((0.003e-3*((10**x)*100.0)**2.3)/(0.480*1000.0*(10**x)**3))},log10(0.18e-2),log10(0.4e-2),lty=1,col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt,c("densely rimed stellar (HK87)"))
    colortype[n] <- rainbowalpha[n-8]
    plot(function(x){log10((7.55e-6*((10**x)*100.0)**3.04)/(0.480*1000.0*(10**x)**3))},log10(1.1e-3),log10(4.7e-3),lty=1,col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt,c("rimed long column (M96)"))
    colortype[n] <- rainbowalpha[n-8]
    plot(function(x){log10((0.00145e-3*((10**x)*100.0)**1.8)/(0.480*1000.0*(10**x)**3))},log10(200.0e-6),log10(2400.0e-6),lty=1,col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt,c("rimed needle crystal (M90)"))
    colortype[n] <- rainbowalpha[n-8]
    plot(function(x){log10((0.006e-6*((10**x)*1003.0)**2.1)/(0.480*1000.0*(10**x)**3))},log10(0.5e-3),log10(2.8e-3),lty=1,col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt,c("aggregate of unrimed side planes (LH74)"))
    colortype[n] <- rainbowalpha[2*(n-17)+15]
    plot(function(x){log10((0.04e-6*((10**x)*1000.0)**1.4)/(0.480*1000.0*(10**x)**3))},log10(0.05e-2),log10(0.4e-2),lty=1,col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt,c("aggregate of side planes (M96)"))
    colortype[n] <- rainbowalpha[2*(n-17)+15]
    plot(function(x){log10((0.0033e-3*((10**x)*100.0)**2.2)/(0.480*1000.0*(10**x)**3))},log10(600.0e-6),log10(4100.0e-6),lty=1,col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt,c("aggregate of side planes, columns, bullets (M96)"))
    colortype[n] <- rainbowalpha[2*(n-17)+15]
    plot(function(x){log10((0.0028e-3*((10**x)*100.0)**2.1)/(0.480*1000.0*(10**x)**3))},log10(0.08e-2),log10(0.45e-2),lty=1,col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt,c("aggregate of dendrites (LH74)"))
    colortype[n] <- rainbowalpha[2*(n-17)+15]
    plot(function(x){log10((0.037e-6*((10**x)*1000.0)**1.4)/(0.480*1000.0*(10**x)**3))},log10(2.0e-3),log10(12.0e-3),lty=1,col=colortype[n],lwd=3,add=T)

    ### draw the line at the center
    n <- 0 
    n <- n+1
#    legendtxt <- c("hexagonal column (M96)")
    linecenter_col <- c("white")
    plot(function(x){log10((0.1677e-3*((10**x)*100.0)**2.91)/(0.480*1000.0*(10**x)**3))},log10(0.003e-2),log10(0.01e-2),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    plot(function(x){log10((0.00166e-3*((10**x)*100.0)**1.91)/(0.480*1000.0*(10**x)**3))},log10(0.01e-2),log10(0.03e-2),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    plot(function(x){log10((0.000907e-3*((10**x)*100.0)**1.74)/(0.480*1000.0*(10**x)**3))},log10(0.03e-2),log10(0.3e-2),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt,c("hexagonal plate (M96)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10((0.00739e-3*((10**x)*100.0)**2.45)/(0.480*1000.0*(10**x)**3))},log10(0.0015e-2),log10(0.3e-2),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt,c("broad-branched crystal (M96)"))
    linecenter_col <- append(linecenter_col,"white")
    plot(function(x){log10((0.00583e-3*((10**x)*100.0)**2.42)/(0.480*1000.0*(10**x)**3))},log10(0.001e-2),log10(0.01e-2),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    plot(function(x){log10((0.000516e-3*((10**x)*100.0)**1.8)/(0.480*1000.0*(10**x)**3))},log10(0.01e-2),log10(0.1e-2),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt,c("crystal with sector like branches (M96)"))
    linecenter_col <- append(linecenter_col,c("black"))
    plot(function(x){log10((0.00614e-3*((10**x)*100.0)**2.42)/(0.480*1000.0*(10**x)**3))},log10(10.0e-6),log10(40.0e-6),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    plot(function(x){log10((0.00142e-3*((10**x)*100.0)**2.02)/(0.480*1000.0*(10**x)**3))},log10(40.0e-6),log10(2000.0e-6),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt,c("stellar crystal with broad arms (M96)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10((0.00583e-3*((10**x)*100.0)**2.42)/(0.480*1000.0*(10**x)**3))},log10(10.0e-6),log10(90.0e-6),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    plot(function(x){log10((0.000270e-3*((10**x)*100.0)**1.67)/(0.480*1000.0*(10**x)**3))},log10(90.0e-6),log10(1500.0e-6),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt,c("side plane (M96)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10((0.00419e-3*((10**x)*100.0)**2.3)/(0.480*1000.0*(10**x)**3))},log10(300.0e-6),log10(2500.0e-6),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt,c("dendritic crystal (HK87)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10((6.12e-7*((10**x)*100.0)**2.29)/(0.480*1000.0*(10**x)**3))},log10(0.6e-3),log10(5.3e-3),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt,c("dendritic crystal (K89)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10((4.82e-7*((10**x)*100.0)**1.97)/(0.480*1000.0*(10**x)**3))},log10(0.8e-3),log10(6.8e-3),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt,c("elementary needle (M90)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10((0.0049e-6*((10**x)*1000.0)**1.8)/(0.480*1000.0*(10**x)**3))},log10(0.6e-3),log10(2.7e-3),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt,c("hail (M96)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10((0.466*1000.0*(10**x)**3)/(0.480*1000.0*(10**x)**3))},log10(0.5e-2),log10(2.5e-2),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt,c("lump graupel (M96)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10((0.049e-3*((10**x)*100.0)**2.8)/(0.480*1000.0*(10**x)**3))},log10(0.05e-2),log10(0.3e-2),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt,c("lump graupel (HK87)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10((1.07e-4*((10**x)*100.0)**3.10)/(0.480*1000.0*(10**x)**3))},log10(0.4e-3),log10(9.0e-3),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt,c("densely rimed plate or sector (HK87)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10((9.53e-5*((10**x)*100.0)**3.8)/(0.480*1000.0*(10**x)**3))},log10(0.7e-3),log10(2.2e-3),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt,c("densely rimed dendrite (M96)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10((0.003e-3*((10**x)*100.0)**2.3)/(0.480*1000.0*(10**x)**3))},log10(0.18e-2),log10(0.4e-2),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt,c("densely rimed stellar (HK87)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10((7.55e-6*((10**x)*100.0)**3.04)/(0.480*1000.0*(10**x)**3))},log10(1.1e-3),log10(4.7e-3),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt,c("rimed long column (M96)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10((0.00145e-3*((10**x)*100.0)**1.8)/(0.480*1000.0*(10**x)**3))},log10(200.0e-6),log10(2400.0e-6),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt,c("rimed needle crystal (M90)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10((0.006e-6*((10**x)*1003.0)**2.1)/(0.480*1000.0*(10**x)**3))},log10(0.5e-3),log10(2.8e-3),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt,c("aggregate of side planes (LH74)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10((0.04e-6*((10**x)*1000.0)**1.4)/(0.480*1000.0*(10**x)**3))},log10(0.05e-2),log10(0.4e-2),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt,c("aggregate of side planes (M96)"))
    linecenter_col <- append(linecenter_col,c("black"))
    plot(function(x){log10((0.0033e-3*((10**x)*100.0)**2.2)/(0.480*1000.0*(10**x)**3))},log10(600.0e-6),log10(4100.0e-6),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt,c("aggregate of side planes, columns,bullets (M96)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10((0.0028e-3*((10**x)*100.0)**2.1)/(0.480*1000.0*(10**x)**3))},log10(0.08e-2),log10(0.45e-2),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt,c("aggregate of dendrites (LH74)"))
    linecenter_col <- append(linecenter_col,c("black"))
    plot(function(x){log10((0.037e-6*((10**x)*1000.0)**1.4)/(0.480*1000.0*(10**x)**3))},log10(2.0e-3),log10(12.0e-3),lty=linetype[n],col=alpha(linecenter_col[n],line_alpha),lwd=1,add=T)

    op <- par(cex = 0.8)
    legend("bottomleft",
    lty=1,col=colortype[1:n],lwd=3,legend=legendtxt[1:n],
    bty="n")

    legend("bottomleft",
    lty=linetype[1:n],col=linecenter_col[1:n],lwd=1,legend=legendtxt[1:n],
    bty="n")
    op <- par(cex = 1.0)

    # plot the mass density distribution of the "terminal velocity"-dimension relation
    dev.set(dev_vel)
    # plot fall velocities of typical ice particle types
    linetype <- rep(1,50)
    line_alpha <- 1.0
    rainbowalpha <- alpha(rainbow(50),line_alpha)
    colortype <- rep(1,50)

    n <- 0
    legendtxt <- character()

    n <- n+1
    legendtxt <- append(legendtxt, c("Stokes' law for ice sphere"))
    colortype[n] <- rainbowalpha[n+28]
    plot(function(x){log10((9.81*916.8/(18.0*1.527e-5))*(10.0**x)**2)},log10(0.5e-6),log10(20.0e-6),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt, c("column (SC85)"))
    colortype[n] <- rainbowalpha[n+28]
    plot(function(x){log10(8.114e-5*((10.0**x)*1.0e6)**1.585)},log10(1.0e-6),log10(200.0e-6),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt, c("hexagonal plate (L/(2a)=0.05) (W08)"))
    colortype[n] <- rainbowalpha[n+28]
    Westbrook2008_plate <- function(x){# x: maximum dimension in log10[m]
        D <- 10.0**x
	a <- D/2.0          # maximal radius of hexagon [m] 
	aspect_ratio <- 5.0e-2 # := L/2a
	L <- 2.0*a*aspect_ratio # height of the hexagonal prism [m]
    	grav <- 9.81     # [m/s^2]
	eta  <- 1.63e-5 # dynamics viscosity (T=-20C) [kg/m/s]
    	rho  <- 916.8    # [kg/m^3]
	mass <- rho*1.5*sqrt(3.0)*a**2*L # [kg]
	capaci <- 0.58*a*(1.0+0.95*aspect_ratio**0.75) # [m]
	eff_radi <- (4.0/3.0)*capaci # effective radius [m] (Roscoe horizontal orientation)
	v <- (grav/(6.0*pi*eta))*(mass/eff_radi) # terminal velocity [m/s]
	return(log10(v))
    }
    plot(Westbrook2008_plate,log10(1.0e-6),log10(200.0e-6),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt, c("hexagonal column (L/(2a)=20) (W08)"))
    colortype[n] <- rainbowalpha[n+28]
    Westbrook2008_column <- function(x){# x: maximum dimension in log10[m]
        D <- 10.0**x
	L <- D          # height of the hexagonal prism [m]
	aspect_ratio <- 0.2e2 # := L/2a
	a <- L/2.0/aspect_ratio # maximal radius of hexagon [m] 
    	grav <- 9.81     # [m/s^2]
	eta  <- 1.63e-5 # dynamics viscosity (T=-20C) [kg/m/s]
    	rho  <- 916.8    # [kg/m^3]
	mass <- rho*1.5*sqrt(3.0)*a**2*L # [kg]
	capaci <- 0.58*a*(1.0+0.95*aspect_ratio**0.75) # [m]
	eff_radi <- (1.0)*capaci # effective radius [m] (Hubbard-Douglas random orientation)
	v <- (grav/(6.0*pi*eta))*(mass/eff_radi) # terminal velocity [m/s]
	return(log10(v))
    }
    plot(Westbrook2008_column,log10(1.0e-6),log10(200.0e-6),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt, c("plate (H72)"))
    colortype[n] <- rainbowalpha[n+28]
    plot(function(x){log10((-0.55+90.69*((10.0**x)*1000.0)-23.44*((10.0**x)*1000.0)**2-5.26*((10.0**x)*1000.0)**3)/100.0)},log10(1.0e-4),log10(1.0e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt, c("dendrite (H72)"))
    colortype[n] <- rainbowalpha[n+28]
    plot(function(x){log10((25.9/100.0)*((10.0**x)*1000.0)**0.339)},log10(1e-3),log10(6e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt, c("column (T<-22C) (H72)"))
    colortype[n] <- rainbowalpha[n+28]
    plot(function(x){log10((-0.528+157.33*((10.0**x)*1000.0)-85.27*((10.0**x)*1000.0)**2)/100.0)},log10(0.1e-3),log10(1.0e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt, c("needle (H72)"))
    colortype[n] <- rainbowalpha[n+28]
    plot(function(x){log10((0.06+27.96*((10.0**x)*1000.0)-4.97*((10.0**x)*1000.0)**2+0.41*((10.0**x)*1000.0)**3)/100.0)},log10(1.0e-3),log10(5.0e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt, c("hexgonal plate (HK87)"))
    colortype[n] <- rainbowalpha[n+28]
    plot(function(x){log10(297.0e-2*((10.0**x)*100.0)**0.86)},log10(0.3e-3),log10(1.5e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt, c("stellar (HK87)"))
    colortype[n] <- rainbowalpha[n+28]
    plot(function(x){log10(58.0e-2*((10.0**x)*100.0)**0.55)},log10(0.4e-3),log10(2.4e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt, c("side plane (LH74)"))
    colortype[n] <- rainbowalpha[n+28]
    plot(function(x){log10((0.81)*((10.0**x)*1000.0)**0.99)},log10(0.4e-3),log10(1.2e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt,c("hail (M96)"))
    colortype[n] <- rainbowalpha[n-11]
    plot(function(x){log10((1184.0/100.0)*((10.0**x)*100.0)**0.5)},log10(5e-3),log10(25e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt, c("fresh hail (KH83)"))
    colortype[n] <- rainbowalpha[n-11]
    plot(function(x){log10(8.445*((10.0**x)*100.0)**0.553)},log10(0.5e-2),log10(2.0e-2),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt, c("hail and large graupel (A72)"))
    colortype[n] <- rainbowalpha[n-11]
    plot(function(x){log10(9.0*((10.0**x)*100.0)**0.8)},log10(0.1e-2),log10(8.0e-2),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt, c("lump graupel (LH74)"))
    colortype[n] <- rainbowalpha[n-11]
    plot(function(x){log10(1.3*((10.0**x)*1000.0)**0.66)},log10(0.5e-3),log10(3e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt, c("lump graupel (HK87)"))
    colortype[n] <- rainbowalpha[n-11]
    plot(function(x){log10(733.0e-2*((10.0**x)*100.0)**0.89)},log10(0.4e-3),log10(9.0e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

#    n <- n+1
#    legendtxt <- append(legendtxt, c("densely rimed plates?? (HK87)"))
#    plot(function(x){log10(92.0e-2*((10.0**x)*100.0)**0.73)},log10(0.7e-3),log10(2.2e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

#    n <- n+1
#    legendtxt <- append(legendtxt, c("densely rimed stellar?? (HK87)"))
#    plot(function(x){log10(162.0e-2*((10.0**x)*100.0)**0.53)},log10(1.1e-3),log10(4.7e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

#    n <- n+1
#    legendtxt <- append(legendtxt, c("rimed plates?? (HK87)"))
#    plot(function(x){log10(92.0e-2*((10.0**x)*100.0)**0.27)},log10(0.8e-3),log10(2.7e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

#    n <- n+1
#    legendtxt <- append(legendtxt, c("rimed stellar?? (HK87)"))
#    plot(function(x){log10(79.0e-2*((10.0**x)*100.0)**0.36)},log10(0.7e-3),log10(5.3e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt, c("densely rimed column (LH74)"))
    colortype[n] <- rainbowalpha[n-11]
    plot(function(x){log10(1.1*((10.0**x)*1000.0)**0.56)},log10(0.8e-3),log10(2.0e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt, c("densely rimed dendrite (LH74)"))
    colortype[n] <- rainbowalpha[n-11]
    plot(function(x){log10(0.62*((10.0**x)*1000.0)**0.33)},log10(1.8e-3),log10(4.0e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt, c("densely rimed radiating assemblag of dendrites (LH74)"))
    colortype[n] <- rainbowalpha[n-11]
    plot(function(x){log10(1.1*((10.0**x)*1000.0)**0.12)},log10(0.8e-3),log10(2.8e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt, c("aggregate of planar crystals (H02)"))
    colortype[n] <- rainbowalpha[2*(n-19)+14]
    plot(function(x){log10((136.0/100.0)*((10.0**x)*100.0)**0.4)},log10(6.0e-3),log10(36.0e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt, c("aggregate of unrimed dendrites (LH74)"))
    colortype[n] <- rainbowalpha[2*(n-19)+14]
    plot(function(x){log10((0.8)*((10.0**x)*1000.0)**0.16)},log10(2.0e-3),log10(10.0e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt, c("aggregate of densely rimed dendrites (LH74)"))
    colortype[n] <- rainbowalpha[2*(n-19)+14]
    plot(function(x){log10((0.79)*((10.0**x)*1000.0)**0.27)},log10(2.0e-3),log10(12.0e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt, c("aggregate of unrimed radiating assemblages of plates, side planes, bullets, and columns (LH74)"))
    colortype[n] <- rainbowalpha[2*(n-19)+14]
    plot(function(x){log10((0.69)*((10.0**x)*1000.0)**0.41)},log10(0.2e-3),log10(3.0e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    n <- n+1
    legendtxt <- append(legendtxt, c("aggregate of unrimed side planes (LH74)"))
    colortype[n] <- rainbowalpha[2*(n-19)+14]
    plot(function(x){log10((0.82)*((10.0**x)*1000.0)**0.12)},log10(0.5e-3),log10(4.0e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

    ## draw the center line
    linecenter_linetype <- seq(1,50)

    n <- 0
    linecenter_col <- character()

    n <- n+1
#    legendtxt <- append(legendtxt, c("Stokes' law for ice sphere"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10((9.81*916.8/(18.0*1.527e-5))*(10.0**x)**2)},log10(0.5e-6),log10(20.0e-6),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt, c("column (SC85)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10(8.114e-5*((10.0**x)*1.0e6)**1.585)},log10(1.0e-6),log10(200.0e-6),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt, c("hexagonal plate (L/(2a)=0.1) (W08)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(Westbrook2008_plate,log10(1.0e-6),log10(200.0e-6),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt, c("hexagonal column (L/(2a)=10) (W08)"))
    linecenter_col <- append(linecenter_col,c("black"))
    plot(Westbrook2008_column,log10(1.0e-6),log10(200.0e-6),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt, c("plate (H72)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10((-0.55+90.69*((10.0**x)*1000.0)-23.44*((10.0**x)*1000.0)**2-5.26*((10.0**x)*1000.0)**3)/100.0)},log10(1.0e-4),log10(1.0e-3),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt, c("dendrite (H72)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10((25.9/100.0)*((10.0**x)*1000.0)**0.339)},log10(1e-3),log10(6e-3),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt, c("column (T<-22C) (H72)"))
    linecenter_col <- append(linecenter_col,c("black"))
    plot(function(x){log10((-0.528+157.33*((10.0**x)*1000.0)-85.27*((10.0**x)*1000.0)**2)/100.0)},log10(0.1e-3),log10(1.0e-3),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt, c("needle (H72)"))
    linecenter_col <- append(linecenter_col,c("black"))
    plot(function(x){log10((0.06+27.96*((10.0**x)*1000.0)-4.97*((10.0**x)*1000.0)**2+0.41*((10.0**x)*1000.0)**3)/100.0)},log10(1.0e-3),log10(5.0e-3),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt, c("hexgonal plate (HK87)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10(297.0e-2*((10.0**x)*100.0)**0.86)},log10(0.3e-3),log10(1.5e-3),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt, c("stellar (HK87)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10(58.0e-2*((10.0**x)*100.0)**0.55)},log10(0.4e-3),log10(2.4e-3),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt, c("side plane (LH74)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10((0.81)*((10.0**x)*1000.0)**0.99)},log10(0.4e-3),log10(1.2e-3),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt,c("hail (M96)"))
    linecenter_col <- append(linecenter_col,c("black"))
    plot(function(x){log10((1184.0/100.0)*((10.0**x)*100.0)**0.5)},log10(5e-3),log10(25e-3),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt, c("fresh hail (KH83)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10(8.445*((10.0**x)*100.0)**0.553)},log10(0.5e-2),log10(2.0e-2),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt, c("hail and large graupel (A72)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10(9.0*((10.0**x)*100.0)**0.8)},log10(0.1e-2),log10(8.0e-2),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt, c("lump graupel (LH74)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10(1.3*((10.0**x)*1000.0)**0.66)},log10(0.5e-3),log10(3e-3),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt, c("lump graupel (HK87)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10(733.0e-2*((10.0**x)*100.0)**0.89)},log10(0.4e-3),log10(9.0e-3),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt, c("densely rimed column (LH74)"))
    linecenter_col <- append(linecenter_col,c("black"))
    plot(function(x){log10(1.1*((10.0**x)*1000.0)**0.56)},log10(0.8e-3),log10(2.0e-3),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt, c("densely rimed dendrite (LH74)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10(0.62*((10.0**x)*1000.0)**0.33)},log10(1.8e-3),log10(4.0e-3),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt, c("densely rimed radiating assemblage of dendrites (LH74)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10(1.1*((10.0**x)*1000.0)**0.12)},log10(0.8e-3),log10(2.8e-3),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt, c("aggregate of planar crystals (H02)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10((136.0/100.0)*((10.0**x)*100.0)**0.4)},log10(6.0e-3),log10(36.0e-3),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt, c("aggregate of unrimed dendrites (LH74)"))
    linecenter_col <- append(linecenter_col,c("black"))
    plot(function(x){log10((0.8)*((10.0**x)*1000.0)**0.16)},log10(2.0e-3),log10(10.0e-3),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt, c("aggregate of densely rimed dendrites (LH74)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10((0.79)*((10.0**x)*1000.0)**0.27)},log10(2.0e-3),log10(12.0e-3),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt, c("aggregate of unrimed radiating assemblages of plates, side planes, bullets, and columns (LH74)"))
    linecenter_col <- append(linecenter_col,c("black"))
    plot(function(x){log10((0.69)*((10.0**x)*1000.0)**0.41)},log10(0.2e-3),log10(3.0e-3),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)

    n <- n+1
#    legendtxt <- append(legendtxt, c("aggregate of unrimed side planes (LH74)"))
    linecenter_col <- append(linecenter_col,c("white"))
    plot(function(x){log10((0.82)*((10.0**x)*1000.0)**0.12)},log10(0.5e-3),log10(4.0e-3),lty=linecenter_linetype[n],col=linecenter_col[n],lwd=1,add=T)


    op <- par(cex = 0.6)
    legend("topleft",	
    lty=linetype[1:n],col=colortype[1:n],lwd=3,legend=legendtxt[1:n],
    bty="n")

    legend("topleft",
    lty=linecenter_linetype[1:n],col=linecenter_col[1:n],lwd=1,legend=legendtxt[1:n],
    bty="n")
    op <- par(cex = 1.0)

    graphics.off()

}
