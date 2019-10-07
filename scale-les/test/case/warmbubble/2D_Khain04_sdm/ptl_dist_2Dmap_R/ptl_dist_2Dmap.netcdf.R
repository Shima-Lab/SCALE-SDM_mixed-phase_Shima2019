## This is a R script to plot 2D histogram of ice particles
# Load required libraries
library(fields)
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

# set parameters
VALID2INVALID <- -999.0

# Setting of bins
nbins = 100
min_radius = 1.0e-7 # [m]
max_radius = 2.0e-1 # [m]
min_maxdim = 2.0*min_radius
max_maxdim = 2.0*max_radius
min_density= 1.0e0  # [kg/m^3]
max_density= 1.1e+3 # [kg/m^3]
min_mass = 1.0e-16  # [kg]
max_mass = 1.0e-1   # [kg]
min_vel  = 1.0e-3   # [m/s]
max_vel  = 20.0e0   # [m/s]

equr.bin = seq(log10(min_radius), log10(max_radius), length=nbins)
polr.bin = seq(log10(min_radius), log10(max_radius), length=nbins)
maxr.bin = seq(log10(min_radius), log10(max_radius), length=nbins)
maxdim.bin  = seq(log10(min_maxdim), log10(max_maxdim), length=nbins)
dens.bin = seq(log10(min_density),log10(max_density),length=nbins)
mass.bin = seq(log10(min_mass),log10(max_mass),length=nbins)
vel.bin = seq(log10(min_vel),log10(max_vel),length=nbins)

d_equr.bin = equr.bin[2] - equr.bin[1]
d_polr.bin = polr.bin[2] - polr.bin[1]
d_maxr.bin = maxr.bin[2] - maxr.bin[1]
d_maxdim.bin  = maxdim.bin[2]  - maxdim.bin[1]
d_dens.bin = dens.bin[2] - dens.bin[1]
d_mass.bin = mass.bin[2] - mass.bin[1]
d_vel.bin = vel.bin[2] - vel.bin[1]

# Make the list of files
tmp_files = dir("../",pattern="^SD_output_NetCDF_")
tmp = strsplit(tmp_files,"\\SD_output_NetCDF_|\\.000.pe")
alltimes = unique(matrix(unlist(tmp),nrow=3)[2,])
allmpiranks = unique(matrix(unlist(tmp),nrow=3)[3,])

allfiles = matrix(tmp_files,ncol=length(alltimes))
rownames(allfiles) = allmpiranks
colnames(allfiles) = alltimes

# loop of time
for(time in alltimes){
#for(time in alltimes[which(as.numeric(alltimes)==3600)]){
	 cat(sprintf("processing the time = %s [sec]\n",time))

	 aspect_dist <- diag(nbins)*0.0
	 density_dist <- diag(nbins)*0.0
	 mass_dist <- diag(nbins)*0.0
	 vel_dist <- diag(nbins)*0.0

	 # loop of MPI rank
	 for(mpirank in allmpiranks){
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

		equr <- equr[typed_ice_id]
		polr <- polr[typed_ice_id]
		dens <- dens[typed_ice_id]
		mult <- mult[typed_ice_id]
		icetvel <- icetvel[typed_ice_id]

		# calculate the mass density distribution of aspect ratio
		ptl_num <- length(mult)
		iequr <- numeric(ptl_num)
		ipolr <- numeric(ptl_num)
		iequr[1:ptl_num] <- findInterval(log10(equr[1:ptl_num]), equr.bin, all.inside=TRUE)
		ipolr[1:ptl_num] <- findInterval(log10(polr[1:ptl_num]), polr.bin, all.inside=TRUE)
		for (n in 1:ptl_num){
    		    aspect_dist[iequr[n],ipolr[n]] <- aspect_dist[iequr[n],ipolr[n]] +
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

		# calculate the mass density distribution of mass-dimension
		ptl_num <- length(mult)
		maxdim <- numeric(ptl_num)
		mass <- numeric(ptl_num)
		imaxdim <- numeric(ptl_num)
		imass <- numeric(ptl_num)
		maxdim[1:ptl_num] <- 2.0*pmax(equr[1:ptl_num],polr[1:ptl_num])
		mass[1:ptl_num] <- dens[1:ptl_num]*(4.0/3.0)*pi*equr[1:ptl_num]**2*polr[1:ptl_num]
		imaxdim[1:ptl_num] <- findInterval(log10(maxdim[1:ptl_num]), maxdim.bin, all.inside=TRUE)
		imass[1:ptl_num] <- findInterval(log10(mass[1:ptl_num]), mass.bin, all.inside=TRUE)
		for (n in 1:ptl_num){
    		    mass_dist[imaxdim[n],imass[n]] <- mass_dist[imaxdim[n],imass[n]] +
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
	aspect_dist[,] <- aspect_dist[,]/d_equr.bin/d_polr.bin
	density_dist[,] <- density_dist[,]/d_maxdim.bin/d_dens.bin
	mass_dist[,] <- mass_dist[,]/d_maxdim.bin/d_mass.bin
	vel_dist[,] <- vel_dist[,]/d_maxdim.bin/d_vel.bin

	# plot the equatorial and polar radius distribution
	pdf(paste("aspect.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
	image.plot(equr.bin, polr.bin, log10(aspect_dist), asp=1,
		   zlim=c(-3.0,10.0), nlevel=100,
		   main=paste("Aspect Ratio Distribution (T=",
		              sprintf("%05d",as.numeric(time)),
		              "s)\n(Mass Density log([kg/unit_log10(eq_r)/unit_log10(pol_r)]))"),
		   xlab="Equatorial radius log10([m])",
		   ylab="Polar radius log10([m])")
	dev.off()

	# plot the mass density distribution of size-density relation
	pdf(paste("density.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
	image.plot(maxdim.bin, dens.bin, log10(density_dist),
		   zlim=c(-3.0,10.0), nlevel=100,
		   main=paste("Size-Density Distribution (T=",
		              sprintf("%05d",as.numeric(time)),
		              "s)\n(Mass Density log([kg/unit_log10(max_D)/unit_log10(density)]))"),
		   xlab="Size (maximum dimension) log10([m])",
		   ylab="Density of ice particle log10([kg/m^3])")
	dev.off()

	# plot the mass density distribution of mass-dimension relation
	pdf(paste("mass.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
	image.plot(maxdim.bin, mass.bin, log10(mass_dist),
		   zlim=c(-3.0,10.0), nlevel=100,
		   main=paste("Mass-Dimension Distribution (T=",
		              sprintf("%05d",as.numeric(time)),
                              "s)\n(Mass Density log([kg/unit_log10(max_D)/unit_log10(density)]))"),
		   xlab="Size (maximum dimension) log10([m])",
		   ylab="Mass of ice particle log10([kg])")

	# plot various typical mass-dimension relations
	linetype <- rep(1,20)  # all solid lines
	colortype <- seq(2,21) # skip 1 (black)
	n <- 0 
	n <- n+1
	legendtxt <- c("ice spheres")
        plot(function(x){log10(0.4765*1000.0*(10**x)**3)},log10(min_maxdim),log10(max_maxdim),lty=linetype[n],col=colortype[n],lwd=2,add=T)

	n <- n+1
	legendtxt <- append(legendtxt,c("hail (M96)"))
        plot(function(x){log10(0.466*1000.0*(10**x)**3)},log10(0.5e-2),log10(2.5e-2),lty=linetype[n],col=colortype[n],lwd=2,add=T)

	n <- n+1
	legendtxt <- append(legendtxt,c("lump graupel (M96)"))
        plot(function(x){log10(0.049e-3*((10**x)*100.0)**2.8)},log10(0.05e-2),log10(0.3e-2),lty=linetype[n],col=colortype[n],lwd=2,add=T)

	n <- n+1
	legendtxt <- append(legendtxt,c("densely rimed dendrites (M96)"))
        plot(function(x){log10(0.003e-3*((10**x)*100.0)**2.3)},log10(0.18e-2),log10(0.4e-2),lty=linetype[n],col=colortype[n],lwd=2,add=T)

#	n <- n+1
#	legendtxt <- append(legendtxt,c("aggregates (SK10)"))
#        plot(function(x){log10(0.00183e-3*((10**x)*100.0)**2.04)},log10(min_maxdim),log10(max_maxdim),lty=linetype[n],col=colortype[n],lwd=2,add=T)

#	n <- n+1
#	legendtxt <- append(legendtxt,c("aggregates of unrimed side planes (LH74)"))
#        plot(function(x){log10(0.04e-6*((10**x)*1000.0)**1.4)},log10(0.05e-2),log10(0.4e-2),lty=linetype[n],col=colortype[n],lwd=2,add=T)

	n <- n+1
	legendtxt <- append(legendtxt,c("aggregates (M96)"))
        plot(function(x){log10(0.0028e-3*((10**x)*100.0)**2.1)},log10(0.08e-2),log10(0.45e-2),lty=linetype[n],col=colortype[n],lwd=2,add=T)

	n <- n+1
	legendtxt <- append(legendtxt,c("hexagonal columns (small) (M96)"))
        plot(function(x){log10(0.1677e-3*((10**x)*100.0)**2.91)},log10(0.003e-2),log10(0.01e-2),lty=linetype[n],col=colortype[n],lwd=2,add=T)

	n <- n+1
	legendtxt <- append(legendtxt,c("hexagonal columns (middle) (M96)"))
        plot(function(x){log10(0.00166e-3*((10**x)*100.0)**1.91)},log10(0.01e-2),log10(0.03e-2),lty=linetype[n],col=colortype[n],lwd=2,add=T)

	n <- n+1
	legendtxt <- append(legendtxt,c("hexagonal columns (big) (M96)"))
        plot(function(x){log10(0.000907e-3*((10**x)*100.0)**1.74)},log10(0.03e-2),log10(0.3e-2),lty=linetype[n],col=colortype[n],lwd=2,add=T)

	n <- n+1
	legendtxt <- append(legendtxt,c("hexagonal plates (M96)"))
        plot(function(x){log10(0.00739e-3*((10**x)*100.0)**2.45)},log10(0.0015e-2),log10(0.3e-2),lty=linetype[n],col=colortype[n],lwd=2,add=T)

	n <- n+1
	legendtxt <- append(legendtxt,c("broad-branched crystals (small) (M96)"))
        plot(function(x){log10(0.00583e-3*((10**x)*100.0)**2.42)},log10(0.001e-2),log10(0.01e-2),lty=linetype[n],col=colortype[n],lwd=2,add=T)

	n <- n+1
	legendtxt <- append(legendtxt,c("broad-branched crystals (big) (M96)"))
        plot(function(x){log10(0.000516e-3*((10**x)*100.0)**1.8)},log10(0.01e-2),log10(0.1e-2),lty=linetype[n],col=colortype[n],lwd=2,add=T)

	n <- n+1
	legendtxt <- append(legendtxt,c("vapor spheres"))
        plot(function(x){log10(0.8*(pi/6.0)*(10**x)**3)},log10(min_maxdim),log10(max_maxdim),lty=linetype[n],col=colortype[n],lwd=2,add=T)

	op <- par(cex = 0.8)
	legend("topleft",
	lty=linetype[1:n],col=colortype[1:n],legend=legendtxt[1:n],
	bty="n")

	dev.off()

	# plot the mass density distribution of the "terminal velocity"-dimension relation
	pdf(paste("term_vel_ice.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
	image.plot(maxdim.bin, vel.bin, log10(vel_dist),
		   zlim=c(-3.0,10.0), nlevel=100,
		   main=paste("Terminal velocity of ice particles (T=",
				     sprintf("%05d",as.numeric(time)),
				     "s)\n(Mass Density log([kg/unit_log10(max_D)/unit_log10(velocity)]))"),
		   xlab="Dimension of ice particles log10([m])",
		   ylab="Terminal velocity of ice particles log10([m/s])")

        # plot fall velocities of typical ice particle types
        linetype <- rep(1,10)  # all solid lines
        colortype <- seq(2,11) # skip 1 (black)
        n <- 0

        n <- n+1
        legendtxt <- c("hail (M96)")
        plot(function(x){log10((1184.0/100.0)*((10.0**x)*100.0)**0.5)},log10(5e-3),log10(25e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

        n <- n+1
        legendtxt <- append(legendtxt, c("lump graupel (LH74)"))
        plot(function(x){log10(1.3*((10.0**x)*1000.0)**0.66)},log10(0.5e-3),log10(3e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

        n <- n+1
        legendtxt <- append(legendtxt, c("side planes (LH74)"))
        plot(function(x){log10((0.81)*((10.0**x)*1000.0)**0.99)},log10(0.4e-3),log10(1.2e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

        n <- n+1
        legendtxt <- append(legendtxt, c("plate (H72)"))
        plot(function(x){log10((-0.55+90.69*((10.0**x)*1000.0)-23.44*((10.0**x)*1000.0)**2-5.26*((10.0**x)*1000.0)**3)/100.0)},log10(1.0e-5),log10(1.0e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

        n <- n+1
        legendtxt <- append(legendtxt, c("aggregates of\nplanar crystals (H02)"))
        plot(function(x){log10((136.0/100.0)*((10.0**x)*100.0)**0.4)},log10(6.0e-3),log10(36.0e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

        n <- n+1
        legendtxt <- append(legendtxt, c("dendrite (H72)"))
        plot(function(x){log10((25.9/100.0)*((10.0**x)*1000.0)**0.339)},log10(1e-3),log10(6e-3),lty=linetype[n],col=colortype[n],lwd=3,add=T)

        legend("topleft",
        lty=linetype[1:n],col=colortype[1:n],legend=legendtxt[1:n],cex=0.8)

        dev.off()

}
