###########################################################################
###
### Program to plot the x-z distribution of hydrometeor mixing ratio
###
############################################################################

############################################################################
########## Loading Libraries
library(ncdf4)

############################################################################
########## Parameters
pal <- c(colorRampPalette(c(rgb(     0,   0,      1, 0.0), rgb(     0,   0,      1,0.6)), alpha = TRUE),
         colorRampPalette(c(rgb(     0,   1,      0, 0.0), rgb(     0,   1,      0,0.6)), alpha = TRUE),
         colorRampPalette(c(rgb(     1,   0,      0, 0.0), rgb(     1,   0,      0,0.6)), alpha = TRUE),
         colorRampPalette(c(rgb(     1,   1,      1, 0.0), rgb(     1,   1,      1,0.6)), alpha = TRUE),
         colorRampPalette(c(rgb(     1,   1,      0, 0.0), rgb(     1,   1,      0,0.6)), alpha = TRUE)
	 )
names(pal)=c("ice","snow","graupel","cloud water","rain")

zlim_min <- 1e-7
zlim_max <- 1e-2
xlim_min <- 0.0e3  # horizontal range [m]
xlim_max <- 60.0e3 # [m]
ylim_min <- 0.0e3  # vertical range [m]
ylim_max <- 16.0e3 # [m]
w_plot_area <- (xlim_max-xlim_min)/1.0e4
h_plot_area <- (ylim_max-ylim_min)/1.0e4

############################################################################
########## Read Data from Files
##### Make a list of files,  mpiranks
allfiles = dir("../",pattern="^history.")
tmp = strsplit(allfiles,"\\history.pe|\\.nc")
allmpiranks = unique(matrix(unlist(tmp),nrow=2)[2,])
names(allfiles) = allmpiranks

##### Open the first file
ncin <- nc_open(paste("../",allfiles[1],sep=""))

##### Times
alltimes <- ncvar_get(ncin,"time")
time_units <- ncatt_get(ncin,"time","units")

##### Close
nc_close(ncin)

############################################################################
########## Prepare a Vector for Devices
dev_num <- rep(0,length(alltimes))
names(dev_num) <- alltimes

############################################################################
########## Plot 50 files at once
skip <- as.integer(length(alltimes)/50)+1
for(i0 in 1:skip){
    #### Prepare Devices (=Files)
    for(time in alltimes[seq(i0,length(alltimes),by=skip)]){
	cat(sprintf("processing time = %s [s]	\n",time))

	pdf(paste("QHYD_overlay.",sprintf("%05d",as.numeric(time)),".pdf",sep=""),width=w_plot_area*2+2,height=h_plot_area*2+2)
    	dev_num[as.character(time)] <- dev.cur()
	par(pin=c(w_plot_area*2,h_plot_area*2))
	par(bg="black", col="white", col.axis="white", 
	       col.lab="white", fg="white", col.main="white", col.sub="white")
	image(c(0,1), c(0,1), diag(2)*0,
                   xlim=c(xlim_min,xlim_max),
                   ylim=c(ylim_min,ylim_max),
		   asp = 1,
                   zlim=c(log(zlim_min),log(zlim_max)),
                   col = rgb(1,1,1,0),
                   main=paste("Mixing Ratio of Hydrometeors (T=",
                              sprintf("%05d",as.numeric(time)),
                              "s)"),
                   xlab="x [m]",
                   ylab="z [m]")
    }

    ##### loop of MPI rank
    for(mpirank in allmpiranks){
	cat(sprintf("processing the rank = %s \n",mpirank))

	##### Open
	ncin <- nc_open(paste("../",allfiles[mpirank],sep=""))

	##### Grids
	IMAX <- length(ncvar_get(ncin,"CX"))-4
	JMAX <- length(ncvar_get(ncin,"CY"))-4
	KMAX <- length(ncvar_get(ncin,"CZ"))-4
	CX <- ncvar_get(ncin,"CX")
	CY <- ncvar_get(ncin,"CY")
	CZ <- ncvar_get(ncin,"CZ")


	##### Read QR
	QR <- ncvar_get(ncin,"QR_sd")
	QR_units <- ncatt_get(ncin,"QR_sd","units")
	dimnames(QR)[[4]] <- alltimes

	##### Read QC
	QC <- ncvar_get(ncin,"QC_sd")
	QC_units <- ncatt_get(ncin,"QC_sd","units")
	dimnames(QC)[[4]] <- alltimes

	##### Read QG
	QG <- ncvar_get(ncin,"QG_sd")
	QG_units <- ncatt_get(ncin,"QG_sd","units")
	dimnames(QG)[[4]] <- alltimes

	##### Read QI
	QI <- ncvar_get(ncin,"QI_sd")
	QI_units <- ncatt_get(ncin,"QI_sd","units")
	dimnames(QI)[[4]] <- alltimes

	##### Read QS
	QS <- ncvar_get(ncin,"QS_sd")
	QS_units <- ncatt_get(ncin,"QS_sd","units")
	dimnames(QS)[[4]] <- alltimes
	
	##### Close
	nc_close(ncin)

	#### Plot Data
#	for(time in alltimes[which(as.numeric(alltimes)==3600)]){
	for(time in alltimes[seq(i0,length(alltimes),by=skip)]){
		 #### set the device number
		 dev.set(dev_num[as.character(time)])

		 for(hyd_category in c("graupel", "rain", "cloud water", "ice", "snow")){
		 	if(hyd_category=="graupel"){
				makeActiveBinding("qhyd", function() QG,.GlobalEnv)
			}else if(hyd_category=="ice"){
				makeActiveBinding("qhyd", function() QI,.GlobalEnv)
			}else if(hyd_category=="snow"){
				makeActiveBinding("qhyd", function() QS,.GlobalEnv)
			}else if(hyd_category=="rain"){
				makeActiveBinding("qhyd", function() QR,.GlobalEnv)
			}else if(hyd_category=="cloud water"){
				makeActiveBinding("qhyd", function() QC,.GlobalEnv)
			}

		 	par(pin=c(w_plot_area*2,h_plot_area*2))
		 	par(new=T)
		 	image(CX[3:(IMAX+2)], CZ[3:(KMAX+2)],
		   		log(pmax(pmin(qhyd[1:IMAX,1:1,1:KMAX,as.character(time)],zlim_max),zlim_min)),
                   		xlim=c(xlim_min,xlim_max),
                   		ylim=c(ylim_min,ylim_max),
		   		asp = 1,
				zlim=c(log(zlim_min),log(zlim_max)),
                   		col=pal[[hyd_category]](100),
                   		main=paste("Mixing Ratio of Hydrometeors (T=",
                              			   sprintf("%05d",as.numeric(time)),
                              			   "s)"),
                   		xlab="x [m]",
                   		ylab="z [m]",
				ann=F,
				useRaster=TRUE)
		}
	}
    }

    #### Close all the files
    graphics.off()

}
