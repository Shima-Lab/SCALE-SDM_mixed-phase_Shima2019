###########################################################################
###
### Program to calculate the time series of domain avaraged variables (0D)
###
############################################################################
###
### History: 171218, S.Shima, initital version
###          190113, S.Shima, also plot CIWP, GWP, SWP
###
############################################################################


############################################################################
########## Loading Libraries
library(ncdf4)

############################################################################
########## Parameters

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
########## LWP,CWP,RWP
cat(sprintf("start plotting water path\n"))

TOT_CLOUD_MASS <- rep(0.0,length(alltimes))
TOT_RAIN_MASS  <- rep(0.0,length(alltimes))
names(TOT_CLOUD_MASS) <- alltimes
names(TOT_RAIN_MASS)  <- alltimes
AREA <- 0.0

##### loop of MPI rank
for(mpirank in allmpiranks){
	cat(sprintf("processing the rank = %s \n",mpirank))

	##### Open
	ncin <- nc_open(paste("../",allfiles[mpirank],sep=""))

	##### Grids
	IMAX <- length(ncvar_get(ncin,"CX"))-4
	JMAX <- length(ncvar_get(ncin,"CY"))-4
	KMAX <- length(ncvar_get(ncin,"CZ"))-4
	CDX <- ncvar_get(ncin,"CDX")
	CDY <- ncvar_get(ncin,"CDY")
	CDZ <- ncvar_get(ncin,"CDZ")

	##### Read DENS
	DENS <- ncvar_get(ncin,"DENS")
	DENS_units <- ncatt_get(ncin,"DENS","units")
	dimnames(DENS)[[4]] <- alltimes
	if(!setequal(dim(DENS[,,,1]),c(IMAX,JMAX,KMAX))){
		cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
		return()
	}

	##### Read QC
	QC <- ncvar_get(ncin,"QC_sd")
	QC_units <- ncatt_get(ncin,"QC_sd","units")
	dimnames(QC)[[4]] <- alltimes

	##### Read QR
	QR <- ncvar_get(ncin,"QR_sd")
	QR_units <- ncatt_get(ncin,"QR_sd","units")
	dimnames(QR)[[4]] <- alltimes

	##### Close
	nc_close(ncin)

	##### loop of time
	for(time in alltimes){
		CLOUD_MASS <- QC[1:IMAX,1:JMAX,1:KMAX,as.character(time)]*
                              DENS[1:IMAX,1:JMAX,1:KMAX,as.character(time)]*
			      (CDX[3:(IMAX+2)]%o%CDY[3:(JMAX+2)]%o%CDZ[3:(KMAX+2)]) # volume of each grid box
		TOT_CLOUD_MASS[as.character(time)] <- TOT_CLOUD_MASS[as.character(time)] + 
						      sum(CLOUD_MASS[1:IMAX,1:JMAX,1:KMAX])

		RAIN_MASS  <- QR[1:IMAX,1:JMAX,1:KMAX,as.character(time)]*
                              DENS[1:IMAX,1:JMAX,1:KMAX,as.character(time)]*
			      (CDX[3:(IMAX+2)]%o%CDY[3:(JMAX+2)]%o%CDZ[3:(KMAX+2)])
		TOT_RAIN_MASS[as.character(time)]  <- TOT_RAIN_MASS[as.character(time)] + 
						      sum(RAIN_MASS[1:IMAX,1:JMAX,1:KMAX])
	}

	##### Horizontal area 
	AREA <- AREA + sum(CDX[3:(IMAX+2)])*sum(CDY[3:(JMAX+2)])
}

CWP_AVE <- TOT_CLOUD_MASS/AREA
RWP_AVE <- TOT_RAIN_MASS/AREA
LWP_AVE <- CWP_AVE + RWP_AVE

############################################################################
########## IWP
TOT_CLOUDICE_MASS <- rep(0.0,length(alltimes))
TOT_SNOW_MASS  <- rep(0.0,length(alltimes))
TOT_GRAUPEL_MASS  <- rep(0.0,length(alltimes))
names(TOT_CLOUDICE_MASS) <- alltimes
names(TOT_SNOW_MASS) <- alltimes
names(TOT_GRAUPEL_MASS) <- alltimes
AREA <- 0.0

##### loop of MPI rank
for(mpirank in allmpiranks){
	cat(sprintf("processing the rank = %s \n",mpirank))

	##### Open
	ncin <- nc_open(paste("../",allfiles[mpirank],sep=""))

	##### Grids
	IMAX <- length(ncvar_get(ncin,"CX"))-4
	JMAX <- length(ncvar_get(ncin,"CY"))-4
	KMAX <- length(ncvar_get(ncin,"CZ"))-4
	CDX <- ncvar_get(ncin,"CDX")
	CDY <- ncvar_get(ncin,"CDY")
	CDZ <- ncvar_get(ncin,"CDZ")

	##### Read DENS
	DENS <- ncvar_get(ncin,"DENS")
	DENS_units <- ncatt_get(ncin,"DENS","units")
	dimnames(DENS)[[4]] <- alltimes
	if(!setequal(dim(DENS[,,,1]),c(IMAX,JMAX,KMAX))){
		cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
		return()
	}

	##### Read QI
	QI <- ncvar_get(ncin,"QI_sd")
	QI_units <- ncatt_get(ncin,"QI_sd","units")
	dimnames(QI)[[4]] <- alltimes

	##### Read QS
	QS <- ncvar_get(ncin,"QS_sd")
	QS_units <- ncatt_get(ncin,"QS_sd","units")
	dimnames(QS)[[4]] <- alltimes

	##### Read QG
	QG <- ncvar_get(ncin,"QG_sd")
	QG_units <- ncatt_get(ncin,"QG_sd","units")
	dimnames(QG)[[4]] <- alltimes

	##### Close
	nc_close(ncin)

	##### loop of time
	for(time in alltimes){
		CLOUDICE_MASS <- QI[1:IMAX,1:JMAX,1:KMAX,as.character(time)]*
                              DENS[1:IMAX,1:JMAX,1:KMAX,as.character(time)]*
			      (CDX[3:(IMAX+2)]%o%CDY[3:(JMAX+2)]%o%CDZ[3:(KMAX+2)])
		TOT_CLOUDICE_MASS[as.character(time)] <- TOT_CLOUDICE_MASS[as.character(time)] + 
						      sum(CLOUDICE_MASS[1:IMAX,1:JMAX,1:KMAX])

		SNOW_MASS  <- QS[1:IMAX,1:JMAX,1:KMAX,as.character(time)]*
                              DENS[1:IMAX,1:JMAX,1:KMAX,as.character(time)]*
			      (CDX[3:(IMAX+2)]%o%CDY[3:(JMAX+2)]%o%CDZ[3:(KMAX+2)])
		TOT_SNOW_MASS[as.character(time)]  <- TOT_SNOW_MASS[as.character(time)] + 
						      sum(SNOW_MASS[1:IMAX,1:JMAX,1:KMAX])

		GRAUPEL_MASS  <- QG[1:IMAX,1:JMAX,1:KMAX,as.character(time)]*
                              DENS[1:IMAX,1:JMAX,1:KMAX,as.character(time)]*
			      (CDX[3:(IMAX+2)]%o%CDY[3:(JMAX+2)]%o%CDZ[3:(KMAX+2)])
		TOT_GRAUPEL_MASS[as.character(time)]  <- TOT_GRAUPEL_MASS[as.character(time)] + 
						      sum(GRAUPEL_MASS[1:IMAX,1:JMAX,1:KMAX])
	}

	##### Horizontal area 
	AREA <- AREA + sum(CDX[3:(IMAX+2)])*sum(CDY[3:(JMAX+2)])
}

CIWP_AVE <- TOT_CLOUDICE_MASS/AREA
SWP_AVE <- TOT_SNOW_MASS/AREA
GWP_AVE <- TOT_GRAUPEL_MASS/AREA
IWP_AVE <- CIWP_AVE + SWP_AVE + GWP_AVE

############################################################################
########## Plot LWP,CWP,RWP,IWP,CIWP,GWP,SWP
plot_var_list <- c(1,2,3,4,5,6,7)

var_name <- c("LWP_AVE","CWP_AVE","RWP_AVE","IWP_AVE","CIWP_AVE","GWP_AVE","SWP_AVE")
cols   <- c("black","dark blue", "dark orange3", "dark green", "blue", "red", "green")
#ltys   <- c(1, 2, 4, 5, 6, 7, 8)
ltys   <- c(1, 1, 1, 1, 1, 1, 1)
lwds <- c(1, 1, 1, 2, 2, 2, 2)
labels <- c("LWP","CWP","RWP","IWP","CIWP","GWP","SWP")
ymax <- 1.1

pdf("water_path.pdf")

# plot the axis
i <- 1
plot(alltimes,get(var_name[i]),
	type="n",lty = ltys[i], col = cols[i], lwd = lwds[i],
	ylim=c(0.0,ymax),
	main="Domain Averaged Water Path",
	xlab="Time [s]",
	ylab="Domain Averaged Water Path [kg/m^2]",
	ann=T)
par(new=T)

# plot the variables
for(i in plot_var_list){
	plot(alltimes,get(var_name[i]),
		type="l",lty = ltys[i], col = cols[i], lwd = lwds[i],
		ylim=c(0.0,ymax),
		main="Domain Averaged Water Path",
		xlab="Time [s]",
		ylab="Domain Averaged Water Path [kg/m^2]",
		ann=F)
	par(new=T)
}

par(new=F)

legend("topleft", legend = labels[plot_var_list], col = cols[plot_var_list],
lty = ltys[plot_var_list], lwd = lwds[plot_var_list])

dev.off()

sink("max_water_path.dat")
for(i in plot_var_list){
      cat(sprintf("%s %e\n",var_name[i],max(get(var_name[i]))))
}
sink()

sink("water_path.dat")
for(i in plot_var_list){
      for(j in 1:length(alltimes)){
	    cat(sprintf("%s %s %e\n", var_name[i], alltimes[j], get(var_name[i])[j]))
      }
}
sink()

cat(sprintf("end plotting water path\n\n"))

############################################################################
########## Domain Averaged Accumulated Precipitation Amount
cat(sprintf("start plotting accumulated precipitation amount\n"))

TOT_ACC_RAIN <- rep(0.0,length(alltimes))
names(TOT_ACC_RAIN) <- alltimes
AREA <- 0.0

##### loop of MPI rank
for(mpirank in allmpiranks){
	cat(sprintf("processing the rank = %s \n",mpirank))

	##### Open
	ncin <- nc_open(paste("../",allfiles[mpirank],sep=""))

	##### Grids
	IMAX <- length(ncvar_get(ncin,"CX"))-4
	JMAX <- length(ncvar_get(ncin,"CY"))-4
	KMAX <- length(ncvar_get(ncin,"CZ"))-4
	CDX <- ncvar_get(ncin,"CDX")
	CDY <- ncvar_get(ncin,"CDY")
	CDZ <- ncvar_get(ncin,"CDZ")

	##### Read RAIN_sd
	RAIN <- ncvar_get(ncin,"RAIN_sd")
	RAIN_units <- ncatt_get(ncin,"RAIN_sd","units")
	dimnames(RAIN)[[3]] <- alltimes
	if(!setequal(dim(RAIN[,,1]),c(IMAX,JMAX))){
		cat(sprintf("ERROR: Size of RAIN is not consistent. Check HALO treatment,\n"))
		return()
	}

	##### Close
	nc_close(ncin)

	##### loop of time
	for(time in alltimes){
		ACC_RAIN <- RAIN[1:IMAX,1:JMAX,as.character(time)]*(CDX[3:(IMAX+2)]%o%CDY[3:(JMAX+2)])
		TOT_ACC_RAIN[as.character(time)] <- TOT_ACC_RAIN[as.character(time)] + sum(ACC_RAIN[1:IMAX,1:JMAX])
	}

	##### Horizontal area 
	AREA <- AREA + sum(CDX[3:(IMAX+2)])*sum(CDY[3:(JMAX+2)])
}

ACC_RAIN_AVE <- TOT_ACC_RAIN/AREA

############################################################################
########## Precipitation Amount
ymax <- 1.5e-3

pdf("precipitation.pdf")

plot(alltimes,ACC_RAIN_AVE,
	type="o",
	ylim=c(0.0,ymax),
	main="Domain Averaged Precipitation Amount [m]",
        xlab="Time [s]",
        ylab="Domain Averaged Precipitation Amount [m]")

dev.off()

tmp <- read.table(textConnection(gsub("="," ",gsub(",", " ", readLines("../run.conf")))),header=FALSE,fill=TRUE)
sdm_inisdnc <- as.numeric(gsub("[dD]","e",tmp[tmp$V1=="sdm_inisdnc",2]))

sink("precipitation.dat")
	for(time in alltimes){
		 cat(sprintf("%s %e %e\n",time,ACC_RAIN_AVE[as.character(time)],sdm_inisdnc))
	}
sink()

cat(sprintf("end plotting accumulated precipitation amount\n\n"))
