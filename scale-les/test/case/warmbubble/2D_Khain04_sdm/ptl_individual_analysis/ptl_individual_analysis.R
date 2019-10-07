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

# Make the list of files
tmp_files = dir("../",pattern="^SD_output_NetCDF_")
tmp = strsplit(tmp_files,"\\SD_output_NetCDF_|\\.000.pe")
alltimes = unique(matrix(unlist(tmp),nrow=3)[2,])
allmpiranks = unique(matrix(unlist(tmp),nrow=3)[3,])

allfiles = matrix(tmp_files,ncol=length(alltimes))
rownames(allfiles) = allmpiranks
colnames(allfiles) = alltimes

# loop of time
#for(time in alltimes){
for(time in alltimes[which(as.numeric(alltimes)==3000)]){
    cat(sprintf("processing the time = %s [sec]\n",time))

    sink(paste("selected_ice.",sprintf("%05d",as.numeric(time)),".dat",sep=""))
    sink()
    sink(paste("selected_ice.",sprintf("%05d",as.numeric(time)),".dat",sep=""))
    cat(sprintf("time[s] x[m] z[m] a[m] c[m] dens[kg/m^3] mult[] tvel[m/s] mrime[kg] nmono[] ice_category\n"))

#    for(ice_category in c("snow")){
    for(ice_category in c("graupel", "ice", "snow")){

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

		# x coordinate
		all_sd_x <- ncvar_get(ncin,"sd_x")
		sd_x     <- all_sd_x[ice_id]
		# z coordinate
		sd_z     <- all_sd_z[ice_id]

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

		# save categorized ice particles
		equr <- equr[typed_ice_id]
		polr <- polr[typed_ice_id]
		dens <- dens[typed_ice_id]
		mult <- mult[typed_ice_id]
		icetvel <- icetvel[typed_ice_id]
		mrime <- mrime[typed_ice_id]
		nmono <- nmono[typed_ice_id]
		sd_x <- sd_x[typed_ice_id]
		sd_z <- sd_z[typed_ice_id]

		### output selected ice		
		selected_ice_id <- which( ((polr/equr)>10.0) & (2.0*polr>0.03) )
		if(length(selected_ice_id)>0){
			for(id in selected_ice_id){
    		    	       cat(sprintf("%d %e %e %e %e %e %e %e %e %d %s\n",as.integer(time),sd_x[id], sd_z[id], equr[id],polr[id],dens[id],mult[id],icetvel[id],mrime[id],nmono[id],ice_category))
			}
		}
    	}

    }

    sink()

}
