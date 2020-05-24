Readme file for SCALE-SDM 0.2.5-2.2.2 to reproduce the results of Shima et al. (2019, GMDD).

Corresponding Author: Shin-ichiro Shima (s_shima[at]sim.u-hyogo.ac.jp)

[//]:#(#################################################################################)
# General description
SCALE is a library of weather and climate models of the Earth and planets (Nishizawa et al., 2015; Sato et al.,
2015, http://r-ccs-climate.riken.jp/scale/). We implemented SDM 2.2.0 into SCALE version 0.2.5, 
and constructed a mixed-phase cloud model, SCALE-SDM 0.2.5-2.2.0. 

Then, by correcting three of the four issues (see Shima et al., 2019, GMDD), we updated the model to SCALE-SDM 0.2.5-2.2.1.

Then, by fixing the implementation of the terminal velocity formula of Bohm (1992), and by imposing a limiter Gamma_star<=1 for the deposition of too long ice particles (see Shima et al., 2019, GMDD), we updated the model to SCALE-SDM 0.2.5-2.2.2.

Detailed instruction of SCALE (without SDM) are avaiable from the SCALE web page: http://r-ccs-climate.riken.jp/scale/doc/index.html.

[//]:#(#################################################################################)
# Required software and supported environment
Fortran and C compiler are required to compile SCALE-SDM. 
MPI, NetCDF4, and HDF5 libraries are are also required.

The numerical experiments were conducted by using Intel Fortran/C compiler 17.0.0, SGI MPT 2.15, HDF5 1.8.12, and NetCDF 4.4.1.
For data analysis, R 3.2.3 was used. 

[//]:#(#################################################################################)
# Set environment variable
```
$ export SCALE_SYS=Linux64-intel-impi
```

[//]:#(#################################################################################)
# Set compiler options
In the top directly, you can find a directory named sysdep.
Edit sysdep/Makedef.${SCALE_SYS}, and modify the compiler options according to your system.

For example, 
```
$ cd sysdep
$ vi Makedef.Linux64-intel-impi
----------------
...
-               -mcmodel=medium -heap-arrays
+               -mcmodel=medium -heap-arrays -lmpi

-               -DDEBUG
+               -DDEBUG -lmpi

-FC     = mpiifort
+FC     = ifort

-CC     = mpiicc
-CFLAGS = -O3 -xHost -ip -ftz -mcmodel=medium -shared-intel
+CC     = icc
+CFLAGS = -O3 -xHost -ip -ftz -mcmodel=medium -shared-intel -lmpi

-NETCDF_INCLUDE ?= -I$(NETCDF4)/include
+NETCDF_INCLUDE ?= -I/apps/netcdf/4.4.1/intel/17.0.0/include
-NETCDF_LIBS    ?= -L$(NETCDF4)/lib -L$(HDF5)/lib64 -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lm -lz
+NETCDF_LIBS    ?= -L/apps/netcdf/4.4.1/intel/17.0.0/lib -lnetcdf -lnetcdff -L/apps/hdf5/1.8.12/intel/17.0.0/lib/ -lhdf5_hl -lhdf5 -lm -lz
...
----------------
```

[//]:#(#################################################################################)
# Clean
```
$ cd scale-les/test/case/warmbubble/2D_Khain04_sdm
$ make allclean ENABLE_SDM=T
```

[//]:#(#################################################################################)
# Compile
```
$ cd scale-les/test/case/warmbubble/2D_Khain04_sdm
$ time make -j5 ENABLE_SDM=T
```

[//]:#(#################################################################################)
# Run a batch job
```
$ cd scale-les/test/case/warmbubble/2D_Khain04_sdm
```

Modify the job script according to your system. 
If you are using PBS as the job scheduler, it will be something like this.
```
$ vi iwaya_run.sh
----------------
#!/bin/bash
#PBS -q S
#PBS -l select=4:ncpus=20:mpiprocs=20
#PBS -l walltime=4:00:00
#PBS -N scale-sdm

source /etc/profile.d/modules.sh
cd ${PBS_O_WORKDIR}
module load intel/17.0.0 mpt hdf5/1.8.12 netcdf/4.4.1

# run
mpiexec_mpt dplace -s1 ./scale-les_init init.conf || exit
mpiexec_mpt dplace -s1 ./scale-les  run.conf || exit
----------------
$ qsub  iwaya_run.sh
$ qstat -a
$ tail -f LOG.pe000000
```

[//]:#(#################################################################################)
# Skew-T log-P Diagram of the Initial Atmospheric Sounding
```
$ cd scale-les/test/case/warmbubble/2D_Khain04_sdm
$ cd analyse-env_R/
$ Rscript analyse-env.R
$ convert -rotate 90 -alpha off -density 400 skewTlogP.eps skewTlogP.png
$ ls
alldata.txt  analyse-env.R  skewTlogP.eps  skewTlogP.png  variable_info.txt
$ less variable_info.txt
$ less alldata.txt
$ evince skewTlogP.eps
```

[//]:#(#################################################################################)
# Plot the Spatial Structure of the Cloud Using R
```
$ cd scale-les/test/case/warmbubble/2D_Khain04_sdm
$ cd QHYD_2Dplot
$ Rscript QHYD_2Dplot.R
$ for i in *.pdf ; do convert -density 400 $i ${i/pdf/png} ; done
$ for i in QHYD_overlay ; do convert -delay 20 -loop 0 ${i}.*.png ${i}.gif ; done
```

[//]:#(#################################################################################)
# Plot the Time Series of Water Path and Accumulated Precipitation
```
$ cd scale-les/test/case/warmbubble/2D_Khain04_sdm
$ cd time_series_0D
$ Rscript time_series_0D.R
$ for i in *.pdf; do convert -density 400 $i ${i/pdf/png}; done
```

[//]:#(#################################################################################)
# Make a 2D Color Map of Ice Particle Distribution
```
$ cd scale-les/test/case/warmbubble/2D_Khain04_sdm
$ cd ptl_dist_2Dmap_R
$ Rscript ptl_dist_2Dmap.netcdf.overlay.R
$ for i in *.pdf; do convert -density 400 $i ${i/pdf/png}; done
$ for i in aspect density massratio term_vel_ice ; do convert  -delay 20 -loop 0 ${i}.*.png ${i}.gif ; done
```

[//]:#(#################################################################################)
# Aerosol Size Distribution (Number Density), Freezing Temperature Distribution, Droplet Size Distribution (Mass Density), and Terminal Velocities of Droplets.
```
$ cd ptl_dist_1D_R
$ module load intel/17.0.0 hdf5/1.8.12 netcdf/4.4.1 fftw/3.3.5 gsl/1.16 R
$ Rscript ptl_dist_1D.netcdf.R
$ less numbers.txt
$ for i in *.pdf; do convert -density 400 $i ${i/pdf/png}; done
$ for i in mass_dens_drop num_dens_amsul prob_dens_freezing_temp term_vel_drop xz_SD; do convert -delay 20 -loop 0 ${i}.*.png ${i}.gif; done
```
