##
filename <- "../env.txt"
#filename <- "../../../squallline/tutorial_SDM/env.txt"

## libraries
library(RadioSonde)

## constants
gas_d <- 287.06   # specific gas constant for dry air [J/kg/K	z] 
grav <- 9.81      # gravitational acceleration [m/s^2]
press_0 <- 1000.0 # reference pressure [hPa]
cp <- 1.004e3     # specific heat capacity of dry air at constant pressure [J/kg/K]
T0 <- 273.15      # 0 deg C in K [K]

## functions
#### calculate saturation vapor pressure from temperature
#### Teten's formula 
tempc.to.satvp <- function (tempc)
{ # tempc : temperature [deg C]
  para_a1 <- 6.11 #[hPa]
  para_a2 <- 17.27
  para_a3 <- 237.3 #[[deg C]]

  # saturation vapor pressure [hPa]
  satvp <- para_a1 * exp((para_a2 * tempc)/(tempc + para_a3))

  return(satvp)
}

#### calculate dew point temperature from vapor pressure
#### inverse of Teten's formula 
vp.to.dewtempc <- function (vp)
{ # vp : vapor pressure [hPa]
  para_a1 <- 6.11 #[hPa]
  para_a2 <- 17.27
  para_a3 <- 237.3 #[[deg C]]

  tmp_A <- log(vp/para_a1)

  # dew point temperature [deg C]
  dewtempc <- (tmp_A * para_a3)/(para_a2-tmp_A)

  return(dewtempc)
}

#### vapor pressure from mixing ratio
qv_press.to.vp <- function (qv, press)
{ # qv : mixing ratio [g/kg]
  # press : pressure [hPa]
  para_eps <- 0.622 # Rd/Rv []

  # vp : vapor pressure [hPa]
  vp <- (qv/1000.0)*press/(para_eps+(qv/1000.0))

  return(vp)
}

## read and convert data
tmp <- read.table(filename, header=FALSE, fill=TRUE)
sonde <- data.frame(press=tmp$V1,height=tmp$V1,temp=tmp$V2,pot=tmp$V2,qv=tmp$V3,dewpt=tmp$V2)
num_row <- length(sonde$press)

#### height [m]
sonde$height[1] <- 0

#### pressure [hPa] and temperature [deg C]
sonde$temp[1] <- sonde$pot[1]*(sonde$press[1]/press_0)**(gas_d/cp)
for (i in 2:num_row) {
  sonde$press[i] <- sonde$press[i-1] * exp(-grav/gas_d/sonde$temp[i-1]*(sonde$height[i]-sonde$height[i-1]))
  sonde$temp[i]  <- sonde$pot[i] * (sonde$press[i]/press_0)**(gas_d/cp)
}
sonde$temp[1:num_row] <- sonde$temp[1:num_row] - T0

#### vapor pressure [hPa]
tmp_vp <- qv_press.to.vp(sonde$qv[1:num_row],sonde$press[1:num_row])

#### dew point temperature [deg C]
tmp_dewpt <- vp.to.dewtempc(tmp_vp[1:num_row])
sonde <- transform(sonde,dewpt=tmp_dewpt[1:num_row])

#### saturation vapor pressure [hPa]
tmp_satvp <- tempc.to.satvp(sonde$temp[1:num_row])

#### relative humidity [%]
tmp_rh <- 100.0*tmp_vp[1:num_row]/tmp_satvp[1:num_row]
sonde <- transform(sonde,rh=tmp_rh[1:num_row])

### record units
attr(sonde,"units")[match("height",names(sonde))] <- "m"
attr(sonde,"units")[match("press",names(sonde))] <- "hPa"
attr(sonde,"units")[match("temp",names(sonde))] <- "deg C"
attr(sonde,"units")[match("dewpt",names(sonde))] <- "deg C"
attr(sonde,"units")[match("rh",names(sonde))] <- "%"
attr(sonde,"units")[match("pot",names(sonde))] <- "K"
attr(sonde,"units")[match("qv",names(sonde))] <- "g/kg"

## output
postscript("skewTlogP.eps")
title <- paste("skew-T log-P diagram of \n",filename,sep="")
plotsonde(sonde,winds=F,title=title)
skewt.points(sonde$temp,sonde$press)
dev.off()

sink("variable_info.txt")
names(sonde)
attr(sonde,"units")

sink("alldata.txt")
sonde
