#!/usr/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 5.0 patchlevel 3    last modified 2016-02-21 
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2016
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')
# set terminal x11 
# set output
set size square
set logscale x
set logscale y
set xrange [ 1.00000e-7 : 1.00000e0 ] 
set yrange [ 1.0e0 : 1.1e3 ]
set key right bottom
set format y "10e%L"
set format x "10e%L"
set term png
set out "density.png"
plot   "1800_ice.dat"every 1 using 3:4 w d, \
       "1200_ice.dat"every 1 using 3:4 w d, \
       "0600_ice.dat"every 1 using 3:4 w d, \
       "0000_ice.dat"every 1 using 3:4 w d
quit
#    EOF
