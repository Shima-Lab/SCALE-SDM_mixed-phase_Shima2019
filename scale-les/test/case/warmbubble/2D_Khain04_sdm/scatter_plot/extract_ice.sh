#! /bin/bash
for i in {0000..1800..600}
do
  rm ${i}_ice.dat
  touch ${i}_ice.dat
  # equatorial radius, polar radius, average radius, density, freezing temperature, multiplicity
  grep ICE ../SD_output_ASCII_0000000${i}.000.pe* | awk '{ if(($8+$9)/2 > 1.0e-6) print $8,$9,($8+$9)/2,$10,$11,$12}' >> ${i}_ice.dat
done
