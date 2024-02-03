# PCN
pcn_algorithm.R : estimates the pcn from a sample (lattice).
simulations.R : generate samples from a given PCN.
scriptcomaplicacaodados.zip: contains files used in the application  on fire data from Pantanal: biome_border.shp is used to determine the wetland border
MCD64monthly.A2020245.Win06.006.burndate.tif has the matrix with data on fires (it was added to dropbox due to its size).

Script: application.R is used to manipulate the data to obtain the wetland sampling matrix and then run pcn to obtain the spatial dependence of the wetland matrix.
PCN_w.R has the PCN functions adapted to deal with water pixels (value 0).
