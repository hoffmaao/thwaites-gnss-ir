c= 299792458 # m/sec
#   GPS frequencies and wavelengths
fL1 = 1575.42 # MegaHz 154*10.23
fL2 = 1227.60 # 120*10.23
fL5 = 115*10.23
#  GPS wavelengths
wL1 = c/(fL1*1e6) # meters wavelength
wL2 = c/(fL2*1e6)
wL5 = c/(fL5*1e6)
#   galileo frequency values
gal_L1 = 1575.420
gal_L5 = 1176.450
gal_L6 = 1278.70
gal_L7 = 1207.140
gal_L8 = 1191.795
#  galileo wavelengths, meters
wgL1 = c/(gal_L1*1e6)
wgL5 = c/(gal_L5*1e6)
wgL6 = c/(gal_L6*1e6)
wgL7 = c/(gal_L7*1e6)
wgL8 = c/(gal_L8*1e6)

#   beidou frequencies and wavelengths
bei_L2 = 1561.098
bei_L7 = 1207.14
bei_L6 = 1268.52
wbL2 = c/(bei_L2*1e6)
wbL7 = c/(bei_L7*1e6)
wbL6 = c/(bei_L6*1e6)

#   Earth rotation rate used in Nav message
omegaEarth = 7.2921151467E-5 #	%rad/sec
mu = 3.986005e14 # Earth GM value