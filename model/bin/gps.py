#!usr/bin/env python
# -*- coding: utf-8 -*-
# kristine larson, june 2017
# will try to use this as a library.  however, i am not very competent
import sys
import os
import wget
import numpy as np
import matplotlib.pyplot as plt
from ftplib import FTP
import datetime
from scipy.interpolate import interp1d
import pickle
# don't think it is being used
import math
# do not know what this is
import re
import scipy.signal as spectral

# i think these should be in a class ...
# various numbers you need in the GNSS world
class constants:
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


class wgs84:
    """
    wgs84 parameters
    """
    a = 6378137. # meters
    f  =  1./298.257223563 # flattening factor
    e = np.sqrt(2*f-f**2) # 

def define_filename(station,year,doy,snr):
    """
    given station name, year, doy, snr type
    returns snr filename
    author: Kristine Larson
    """
    xdir = str(os.environ['REFL_CODE'])
    cdoy = '{:03d}'.format(doy)
    cyy = '{:02d}'.format(year-2000)
    f= station + str(cdoy) + '0.' + cyy + '.snr' + str(snr)
    fname = xdir + '/' + str(year) + '/snr/' + station + '/' + f 
    print('snr filename is ', fname) 
    return fname 

def define_filename_prevday(station,year,doy,snr):
    """
    given station name, year, doy, snr type
    returns snr filename for the PREVIOUS day
    author: Kristine Larson
    """
    xdir = str(os.environ['REFL_CODE'])
    year = int(year)
    doy = int(doy)
    if (doy == 1):
        pyear = year -1
        print('found january 1, so previous day is december 31')
        doyx,cdoyx,cyyyy,cyy = ymd2doy(pyear,12,31)
        pdoy = doyx 
    else:
#       doy is decremented by one and year stays the same
        pdoy = doy - 1
        pyear = year

    cdoy = '{:03d}'.format(pdoy)
    cyy = '{:02d}'.format(pyear-2000)
    f= station + str(cdoy) + '0.' + cyy + '.snr' + str(snr)
    fname = xdir + '/' + str(year) + '/snr/' + station + '/' + f 
    print('snr filename for the previous day is ', fname) 
    return fname 

def read_inputs(station):
    """
    given station name, read LSP parameters for strip_snrfile.py
    author: Kristine M Larson
    """
#   directory name is currently defined using REFL_CODE
    xdir = str(os.environ['REFL_CODE'])
    fname = xdir + '/input/' + station
    print('default inputs: ', fname)
#   default location values - not used now
    lat = 0; long = 0; h = 0;
#   counter variables
    k = 0
    lc = 0
#   elevation angle list, degrees
    elang = []
#   azimuth angle list, degrees
    azval = []
    freqs = [] ; reqAmp = []
#   this limits the reflector height values you compute LSP for
    Hlimits = []
#   default polyfit, elevation angle limits for DC removal
    polyFit = 3; pele = [5, 30]
#   desired reflector height precision, in meters
    desiredP = 0.01
    ediff = 2
    noiseRegion = [0.25, 6]
    try:
        f = open(fname, 'r')
#       read all the lines
        l=f.readlines()
#       figure out what is on each
        for line in l:
#           print(line)
            if line[0] == '#':
                lc = lc + 1
            else:
                k=k+1
                nfo=line.split()
                if k==1:
# read in the latitude, longitude, and height. currently this information is not used
                    lat = float(nfo[0])
                    long = float(nfo[1])
                    h = float(nfo[2])
                    print("latitude/longitude/height:", lat,long,h)
# read in the elevation angle ranges
                if k==2:
                    elang.append(float(nfo[0]))
                    elang.append(float(nfo[1]))
                    print('elevation angle limits', elang[0], elang[1])
                    # these are elevation angle limits for polyfit
                    if (len(nfo)==4):
                        pele = [float(nfo[2]), float(nfo[3]) ]
# read in the azimuth ranges
                if k==3:
                    naz = len(nfo)
                    for j in range(naz):
                        azval.append(float(nfo[j]))
# read in the frequencies and the required spectral amplitudes
                if k==4:
                    nnp = int(len(nfo)/2)
                    for j in range(nnp):
                        ni = 2*j; nj = 2*j+1
                        freqs.append(int(nfo[ni]))
                        reqAmp.append(float(nfo[nj]))
                        print('Frequency ', int(nfo[ni]), ' Amplitude ', float(nfo[nj]))
# read in reflector height restrictions etc
                if k==5:
                    np = len(nfo)
                    polyFit = int(nfo[0])
# this is desired precision of the LSP
                    desiredP = float(nfo[1])
                    H1 = Hlimits.append(float(nfo[2]))
                    H2 = Hlimits.append(float(nfo[3]))
                    ediff = float(nfo[4])
                    print('Polynomial Fit ', polyFit, ' Precision: ', desiredP,' ediff ', ediff)
                    if (np > 5):
# range of reflector heights used for noise calculation
                        ns1 = float(nfo[5])
                        ns2 = float(nfo[6])
                        noiseRegion = [ns1,ns2]
                        
        f.close
    except:
        print('some kind of problem reading input file')
        sys.exit()

    return lat,long,h,elang, azval, freqs, reqAmp,polyFit, desiredP, Hlimits, ediff, pele,noiseRegion

def satclock(week, epoch, prn, closest_ephem):
    """
    inputs: gps week, second of week, satellite number (PRN)
    and ephemeris. returns clock clorrection in  meters
    note: although second order correction exists, it is not used  
    """
    # what is sent should be the appropriate ephemeris for given
    # satellite and time
    prn, week, Toc, Af0, Af1, Af2, IODE, Crs, delta_n, M0, Cuc,\
    ecc, Cus, sqrta, Toe, Cic, Loa, Cis, incl, Crc, perigee, radot, idot,\
    l2c, week, l2f, sigma, health, Tgd, IODC, Tob, interval = closest_ephem

    correction = (Af0+Af1*(epoch-Toc))*constants.c
    return correction[0]


def ionofree(L1, L2):
    """
    input are L1 and L2 observables (either phase or pseudorange, in meters)
    output is L3 (meters)
    author: kristine larson
    """
    f1 = constants.fL1
    f2 = constants.fL2
    
    P3 = f1**2/(f1**2-f2**2)*L1-f2**2/(f1**2-f2**2)*L2
    return P3

def azimuth_angle(RecSat, East, North):
    """
    kristine larson
    inputs are receiver satellite vector (meters)
    east and north unit vectors, computed with the up vector
    returns azimuth angle in degrees
    """
    staSatE = East[0]*RecSat[0] + East[1]*RecSat[1] + East[2]*RecSat[2]
    staSatN = North[0]*RecSat[0] + North[1]*RecSat[1] + North[2]*RecSat[2]
#    azangle = 0
    azangle = np.arctan2(staSatE, staSatN)*180/np.pi
    if azangle < 0:
        azangle = 360 + azangle
# 
    return azangle

def rot3(vector, angle):
    """
    input a vector (3) and output the same vector rotated by angle
    in radians apparently.
    from ryan hardy
    """
    rotmat = np.matrix([[ np.cos(angle), np.sin(angle), 0],
                        [-np.sin(angle), np.cos(angle), 0],
                        [             0,             0, 1]])
    vector2 = np.array((rotmat*np.matrix(vector).T).T)[0]
    return vector2

def xyz2llh(xyz, tol):
    """
    inputs are station coordinate vector xyz (x,y,z in meters), tolerance for convergence
    outputs are lat, lon in radians and wgs84 ellipsoidal height in meters
    kristine larson
    """
    x=xyz[0]
    y=xyz[1]
    z=xyz[2]
    lon = np.arctan2(y, x)
    p = np.sqrt(x**2+y**2)
    lat0 = np.arctan((z/p)/(1-wgs84.e**2))
    b = wgs84.a*(1-wgs84.f)
    error = 1
    a2=wgs84.a**2
    i=0 # make sure it doesn't go forever
    while error > tol and i < 10:
        n = a2/np.sqrt(a2*np.cos(lat0)**2+b**2*np.sin(lat0)**2)
        h = p/np.cos(lat0)-n
        lat = np.arctan((z/p)/(1-wgs84.e**2*n/(n+h)))
        error = np.abs(lat-lat0)
        lat0 = lat
        i+=1
    return lat, lon, h

def xyz2llhd(xyz):
    """
    inputs are station vector xyz (x,y,z in meters), tolerance for convergence is hardwired
    outputs are lat, lon in degrees and wgs84 ellipsoidal height in meters
    kristine larson
    """
    x=xyz[0]
    y=xyz[1]
    z=xyz[2]
    lon = np.arctan2(y, x)
    p = np.sqrt(x**2+y**2)
    lat0 = np.arctan((z/p)/(1-wgs84.e**2))
    b = wgs84.a*(1-wgs84.f)
    error = 1
    a2=wgs84.a**2
    i=0 # make sure it doesn't go forever
    tol = 1e-10
    while error > tol and i < 10:
        n = a2/np.sqrt(a2*np.cos(lat0)**2+b**2*np.sin(lat0)**2)
        h = p/np.cos(lat0)-n
        lat = np.arctan((z/p)/(1-wgs84.e**2*n/(n+h)))
        error = np.abs(lat-lat0)
        lat0 = lat
        i+=1
    return lat*180/np.pi, lon*180/np.pi, h



def zenithdelay(h):
    """
    author: kristine larson
    input the station ellipsoidal (height) in meters
    the output is a very simple zenith troposphere delay in meters
    """

    zd = 0.1 + 2.31*np.exp(-h/7000.0)
    return zd

def up(lat,lon):
    """
    author: kristine larson
    inputs rae latitude and longitude in radians
    returns the up unit vector, and local east and north used for azimuth calc.
    """
    xo = np.cos(lat)*np.cos(lon)
    yo = np.cos(lat)*np.sin(lon)
    zo = np.sin(lat)
    u= np.array([xo,yo,zo])    
#    c ... also define local east/north for station: took these from fortran
    North = np.zeros(3)
    East = np.zeros(3)
    North[0] = -np.sin(lat)*np.cos(lon)
    North[1] = -np.sin(lat)*np.sin(lon)
    North[2] = np.cos(lat)
    East[0] = -np.sin(lon)
    East[1] = np.cos(lon)
    East[2] = 0
    return u, East, North

def norm(vect):
    """
    given a three vector - return its norm
    """  
    nv = np.sqrt(np.dot(vect,vect))
    return nv

def elev_angle(up, RecSat):
    """
    inputs:
    up - unit vector in up direction
    RecSat is the Cartesian vector that points from receiver 
    to the satellite in meters
    the output is elevation angle in radians
    author: kristine larson
    """
    ang = np.arccos(np.dot(RecSat,up) / (norm(RecSat)))
    angle = np.pi/2.0 - ang
    return angle

def sp3_interpolator(t, tow, x0, y0, z0, clock0):
    """
    got this from ryan hardy - who coded it for my class?
    inputs are??? tow is GPS seconds
    xyz are the precise satellite coordinates (in meters)
    clocks are likely satellite clock corrections (microseconds)
    i believe n is the order fit, based on what i recall.
    these values do not agree with my test cases in matlab or fortran
    presumably there is an issue with the estimation of the coefficients.
    they are good enough for calculating an elevation angle used in reflectometry
    """
    # ryan set it to 7 - which is not recommended by the paper
    # ryan confirmed that he doesn't know why this doesn't work ...
    n = 7 # don't know why this was being sent before
    coeffs = np.zeros((len(t), 3, n))
    # whatever ryan was doing was not allowed here.  had to make
    # sure these are treated as integers
    s1 = int(-(n-1)/2)
    s2 = int((n-1)/2+1)
#    print(s1,s2)
    
    omega = 2*2*np.pi/(86164.090530833)
    x = np.zeros(len(t))
    y = np.zeros(len(t))
    z = np.zeros(len(t))
    clockf = interp1d(tow, clock0, bounds_error=False, fill_value=clock0[-1])
    clock = clockf(t)
    # looks like it computes it for a number of t values?
    for i in range(len(t)):
        # sets up a matrix with zeros in it - 7 by 7
        independent = np.matrix(np.zeros((n, n)))
        # no idea what this does ...
        m = np.sort(np.argsort(np.abs(tow-t[i]))[:n])
        tinterp = tow[m]-np.median(tow[m])
        # set x, y, and z to zeros
        xr = np.zeros(n)
        yr = np.zeros(n)
        zr = np.zeros(n)
        # coefficients are before and after the time
        for j in range(s1, s2):
            independent[j] = np.cos(np.abs(j)*omega*tinterp-(j > 0)*np.pi/2)
        for j in range(n):
            xr[j], yr[j], zr[j] = rot3(np.array([x0[m], y0[m], z0[m]]).T[j], omega/2*tinterp[j])
 #           print(j, xr[j], yr[j], zr[j])
			
        independent = independent.T
        eig =  np.linalg.eig(independent)
        iinv  = (eig[1]*1/eig[0]*np.eye(n)*np.linalg.inv(eig[1]))
# set up the coefficients
        coeffs[i, 0] = np.array(iinv*np.matrix(xr).T).T[0]
        coeffs[i, 1] = np.array(iinv*np.matrix(yr).T).T[0]
        coeffs[i, 2] = np.array(iinv*np.matrix(zr).T).T[0]
        
        j = np.arange(s1, s2)
        # time since median of the values?
        tx = (t[i]-np.median(tow[m]))
        r_inertial =  np.sum(coeffs[i][:, j]*np.cos(np.abs(j)*omega*tx-(j > 0)*np.pi/2), -1)

        x[i], y[i], z[i] = rot3(r_inertial, -omega/2*tx)
        # returns xyz, in meters ? and satellite clock in microseconds
        return x*1e3, y*1e3, z*1e3, clock
 
    
def readPreciseClock(filename):     
    """
    author: kristine larson
    filename of precise clocks
    returns prn, time (gps seconds of the week), and clock corrections (in meters)
    """          
    StationNFO=open(filename).readlines()
    c= 299792458 # m/sec
 
    nsat = 32
    k=0
# this is for 5 second clocks
    nepochs = 17280
    prn = np.zeros(nepochs*nsat)
    t = np.zeros(nepochs*nsat)
    clockc = np.zeros(nepochs*nsat)
# reads in the high-rate clock file 
    for line in StationNFO:
        if line[0:4] == 'AS G':
            lines= line[7:60]
            sat = int(line[4:6])
            year = int(lines.split()[0])
            month = int(lines.split()[1])
            day = int(lines.split()[2])
            hour = int(lines.split()[3])
            minutes = int(lines.split()[4])
            second = float(lines.split()[5])
            [gw, gpss] = kgpsweek(year, month, day, hour, minutes, second)
            clock = float(lines.split()[7])
            prn[k] = sat
            t[k]=int(gpss)
            clockc[k]=c*clock
            k += 1
    return prn, t, clockc


def ymd2doy(year,month,day):
    """
    takes in year, month, day and returns day of year (doy)
    and the character version of that day of year
    19mar20: and now it returns character versions of 4 and 2 character years
    """
    today=datetime.datetime(year,month,day)
    doy = (today - datetime.datetime(today.year, 1, 1)).days + 1
    cdoy = '{:03d}'.format(doy)
    cyyyy = '{:04d}'.format(year)
    cyy = '{:02d}'.format(year-2000)
    return doy, cdoy, cyyyy, cyy

def rinex_unavco(station, year, month, day):
    """
    author: kristine larson
    picks up a RINEX file from unavco.  it tries to pick up an o file,
    but if it does not work, it tries the "d" version, which must be
    decompressed.  the location of this executable is defined in the crnxpath
    variable. This is from the main unavco directory - not the highrate directory.

    WARNING: only rinex version 2 in this world
    """
    exedir = os.environ['EXE']
    crnxpath ='/usr/local/bin/CRX2RNX '
    doy,cdoy,cyyyy,cyy = ymd2doy(year,month,day)
    rinexfile,rinexfiled = rinex_name(station, year, month, day)
    unavco= 'ftp://data-out.unavco.org'
    filename1 = rinexfile + '.Z'
    filename2 = rinexfiled + '.Z'
    # URL path for the o file and the d file
    url1 = unavco+ '/pub/rinex/obs/' + cyyyy + '/' + cdoy + '/' + filename1
    url2 = unavco+ '/pub/rinex/obs/' + cyyyy + '/' + cdoy + '/' + filename2
    print(url1)
    print(url2)
    try:
        print('try to get o file')
        wget.download(url1,filename1)
        cmd = 'uncompress ' + filename1; os.system(cmd)
        print('found it ')
    except:
        print('did not find o file')
        try:
            wget.download(url2,filename2)
            cmd = 'uncompress ' + filename2; os.system(cmd)
            #convert
            cmd = crnxpath + rinexfiled; os.system(cmd)
            #remove compressed file
            cmd = 'rm -f ' + rinexfiled; os.system(cmd)
            print('found d file and converted to o file')
        except:
            print('failed to find either RINEX file at unavco')


def rinex_sopac(station, year, month, day):
    """
    author: kristine larson
    inputs: station name, year, month, day
    picks up a hatanaka RINEX file from SOPAC - converts to o
    hatanaka exe hardwired  for my machine
    """
    exedir = os.environ['EXE']
    crnxpath = '/usr/local/bin/CRX2RNX '
    doy,cdoy,cyyyy,cyy = ymd2doy(year,month,day)
    sopac = 'ftp://garner.ucsd.edu'
    oname,fname = rinex_name(station, year, month, day) 
    file1 = fname + '.Z'
    path1 = '/pub/rinex/' + cyyyy + '/' + cdoy + '/' 
    url = sopac + path1 + file1 
    print(url)
    try:
        wget.download(url,file1)
        cmd = 'uncompress ' + file1 ; os.system(cmd)
        cmd = crnxpath + fname; os.system(cmd)
        #remove compressed file
        cmd = 'rm -f ' + fname;  os.system(cmd)
        print('successful download from SOPAC ')
    except:
        print('some kind of problem with download',file1)
        cmd = 'rm -f ' + file1
        os.system(cmd)

def getnavfile(year, month, day):
    """
    author: kristine larson
    given year, month, day it picks up a GPS nav file from SOPAC
    and stores it
    returns the name of the file and its directory

    """
    doy,cdoy,cyyyy,cyy = ymd2doy(year,month,day)
    sopac = 'ftp://garner.ucsd.edu'
    navname,navdir = nav_name(year, month, day)
    file1 = navname + '.Z'
    path1 = '/pub/rinex/' + cyyyy + '/' + cdoy + '/'
    url = sopac + path1 + file1
    if (os.path.isfile(navdir + '/' + navname ) == True):
        print('nav file already exists')
    else:
        print('pick up the nav file ')
        try:
            wget.download(url,file1)
            cmd = 'uncompress ' + file1
            os.system(cmd)
            store_orbitfile(navname,year,'nav') 
        except:
            print('some kind of problem with nav download',navname)
            cmd = 'rm -f ' + file1
            os.system(cmd)

    return navname,navdir

def getsp3file(year,month,day):
    """
    author: kristine larson
    retrieves IGS sp3 precise orbit file from CDDIS
    inputs are year, month, and day 
    modified in 2019 to use wget 
    returns the name of the file and its directory
    """
    name, fdir = sp3_name(year,month,day,'igs') 
    print(name)
    print(fdir)
    cddis = 'ftp://cddis.nasa.gov'
    if (os.path.isfile(fdir + '/' + name ) == True):
        print('sp3file already exists')
    else:
        gps_week = name[3:7]
        file1 = name + '.Z'
        filename1 = '/gnss/products/' + str(gps_week) + '/' + file1
        url = cddis + filename1 
        print(url)
        try:
            wget.download(url,file1)
            cmd = 'uncompress ' + file1
            os.system(cmd)
            store_orbitfile(name,year,'sp3') 
        except:
            print('some kind of problem-remove empty file')
            cmd = 'rm -f ' + file1
            os.system(cmd)

#   return the name of the file so that if you want to store it
    return name, fdir

def getsp3file_flex(year,month,day,pCtr):
    """
    author: kristine larson
    retrieves sp3 orbit files from CDDIS
    inputs are year, month, and day  (integers), and 
    pCtr, the processing center  (3 characters)
    returns the name of the file and its directory
    """
    # returns name and the directory
    name, fdir = sp3_name(year,month,day,pCtr) 
    print(name)
    print(fdir)
    gps_week = name[3:7]
    file1 = pCtr + name[3:8] + '.sp3.Z'
    name = pCtr + name[3:8] + '.sp3'
    if (os.path.isfile(fdir + '/' + name ) == True):
        print('sp3file already exists')
    else:
        filename1 = '/gnss/products/' + str(gps_week) + '/' + file1
        cddis = 'ftp://cddis.nasa.gov'
        url = cddis + filename1 
        print(url)
        try:
            wget.download(url,file1)
            cmd = 'uncompress ' + file1
            os.system(cmd)
            store_orbitfile(name,year,'sp3') 
        except:
            print('some kind of problem-remove empty file')
            cmd = 'rm -f ' + file1
            os.system(cmd)
#   return the name of the file so that if you want to store it
    return name, fdir

def getsp3file_mgex(year,month,day,pCtr):
    """
    author: kristine larson
    retrieves MGEX sp3 orbit files 
    inputs are year, month, and day  (integers), and 
    pCtr, the processing center  (3 characters)
    right now it checks for the "new" name, but in reality, it 
    assumes you are going to use the GFZ product
    """
    # this returns sp3 orbit product name
    name, fdir = sp3_name(year,month,day,pCtr) 
    gps_week = name[3:7]
    file1 = name + '.Z'
    print(file1)

    # get the sp3 filename for the new format
    doy,cdoy,cyyyy,cyy = ymd2doy(year,month,day)
    file2 = 'GFZ0MGXRAP_' + cyyyy + cdoy + '0000_01D_05M_ORB.SP3.gz'
    print(file2)
    name2 = file2[:-3] 

    # where the files live at CDDIS
    cddis = 'ftp://cddis.nasa.gov'
    dirlocation = '/gps/products/mgex/' + str(gps_week) + '/'
    url = cddis + dirlocation  + file1; print(url)
    url2 = cddis + dirlocation + file2; print(url2)
    mgex = 0
    if (os.path.isfile(fdir + '/' + name ) == True):
        print('first kind of MGEX sp3file already exists')
        mgex = 1
    if (os.path.isfile(fdir + '/' + name2 ) == True):
        print('second kind of MGEX sp3file already exists')
        mgex = 2
# there has to be a better way ... but for now  this works
# only try to download if neither exists
    if (mgex == 2):
        name = name2
    if (mgex == 1):
        name = file1[:-2]
    if (mgex == 0):
        try:
            wget.download(url,file1)
            cmd = 'uncompress ' + file1
            os.system(cmd)
            name = file1[:-2]
        # store the file in its proper place
            store_orbitfile(name,year,'sp3') 
        except:
            print('some kind of problem trying to get first file')
            cmd = 'rm -f ' + file1
            os.system(cmd)
            name = file2[:-3]
        # try the second file
            try:
                wget.download(url2,file2)
                cmd = 'gunzip ' + file2
                os.system(cmd)
                # name to return
                name = file2[:-3]
            # store the file in its proper place
                store_orbitfile(name,year,'sp3') 
            except:
                print('some kind of problem downloading 2nd kind of MGEX file')

    return name, fdir

def codclock(year,month,day):
    """
    author: kristine lasron
    pick up 5 second clocks from the bernese group...
    """
    n,nn=igsname(year,month,day)
    gps_week = n[3:7]     
    print(nn)
    file1 = nn + '.Z'
    filepath1 = file1
    filename1 = '/gnss/products/' + str(gps_week) + '/'+file1
    try:
        ftp = FTP('cddis.gsfc.nasa.gov')
        ftp.login()
    # not sure why they do this
        f1 = open(filepath1,'wb')
        print('Retrieving: '+file1)
        ftp.retrbinary("RETR " + filename1,f1.write)
        f1.close()
        cmd = 'gunzip -f ' + filepath1
        print(cmd)
        os.system(cmd)
    except:
        print('some kind of problem-remove empty file')
        cmd = 'rm -f ' + filepath1
        os.system(cmd)
    

def kgpsweek(year, month, day, hour, minute, second):
    """
    inputs are year (4 char), month, day, hour, minute, second
    outputs: gps week and second of the week
    author: kristine larson
    """

    year = np.int(year)
    M = np.int(month)
    D = np.int(day)
    H = np.int(hour)
    minute = np.int(minute)
    
    UT=H+minute/60.0 + second/3600. 
    if M > 2:
        y=year
        m=M
    else:
        y=year-1
        m=M+12
        
    JD=np.floor(365.25*y) + np.floor(30.6001*(m+1)) + D + (UT/24.0) + 1720981.5
    GPS_wk=np.floor((JD-2444244.5)/7.0);
    GPS_wk = np.int(GPS_wk)
    GPS_sec_wk=np.rint( ( ((JD-2444244.5)/7)-GPS_wk)*7*24*3600)            
     
    return GPS_wk, GPS_sec_wk
def kgpsweekC(z):
    """
    takes in time tag from a RINEX file and converts to gps week/sec
    so the input is a character string of length 26.  in this 
    kind of string, the year is only two characters
    author: kristine larson
    """
    y= np.int(z[1:3])
    m = np.int(z[4:6])
    d=np.int(z[7:9])
    hr=np.int(z[10:12])
    mi=np.int(z[13:15])
    sec=np.float(z[16:26])
    gpsw,gpss = kgpsweek(y+2000,m,d,hr,mi,sec)
    return gpsw, gpss

def igsname(year,month,day):
    """
    take in year, month, day
    returns IGS sp3 filename and COD clockname (5 sec)
    author: kristine larson
    """
    [wk,sec]=kgpsweek(year,month,day,0,0,0)
    x=int(sec/86400)
    dd = str(wk) + str(x) 
    name = 'igs' + str(wk) + str(x) + '.sp3'
    # i think at some point htey changed to lower case?
   # clockname = 'COD' + dd + '.CLK_05S.txt'
    clockname = 'cod' + dd + '.clk_05s'

    return name, clockname
def read_sp3(file):
    """
    borrowed from Ryan Hardy, who got it from David Wiese
    """
    try:      
        f = open(file)
        raw = f.read()
        f.close()
        lines  = raw.splitlines()
        nprn = np.int(lines[2].split()[1])
        lines  = raw.splitlines()[22:-1]
        epochs = lines[::(nprn+1)]
        nepoch =  len(lines[::(nprn+1)])
        week, tow, x, y, z, clock, prn = np.zeros((nepoch*nprn, 7)).T
        for i in range(nepoch):
            year, month, day, hour, minute, second = np.array(epochs[i].split()[1:], dtype=float)
            week[i*nprn:(i+1)*nprn], tow[i*nprn:(i+1)*nprn] = \
				kgpsweek(year, month, day, hour, minute, second)
            for j in range(nprn):
                prn[i*nprn+j] =  int(lines[i*(nprn+1)+j+1][2:4])
                x[i*nprn+j] = np.float(lines[i*(nprn+1)+j+1][4:18])
                y[i*nprn+j] = np.float(lines[i*(nprn+1)+j+1][18:32])
                z[i*nprn+j] = np.float(lines[i*(nprn+1)+j+1][32:46])
                clock[i*nprn+j] = np.float(lines[(i)*(nprn+1)+j+1][46:60])
    except:
        print('sorry - the sp3file does not exist')
        week,tow,x,y,z,prn,clock=[0,0,0,0,0,0,0]
		
    return week, tow, prn, x, y, z, clock

def myreadnav(file):
    """
    input is navfile name
    output is complicated - broadcast ephemeris blocks
    author: Kristine Larson, April 2017
    """
# input is the nav file
    try:
        f = open(file, 'r')
        nav = f.read()
        f.close()
        nephem = (len(nav.split('END OF HEADER')[1].splitlines())-1)/8
        nephem = int(nephem) #    print(nephem)         
        lines = nav.split('END OF HEADER')[1].splitlines()[1:]
        table = np.zeros((nephem, 32))
        print('Total number of ephemeris messages',nephem)
        for i in range(nephem):
            for j in range(8):
                if j == 0:
                    prn = int(lines[i*8+j][:2])
                    year = int(lines[i*8+j].split()[1])
                    if year > 76:
                        year += 1900
                    else:
                        year += 2000
                    month = int(lines[i*8+j].split()[2])
                    day = int(lines[i*8+j].split()[3])
                    hour = int(lines[i*8+j].split()[4])
                    minute = int(lines[i*8+j].split()[5])
                    second = float(lines[i*8+j][17:22])
                    table[i, 0] = prn
#                    print('Ephem for: ', prn, year, month, day, hour, minute)
                    week, Toc = kgpsweek(year, month, day, hour, minute, second)
                    table[i, 1] =  week
                    table[i, 2] = Toc
                    Af0 = np.float(lines[i*8][-3*19:-2*19].replace('D', 'E'))
                    Af1 = np.float(lines[i*8][-2*19:-1*19].replace('D', 'E'))
                    Af2 = np.float(lines[i*8][-19:].replace('D', 'E'))
                    table[i,3:6] = Af0, Af1, Af2
                elif j != 7:
                    for k in range(4):
                        value = np.float(lines[i*8+j][19*k+3:19*(k+1)+3].replace('D', 'E'))
                        table[i,2+4*j+k] = value
                elif j== 7:
                    table[i,-2]= np.float(lines[i*8+j][3:19+3].replace('D', 'E'))
                    if not lines[i*8+7][22:].replace('D', 'E').isalpha():
                        table[i,-1]= 0
                    else:
                        table[i, -1] = np.float(lines[i*8+7][22:41].replace('D', 'E'))
# output is stored as:
#
# 0-10   prn, week, Toc, Af0, Af1, Af2, IODE, Crs, delta_n, M0, Cuc,\
# 11-22    ecc, Cus, sqrta, Toe, Cic, Loa, Cis, incl, Crc, perigee, radot, idot,\
# 23-24?                   l2c, week, l2f, sigma, health, Tgd, IODC, Tob, interval 
# week would be 24 by this scheme?
# Toe would be 14 
#	
        ephem = table
    except:
        print('This ephemeris file does not exist',file)
        ephem = []
    return ephem
def myfindephem(week, sweek, ephem, prn):
    """
# inputs are gps week, seconds of week
# ephemerides and PRN number
# returns the closest ephemeris block after the epoch
# if one does not exist, returns the first one    
    author: kristine larson
"""
    t = week*86400*7+sweek
# defines the TOE in all the ephemeris 
# he is taking the week and adding ToE (14?)
# poorly coded is all i'm gonna say

    teph = ephem[:, 24]*86400*7+ephem[:, 14]     
    prnmask = np.where(ephem[:, 0]== prn)    

    [nr,nc]=np.shape(prnmask)
#    print(nr,nc)
    if nc == 0:
        print('no ephemeris for that PRN number')
        closest_ephem = []
    else:
        try:          
            signmask = np.where(t >= teph[prnmask])
            proxmask =  np.argmin(t-teph[prnmask][signmask])
            closest_ephem = ephem[prnmask][signmask][proxmask]
        except:
#           print('using first ephemeris - but not after epoch')
            closest_ephem = ephem[prnmask][0]
        
  
    return closest_ephem

def findConstell(cc):
    """
    input is one character (from rinex satellite line)
    output is integer added to the satellite number
    0 for GPS, 100 for Glonass, 200 for Galileo, 300 for everything else?
    author: kristine larson, GFZ, April 2017
    """
    if (cc == 'G' or cc == ' '):
        out = 0
    elif (cc == 'R'): # glonass
        out = 100
    elif (cc == 'E'): # galileo
        out = 200
    else:
        out = 300
        
    return out
def myscan(rinexfile):
    """
    stripping the header code came from pyrinex.  
    data are stored into a variable called table
    columns 0,1,2 are PRN, GPS week, GPS seconds, and observables
    rows are the different observations. these should be stored 
    properly - this is a kluge
    """
    f=open(rinexfile,'r')
    lines = f.read().splitlines(True)
    lines.append('')
    # setting up a set or directionary
    # sets must be unique - so that is hwy he checks to see if it already exists
    header={}        
    eoh=0
# looks like it reads all the header lines, so you can extract them as you like
    for i,line in enumerate(lines):
        if "END OF HEADER" in line:
            eoh=i
            break
#        print(line[60:].strip())
        if line[60:].strip() not in header:
            header[line[60:].strip()] = line[:60].strip()
        else:
            header[line[60:].strip()] += " "+line[:60].strip()
    
    header['APPROX POSITION XYZ'] = [float(i) for i in header['APPROX POSITION XYZ'].split()]
    w = header['APPROX POSITION XYZ']
#    print(w)
#    approxpos = [float(i) for i in header['APPROX POSITION XYZ'].split()]
    header['# / TYPES OF OBSERV'] = header['# / TYPES OF OBSERV'].split()
#    typesObs = header['# / TYPES OF OBSERV'].split()
    aa=header['# / TYPES OF OBSERV']
    types = aa[1:] # this means from element 1 to the end
    # these are from the pyrinex verison of hte code
    header['# / TYPES OF OBSERV'][0] = int(header['# / TYPES OF OBSERV'][0])
    header['INTERVAL'] = float(header['INTERVAL'])
    # need to get approx position of the receiver
    x,y,z = header['APPROX POSITION XYZ']
    numobs = int(np.ceil(header['# / TYPES OF OBSERV'][0]))

# initial three columns in the newheader variable
    newheader = 'PRN\tWEEK\tTOW'
# add in the observation types from this file
    for j in range(numobs):
        newheader += '\t'+types[j]
# header using the Ryan Hardy style
#    print(newheader)

#    # set the line reader to after the end of the header
    # this tells it where to start
    i=eoh+1
    # so try to implement ryan hardy's storing procedure, where 0 is prn,
    # 1 is week, 2 is seconds of week, 3-N are the observables
    table = np.zeros((0, numobs+3))
    print('number of observables ', numobs)
    if numobs > 10:
        print('Tooooo many observables. I cannot deal with this')
        return
    print('line number ' , eoh)
    l = 0 # start counter for blocks
    while True:
        if not lines[i]: break
        if not int(lines[i][28]):
#            print(lines[i])
#            z=lines
            [gw,gs]=kgpsweekC(lines[i])
#            print('week and sec', gw,gs)    
            numsvs = int(lines[i][30:32])  # Number of visible satellites at epoch
#            print('number of satellites',numsvs)
            # strictly speaking i don't understand this line.
            table = np.append(table, np.zeros((numsvs, numobs+3)), axis=0)
            table[l:l+numsvs, 1] = gw
            table[l:l+numsvs, 2] = gs
            #headlength.append(1 + numsvs//12)
            sp = []
                
            if(numsvs>12):
                for s in range(numsvs):
                    xv = findConstell(lines[i][32+(s%12)*3:33+(s%12)*3]) 
                    sp.append(xv + int(lines[i][33+(s%12)*3:35+(s%12)*3]))
                    if s>0 and s%12 == 0:
                        i+= 1  # For every 12th satellite there will be a new row with satellite names                sats.append(sp) # attach satellites here
            else:
                for s in range(numsvs):
                    xv = findConstell(lines[i][32+(s%12)*3:33+(s%12)*3])
                    sp.append(xv + int(lines[i][33+(s%12)*3:35+(s%12)*3]))

#            print(len(sp), 'satellites in this block', sp)
            for k in range(numsvs): 
                table[l+k,0]=sp[k]
                if (numobs > 5):
                    for d in range(5):
                        gg = d*16
                        f=lines[i+1+2*k][gg:gg+14]
                        if not(f == '' or f.isspace()):
                            val = np.float(lines[i+1+2*k][gg:gg+14])
                            table[l+k, 3+d] = val
                    for d in range(numobs-5):
                        gg = d*16
                        f=lines[i+2+2*k][gg:gg+14]
                        if not (f == '' or f.isspace()):
                            val = np.float(lines[i+2+2*k][gg:gg+14])
                            table[l+k, 3+5+d] = val
                else:
                    for d in range(numobs):
                        gg = d*16
                        f = lines[i+1+2*k][gg:gg+14]
                        if (f == '' or f.isspace()):
                            val = np.float(lines[i+2+2*k][gg:gg+14])
                            table[l+k, 3+d] = val

            i+=numsvs*int(np.ceil(header['# / TYPES OF OBSERV'][0]/5))+1
            l+=numsvs
        else:
            print('there was a comment or some header info in the rinex')
            flag=int(lines[i][28])
            if(flag!=4):
                print('this is a flag that is not 4', flag)
            skip=int(lines[i][30:32])
            print('skip this many lines',skip)
            i+=skip+1
            print('You are now on line number',i)
    [nr,nc]=np.shape(table)
    print('size of the table variable is ',nr, ' by ', nc)    
    # code provided by ryan hardy - but not needed???
#    format = tuple(np.concatenate((('%02i', '%4i', '%11.7f'), 
#							((0, numobs+3)[-1]-3)*['%14.3f'])))
# I think this takes the newheader and combines it with the information in
# the table variable
    obs =  dict(zip(tuple(newheader.split('\t')), table.T))
    return obs,x,y,z

def geometric_rangePO(secweek, prn, rrec0, sweek, ssec, sprn, sx, sy, sz, sclock):
    """
    Calculates and returns geometric range (in metres) given
    time (week and sec of week), prn, Cartesisan 
    receiver coordinates rrec0(meters)
    using the precise ephemeris instead of the broadcast
    returns clock correction (clockC) in meters
    Author: Kristine Larson, May 2017
    June 21, 2017 returns transmit time (in seconds) so I can calculate
    relatistic correction
    """
    error = 1
    # find the correct sp3 data for prn
    m = [sprn==prn]
    nx,ny,nz,nc = sp3_interpolator(secweek, ssec[m], sx[m], sy[m], sz[m], sclock[m]) 
    SatOrb = np.array([nx[0],ny[0],nz[0]])
    geo=norm(SatOrb-rrec0)
    c=constants.c
    clockC = nc*1e-6*c
    oE = constants.omegaEarth
    deltaT = norm(SatOrb - rrec0)/constants.c
    ij=0
    while error > 1e-16:  
        nx,ny,nz,nc = sp3_interpolator(secweek-deltaT, ssec[m], sx[m], sy[m], sz[m], sclock[m]) 
        SatOrb = np.array([nx[0],ny[0],nz[0]])
#        SatOrb, relcorr = propagate(week, secweek-deltaT, closest_ephem)
        Th = -oE * deltaT
        xs = SatOrb[0]*np.cos(Th)-SatOrb[1]*np.sin(Th)
        ys = SatOrb[0]*np.sin(Th)+SatOrb[1]*np.cos(Th)
        SatOrbn = [xs, ys, SatOrb[2]]
        geo=norm(SatOrbn-rrec0)
        deltaT_new = norm(SatOrbn-rrec0)/constants.c               
        error = np.abs(deltaT - deltaT_new)
        deltaT = deltaT_new
        ij+=1
    #    print(ij)
    return geo,SatOrbn,clockC, deltaT
def read_files(year,month,day,station):
    """
    """   

    doy,cdoy,cyyyy,cyy = ymd2doy(year,month,day)
    # i have a function for this ....
    rinexfile = station + cdoy + '0.' + cyy + 'o'
    navfilename = 'auto'  + cdoy + '0.' + cyy +  'n'
    if os.path.isfile(rinexfile):
        print('rinexfile exists')
    else:
        print(rinexfile)
        print('get the rinex file')
        rinex_unavco(station, year, month, day)
    # organize the file names
    print('get the sp3 and clock file names')
    sp3file, cname = igsname(year,month,day)

    # define some names of files
    if os.path.isfile(navfilename):
        print('nav exists')
    else:
        print('get nav')
        navname,navdir = getnavfile(year,month,day)
    print('read in the broadcast ephemeris')
    ephemdata = myreadnav(navfilename)
    if os.path.isfile(cname):
        print('file exists')
    else:        
        print('get the CODE clock file')
        codclock(year,month,day)
    pname = cname[0:9] + 'pckl'
    print('pickle', pname)
    # if file exists already     
    if os.path.isfile(pname):
        print('read existing pickle file')
        f = open(pname, 'rb')
        [prns,ts,clks] = pickle.load(f)
        f.close()
    else:
        print('read and save as pickle')
        prns, ts, clks = readPreciseClock(cname)
        # and then save them
        f = open(pname, 'wb')
        pickle.dump([prns,ts,clks], f)
        f.close()
    if os.path.isfile(sp3file):
        print('sp3 exsts')
    else:
        print('get sp3')
        getsp3file(year,month,day)
    print('read in the sp3 file', sp3file)
    
    sweek, ssec, sprn, sx, sy, sz, sclock = read_sp3(sp3file)
    
#20    print('len returned data', len(ephemdata), navfilename
    rinexpickle = rinexfile[0:11] + 'pclk'
    if os.path.isfile(rinexpickle):
        print('rinex pickle exists')
        f=open(rinexpickle,'rb')
        [obs,x,y,z]=pickle.load(f)
        f.close()
    else:     
        print('read the RINEX file ', rinexfile)
        obs,x,y,z = myscan(rinexfile)
        print('save as pickle file')
        f=open(rinexpickle,'wb')
        pickle.dump([obs,x,y,z], f)
        f.close()
        
    return ephemdata, prns, ts, clks, sweek, ssec, sprn, sx, sy,sz,sclock,obs,x,y,z
def precise_clock(prn,prns,ts,clks,gps_seconds):
    """
    input prn and gps_seconds, and contents of precise clock file,
    i think.
    units of returned variable are?
    """
    m = [prns == prn]
    preciset = ts[m]
    precisec = clks[m]
#           # trying to interpolate
# could do this much faster ...
    clockf = interp1d(preciset, precisec, bounds_error=False)
    scprecise = clockf(gps_seconds)
    return scprecise
def precise_clock_test(ttags, prns,ts,clks):
    """
    input timetags (gps seconds of the week)
    and precise clock values.returns interpolated values
    units of returned variable are?
    """
    NP = len(ttags)
    print('number of time tags', NP)
    newclks = np.zeros((NP, 32))
    for i in range(32): 
        prn = i+1
 #        print('checking', prn)
        m = [prns == prn]
        [nr,nc]=np.shape(m)
 #       print(nr,nc)
        preciset = ts[m]
        precisec = clks[m]
        clockf = interp1d(preciset, precisec, bounds_error=False)
        k=0
        for t in ttags:
            newclks[k,i]= clockf(t)
            k +=1
#        print(prn, numrows, numcols)
    return newclks

def propagate(week, sec_of_week, ephem):
    """
    inputs are GPS week, seconds of the week, and the appropriate 
    ephemeris block from the navigation message
    returns the x,y,z, coordinates of the satellite 
    and relativity correction (also in meters), so you add,
    not subtract
    Kristine Larson, April 2017

    """

# redefine the ephem variable
    prn, week, Toc, Af0, Af1, Af2, IODE, Crs, delta_n, M0, Cuc,\
    ecc, Cus, sqrta, Toe, Cic, Loa, Cis, incl, Crc, perigee, radot, idot,\
    l2c, week, l2f, sigma, health, Tgd, IODC, Tob, interval = ephem
    sweek = sec_of_week
    # semi-major axis
    a = sqrta**2
    t = week*7*86400+sweek
    tk = t-Toe
    # no idea if Ryan Hardy is doing this correctly - it should be in a function
    tk  =  (tk - 302400) % (302400*2) - 302400
    n0 = np.sqrt(constants.mu/a**3)
    n = n0+ delta_n
    Mk = M0 + n*tk
    i = 0
    Ek = Mk
    E0 = Mk + ecc*np.sin(Mk)
    # solve kepler's equation
    while(i < 15 or np.abs(Ek-E0) > 1e-12):
        i +=1
        Ek = Mk + ecc*np.sin(E0) 
        E0 = Mk + ecc*np.sin(Ek)
    nuk = np.arctan2(np.sqrt(1-ecc**2)*np.sin(Ek),np.cos(Ek)-ecc)
    Phik = nuk + perigee
    duk = Cus*np.sin(2*Phik)+Cuc*np.cos(2*Phik)
    drk = Crs*np.sin(2*Phik)+Crc*np.cos(2*Phik)
    dik = Cis*np.sin(2*Phik)+Cic*np.cos(2*Phik)
    uk = Phik + duk
    rk = a*(1-ecc*np.cos(Ek))+drk
       
    ik = incl+dik+idot*tk
    xkp = rk*np.cos(uk)
    ykp = rk*np.sin(uk)
    Omegak = Loa + (radot-constants.omegaEarth)*tk -constants.omegaEarth*Toe
    xk = xkp*np.cos(Omegak)-ykp*np.cos(ik)*np.sin(Omegak)
    yk = xkp*np.sin(Omegak)+ykp*np.cos(ik)*np.cos(Omegak)
    zk = ykp*np.sin(ik) 
    # using class
    F = -2*np.sqrt(constants.mu)/constants.c
    relcorr = F*ecc*sqrta*np.sin(Ek)
#    return [xk, yk, zk], relcorr
    return [xk[0], yk[0], zk[0]], relcorr

def mygeometric_range(week, secweek, prn, rrec0, closest_ephem):
    """
    Calculates and returns geometric range (in metres) given
    time (week and sec of week), prn, receiver coordinates (cartesian, meters)
    this assumes someone was nice enough to send you the closest ephemeris
    returns the satellite coordinates as well, so you can use htem
    in the A matrix
    Kristine Larson, April 2017
    """
    error = 1

    SatOrb, relcorr = propagate(week, secweek, closest_ephem)
    # first estimate of the geometric range
    geo=norm(SatOrb-rrec0)

    deltaT = norm(SatOrb - rrec0)/constants.c
    while error > 1e-16:     
        SatOrb, relcorr = propagate(week, secweek-deltaT, closest_ephem)
        Th = -constants.omegaEarth * deltaT
        xs = SatOrb[0]*np.cos(Th)-SatOrb[1]*np.sin(Th)
        ys = SatOrb[0]*np.sin(Th)+SatOrb[1]*np.cos(Th)
        SatOrbn = [xs, ys, SatOrb[2]]
        geo=norm(SatOrbn-rrec0)
        deltaT_new = norm(SatOrbn-rrec0)/constants.c               
        error = np.abs(deltaT - deltaT_new)
        deltaT = deltaT_new
    return geo,SatOrbn


def readobs(file,nepochX):
    """
    inputs: filename and number of epochs from the RINEX file you want to unpack
    returns: receiver Cartesian coordinates (meters) and observation blocks
    Kristine Larson, April 2017
    18aug20: updated so that 0,0,0 is returned if there are no receiver coordinates
    """
    f = open(file, 'r')
    obs = f.read()
    f.close()
    testV = obs.split('APPROX POSITION XYZ')[0].split('\n')[-1].split()
#   set default receiver location values,              
    x0=0
    y0=0
    z0=0
    if len(testV) == 3:
        x0, y0, z0 = np.array(obs.split('APPROX POSITION XYZ')[0].split('\n')[-1].split(), 
    
							dtype=float)
    types = obs.split('# / TYPES OF OBSERV')[0].split('\n')[-1].split()[1:]
    print(types)
    tfirst =  obs.split('TIME OF FIRST OBS')[0].split('\n')[-1].split()
    print(tfirst)
    ntypes = len(types)
	#Identify unique epochs
    countstr = '\n '+tfirst[0][2:]+' '+'%2i' % np.int(tfirst[1])+' '+'%2i' % np.int(tfirst[2])
    print(countstr) # so this is 07 10 13, e.g.
    # figures out how many times this string appears between END OF HEADER and end of file
    nepoch = obs.split('END OF HEADER')[1].count(countstr)
    
    table = np.zeros((0, ntypes+3))
	#Identify kinds of observations in the file
    header = 'PRN\tWEEK\tTOW'
    t0 = 0
    for i in range(ntypes):
        header += '\t'+types[i]
    l = 0
#    print('countstr',countstr[0])
    print('header', header)
    epochstr_master = re.split(countstr, obs)
    print('number of data epochs in the file', nepoch)
    # only unpack a limited number of epochs    #nepoch = 5
    print('restricted number of epochs to be read ', nepochX)
    for i in range(nepochX):
        epochstr = epochstr_master[i+1]        
        # this will only work with new rinex files that use G as descriptor
# this only finds GPS!
        prnstr = re.findall(r'G\d\d|G \d', epochstr.split('\n')[0])
        print(prnstr)
#        print('prnstr', prnstr)
        # number of satellites for an observation block
        # this code won't work if there are more than 12 satellites (i think)
        nprn = len(prnstr)
        # append
        table = np.append(table, np.zeros((nprn, ntypes+3)), axis=0)
        # decode year, month, day hour, minute ,second
        # convert 2 char to 4 char year
        # but it appears that the year nad month and day 
        # are being hardwired from the header, which is both  very very odd
        # and stupid
        year = int(tfirst[0])
        month = int(tfirst[1])
        day = int(tfirst[2])
        hour = int(epochstr.split()[0])
        minute = int(epochstr.split()[1])
        second = float(epochstr.split()[2])
        print(year, month, day, hour, minute, second)
        week, tow = kgpsweek(year, month, day, hour, minute, second)
        print('reading week number ', week, 'sec of week', tow)
        table[l:l+nprn, 1] = week
        table[l:l+nprn,2] = tow
        for j in range(nprn):
            table[l+j, 0] = np.int(prnstr[j][1:])
            # split up by end of line markers, which makes sense
#            print('something',2*j+1)
            line0 = epochstr.split('\n')[2*j+1]#.split()
            line = re.findall('.{%s}' % 16, line0+' '*(80-len(line0)))
# this seems to get strings of a certain width, which he will later change 
#through float# also where he made his mistake on reading last two characters 
#            print(j,len(line), line)

            for k in range(len(line)):
                if line[k].isspace():
                    table[l+j, 3+k] = np.nan
                    continue
                table[l+j, 3+k] = np.float(line[k][:-2])
#                print(table[l+j,3+k])
                # if more than 5 observaitons, he has to read the next line
            if ntypes > 5:
                line2 = epochstr.split('\n')[2*j+2].split()
#                print(line2)
                for k in range(len(line2)):
                    table[l+j, 3+len(line)+k] = np.float(line2[k][:-2])
        l += nprn	
    format = tuple(np.concatenate((('%02i', '%4i', '%11.7f'), 
							((0, ntypes+3)[-1]-3)*['%14.3f'])))
    # kinda strange - but ok - format statemenst for PRN, week, sec of week, etc
#    print(format)
# I guess it is making a super observable so that it can later be parsed.
# very strange 
    obs =  dict(zip(tuple(header.split('\t')), table.T))

    return np.array([x0,y0,z0]),obs

def tmpsoln(navfilename,obsfilename):
    """
    kristine larson
    inputs are navfile and obsfile names
    should compute a pseudorange solution
    """
    r2d = 180.0/np.pi
#  elevation mask
    emask = 10 
 
    ephemdata = myreadnav(navfilename)
    if len(ephemdata) == 0:
        print("empty ephmeris or does not exist")
        return
    #cartesian coordinates - from the header
    # number of epochs you want to read
    nep = 2
    recv,obs=readobs2(obsfilename,nep)
    print('A priori receiver coordinates', recv)
    if recv[0] == 0.0:
        print('using new a priori')
        recv[0]=-2715532
        recv[1] = -881995
        recv[2] = 5684286
        print('now using', recv)
    lat, lon, h = xyz2llh(recv,1e-8)
    print("%15.7f %15.7f"% (lat*r2d,  lon*r2d) )
    # zenith delay - meters
    zd = zenithdelay(h)
    u,east,north = up(lat,lon)
    epoch = np.unique(obs['TOW'])
    NN=len(epoch)
    NN = 2 # only do two positions
    for j in range(NN):
        print('Epoch', j+1)
        sats = obs['PRN'][np.where(epoch[j]==obs['TOW'])]
        gps_seconds = np.unique(obs['TOW'][np.where(epoch[j]==obs['TOW'])])
        gps_weeks = np.unique(obs['WEEK'][np.where(epoch[j]==obs['TOW'])])
        # not sure this does anything
        gps_weeks = gps_weeks.tolist()
        gps_seconds = gps_seconds.tolist()

        P1 = obs['C1'][np.where(epoch[j]==obs['TOW'])]
        P2 = obs['P2'][np.where(epoch[j]==obs['TOW'])]
#       do not know why this is here
#       S1 = obs['S1'][np.where(epoch[j]==obs['TOW'])]
        print('WEEK', gps_weeks, 'SECONDS', gps_seconds)
        k=0
        # set up the A matrix with empty stuff
        M=len(sats)
        print('Number of Satellites: ' , M)
        A=np.zeros((M,4))
        Y=np.zeros(M)
        elmask = np.zeros(M, dtype=bool)
        for prn in sats:
            closest = myfindephem(gps_weeks, gps_seconds, ephemdata, prn) #           
            p1=P1[k]
            p2=P2[k]
            p3 = ionofree(p1,p2)
            N=len(closest)
            if N > 0:    
                satv, relcorr = propagate(gps_weeks, gps_seconds, closest)
#               print("sat coor",gps_weeks,gps_seconds,satv)
                r=np.subtract(satv,recv)
                elea = elev_angle(u, r) # 
                tropocorr = zd/np.sin(elea)
                R,satv = mygeometric_range(gps_weeks, gps_seconds, prn, recv, closest)
                A[k]=-np.array([satv[0]-recv[0],satv[1]-recv[1],satv[2]-recv[2],R])/R
                elea = elev_angle(u,np.subtract(satv,recv))
                elmask[k] = elea*180/np.pi > emask
#               satellite clock correction
                satCorr = satclock(gps_weeks, gps_seconds, prn, closest)
                # prefit residual, ionosphere free pseudorange - geometric rnage
                # plus SatelliteClock - relativity and troposphere corrections
                Y[k] = p3-R+satCorr  -relcorr-tropocorr
#               print(int(prn), k,p1,R, p1-R)
                print(" {0:3.0f} {1:15.4f} {2:15.4f} {3:15.4f} {4:10.5f}".format(prn, p3, R, Y[k], 180*elea/np.pi))
                k +=1
# only vaguest notion of what is going on here - code from Ryan Hardy                
#       applying an elevation mask
        Y=np.matrix(Y[elmask]).T
        A=np.matrix(A[elmask])
        soln = np.array(np.linalg.inv(A.T*A)*A.T*Y).T[0]
#       update Cartesian coordinates
        newPos = recv+soln[:3]
        print('New Cartesian solution', newPos)
#       receiver clock solution
#        rec_clock = soln[-1]
        lat, lon, h = xyz2llhd(newPos)
        print("%15.7f %15.7f %12.4f "% (lat,  lon,h) )
# print("%15.5f"% xyz[0])
  
def readobs2(file,nepochX):
    """
    inputs: filename and number of epochs from the RINEX file you want to unpack
    returns: receiver Cartesian coordinates (meters) and observation blocks
    Kristine Larson, April 2017
    18aug20: updated so that 0,0,0 is returned if there are no receiver coordinates
    18aug21: version to include Glonass satellites etc
    """
    f = open(file, 'r')
    obs = f.read()
    f.close()
    testV = obs.split('APPROX POSITION XYZ')[0].split('\n')[-1].split()
#   set default receiver location values,              
    x0=0
    y0=0
    z0=0
    if len(testV) == 3:
        x0, y0, z0 = np.array(obs.split('APPROX POSITION XYZ')[0].split('\n')[-1].split(), 
							dtype=float)
    types = obs.split('# / TYPES OF OBSERV')[0].split('\n')[-1].split()[1:]
    print(types)
    tfirst =  obs.split('TIME OF FIRST OBS')[0].split('\n')[-1].split()
    print(tfirst)
    ntypes = len(types)
	#Identify unique epochs
    countstr = '\n '+tfirst[0][2:]+' '+'%2i' % np.int(tfirst[1])+' '+'%2i' % np.int(tfirst[2])
    print(countstr) # so this is 07 10 13, e.g.
    # figures out how many times this string appears between END OF HEADER and end of file
    nepoch = obs.split('END OF HEADER')[1].count(countstr)
    
    table = np.zeros((0, ntypes+3))
	#Identify kinds of observations in the file
    header = 'PRN\tWEEK\tTOW'
    t0 = 0
    for i in range(ntypes):
        header += '\t'+types[i]
    l = 0
#    print('countstr',countstr[0])
    print('header', header)
    epochstr_master = re.split(countstr, obs)
    print('number of data epochs in the file', nepoch)
    # only unpack a limited number of epochs    #nepoch = 5
    print('restricted number of epochs to be read ', nepochX)
    for i in range(nepochX):
#       clear satarray
        satarray = []
        epochstr = epochstr_master[i+1]        
        print(epochstr[21:23])
#       this is now the number of satellites properly read from the epochstr variable
        nprn = int(epochstr[21:23])
        print('Epochstr number of satellites',nprn)
        # this will only work with new rinex files that use G as descriptor
        for jj in range(nprn):
            i2 = 26+jj*3
            i1 = i2-3 
            satName = epochstr[i1:i2] 
            print(satName)
            if satName[0] == 'R':
#               found glonass
                list.append(satarray, 100+int(satName[1:3]))
            else:
#               assume rest are GPS for now
                list.append(satarray, int(satName[1:3]))
            print(satarray)
        prnstr = re.findall(r'G\d\d|G \d', epochstr.split('\n')[0])
        
        print(nprn, prnstr)
        if nprn > 12:
           # add the commant to read next line
           print('read next line of satellite names')
        # append
        table = np.append(table, np.zeros((nprn, ntypes+3)), axis=0)
        # decode year, month, day hour, minute ,second
        # convert 2 char to 4 char year
        # but it appears that the year and month and day 
        # are being hardwired from the header, which is both  very very odd
        # and stupid
        year = int(tfirst[0])
        month = int(tfirst[1])
        day = int(tfirst[2])
        hour = int(epochstr.split()[0])
        minute = int(epochstr.split()[1])
        second = float(epochstr.split()[2])
        print(year, month, day, hour, minute, second)
        week, tow = kgpsweek(year, month, day, hour, minute, second)
        print('reading week number ', week, 'sec of week', tow)
#       now you are saving the inforamtion to the table variable
        table[l:l+nprn, 1] = week
        table[l:l+nprn,2] = tow
        for j in range(nprn):
#            store this satellite
            print(j, satarray[j])
#            table[l+j, 0] = np.int(prnstr[j][1:])
# try this - using new definition of observed satellite
            table[l+j, 0] = satarray[j]
            # split up by end of line markers, which makes sense
#            print('something',2*j+1)
            line0 = epochstr.split('\n')[2*j+1]#.split()
            line = re.findall('.{%s}' % 16, line0+' '*(80-len(line0)))
# this seems to get strings of a certain width, which he will later change 
#through float# also where he made his mistake on reading last two characters 
#            print(j,len(line), line)

            for k in range(len(line)):
                if line[k].isspace():
                    table[l+j, 3+k] = np.nan
                    continue
                table[l+j, 3+k] = np.float(line[k][:-2])
#                print(table[l+j,3+k])
                # if more than 5 observaitons, he has to read the next line
            if ntypes > 5:
                line2 = epochstr.split('\n')[2*j+2].split()
#                print(line2)
                for k in range(len(line2)):
                    table[l+j, 3+len(line)+k] = np.float(line2[k][:-2])
        l += nprn	
    format = tuple(np.concatenate((('%02i', '%4i', '%11.7f'), 
							((0, ntypes+3)[-1]-3)*['%14.3f'])))
    # kinda strange - but ok - format statemenst for PRN, week, sec of week, etc
#    print(format)
# I guess it is making a super observable so that it can later be parsed.
# very strange 
    obs =  dict(zip(tuple(header.split('\t')), table.T))

    return np.array([x0,y0,z0]),obs

def get_ofac_hifac(elevAngles, cf, maxH, desiredPrec):
    """
    computes two factors - ofac and hifac - that are inputs to the
    Lomb-Scargle Periodogram code.
    We follow the terminology and discussion from Press et al. (1992)
    in their LSP algorithm description.

    INPUT
    elevAngles:  vector of satellite elevation angles in degrees 
    cf:(L-band wavelength/2 ) in meters    
    maxH:maximum LSP grid frequency in meters
    desiredPrec:  the LSP frequency grid spacing in meters
    i.e. how precise you want he LSP reflector height to be estimated
    OUTPUT
    ofac: oversampling factor
    hifac: high-frequency factor
    """
# in units of inverse meters
    X= np.sin(elevAngles*np.pi/180)/cf     

# number of observations
    N = len(X) 
# observing Window length (or span)
# units of inverse meters
    W = np.max(X) - np.min(X)         

# characteristic peak width, meters
    cpw= 1/W

# oversampling factor
    ofac = cpw/desiredPrec 

# Nyquist frequency if the N observed data samples were evenly spaced
# over the observing window span W, in meters
    fc = N/(2*W)

# Finally, the high-frequency factor is defined relative to fc
    hifac = maxH/fc  

    return ofac, hifac

def strip_compute(x,y,cf,maxH,desiredP,pfitV,minH):
    """
    strips snr data
    inputs; max reflector height, desiredP is desired precision in meters
    pfitV is polynomial fit order
    minH - do not allow LSP below this value
    returns 
    max reflector height and its amplitude
    min and max observed elevation angle
    riseSet is 1 for rise and -1 for set
    author: Kristine Larson
    """
    ofac,hifac = get_ofac_hifac(x,cf,maxH,desiredP)
#   min and max observed elevation angles
    eminObs = min(x); emaxObs = max(x)
    if x[0] > x[1]:
        riseSet = -1
    else:
        riseSet = 1

#   change so everything is rising, i.e. elevation angle is increasing
    ij = np.argsort(x)
    x = x[ij]
    y = y[ij]

    x = np.sin(x*np.pi/180)
#   polynomial fit done before

#   scale by wavelength
    x=x/cf
#    y=newy
#   get frequency spacing
    px = freq_out(x,ofac,hifac) 
#   compute spectrum using scipy
    scipy_LSP = spectral.lombscargle(x, y, 2*np.pi*px)

#   find biggest peak
#   scaling required to get amplitude spectrum
    pz = 2*np.sqrt(scipy_LSP/len(x))
#   now window
#    ij = np.argmax(px > minH)
#    new_px = px[ij]
    new_pz = pz[(px > minH)]
    new_px = px[(px > minH)]
    px = new_px
    pz = new_pz
#   find the max
    ij = np.argmax(pz)
    maxF = px[ij]
    maxAmp = np.max(pz)
    return maxF, maxAmp, eminObs, emaxObs,riseSet, px,pz

def window_data(s1,s2,s5,s6,s7,s8, sat,ele,azi,seconds,edot,f,az1,az2,e1,e2,satNu,pfitV,pele):
    """
    author kristine m. larson
    also calculates the scale factor for various GNNS frequencies.  currently
    returns meanTime in UTC hours and mean azimuth in degrees
    cf, which is the wavelength/2
    currently works for GPS, GLONASS, GALILEO, and Beidou
    new: pele are the elevation angle limits for the polynomial fit. these are appplied
    before you start windowing the data
    """
    cunit = 1
    dat = []; x=[]; y=[]
#   get scale factor
#   added glonass, 101 and 102
    if (f == 1) or (f==101) or (f==201):
        dat = s1
    if (f == 2) or (f == 20) or (f == 102) or (f==302):
        dat = s2
    if (f == 5) or (f==205):
        dat = s5
#   these are galileo frequencies (via RINEX definition)
    if (f == 206) or (f == 306):
        dat = s6
    if (f == 207) or (f == 307):
        dat = s7
    if (f == 208):
        dat = s8
#   get the scaling factor for this frequency and satellite number
#   print(f,satNu)
    cf = arc_scaleF(f,satNu)

#   if not, frequency does not exist, will be tripped by Nv
#   remove the direct signal component
    if (cf > 0):
        x,y,sat,azi,seconds,edot  = removeDC(dat, satNu, sat,ele, pele, azi,az1,az2,edot,seconds) 

#
    Nv = len(y); Nvv = 0 ; 
#   some defaults in case there are no data in this region
    meanTime = 0.0; avgAzim = 0.0; avgEdot = 1; Nvv = 0
    avgEdot_fit =1; delT = 0.0
#   no longer have to look for specific satellites. some minimum number of points required 
    if Nv > 30:
        model = np.polyfit(x,y,pfitV)
        fit = np.polyval(model,x)
#       redefine x and y as old variables
        ele = x
        dat = y - fit
#       ok - now figure out what is within the more restricted elevation angles
        x =   ele[(ele > e1) & (ele < e2) & (azi > az1) & (azi < az2)]
        y =   dat[(ele > e1) & (ele < e2) & (azi > az1) & (azi < az2)]
        ed = edot[(ele > e1) & (ele < e2) & (azi > az1) & (azi < az2)]
        a =   azi[(ele > e1) & (ele < e2) & (azi > az1) & (azi < az2)]
        t = seconds[(ele > e1) & (ele < e2) & (azi > az1) & (azi < az2)]
        sumval = np.sum(y)
        if sumval == 0:
            x = []; y=[] ; Nv = 0 ; Nvv = 0
#   since units were changed to volts/volts, the zeros got changed to 1 values
        if sumval == Nv:
            x = []; y=[] ; Nv = 0 ; Nvv = 0
        Nvv = len(y)
#       calculate average time in UTC (actually it is GPS time) in hours and average azimuth
#       this is fairly arbitrary, but can't be so small you can't fit a polymial to it
        if (Nvv > 10):
            dd = np.diff(t)
#           edot, in radians/sec
            model = np.polyfit(t,x*np.pi/180,1)
#  edot in radians/second
            avgEdot_fit = model[0]
            avgAzim = np.mean(a)
            meanTime = np.mean(t)/3600
            avgEdot = np.mean(ed) 
#  delta Time in minutes
            delT = (np.max(t) - np.min(t))/60 
# average tan(elev)
            cunit =np.mean(np.tan(np.pi*x/180))
#           return tan(e)/edot, in units of radians/hour now. used for RHdot correction
    outFact2 = cunit/(avgEdot_fit*3600) 
    outFact1 = cunit/(avgEdot*3600) 
    return x,y,Nvv,cf,meanTime,avgAzim,outFact1, outFact2, delT

def arc_scaleF(f,satNu):
    """
    input a frequency and put out a scale factor cf which is wavelength*0.5 
    """ 
#   default value for w so that if someone inputs an illegal frequency, it does not crash
    w = 0
    if f == 1:
        w = constants.wL1
    if (f == 2) or (f == 20):
        w = constants.wL2
    if f == 5:
        w = constants.wL5
#   galileo satellites
#   must be a smarter way to do this
    if (f > 200) and (f < 210):
        if (f == 201):
            w = constants.wgL1
        if (f == 205):
            w = constants.wgL5
        if (f == 206):
            w = constants.wgL6
        if (f == 207):
            w = constants.wgL7
        if (f == 208):
            w = constants.wgL8
#
#   add beidou 18oct15
    if (f > 300) and (f < 310):
        if (f == 302):
            w = constants.wbL2
        if (f == 307):
            w = constants.wbL7
        if (f == 306):
            w = constants.wbL6

#   glonass satellite frequencies
    if (f == 101) or (f == 102):
        w = glonass_channels(f,satNu) 
    cf = w/2
    return cf 

def freq_out(x,ofac,hifac):
    """
    inputs: x 
    ofac: oversamping factor
    hifac
    outputs: two sets of frequencies arrays
    """
#
# number of points in input array
    n=len(x)
#
# number of frequencies that will be used
    nout=np.int(0.5*ofac*hifac*n)
	 
    xmax = np.max(x) 
    xmin = np.min(x) 
    xdif=xmax-xmin 
# starting frequency 
    pnow=1.0/(xdif*ofac) 
    pstart = pnow
    pstop = hifac*n/(2*xdif)
# 
# output arrays
#    px = np.zeros(nout)
#    for i in range(0,nout):
#        px[i]=pnow
#        pnow=pnow+1.0/(ofac*xdif)
# simpler way
    pd = np.linspace(pstart, pstop, nout)
    return pd

def read_snr_file(obsfile):
    """
    input: observation filename
    output: contents of the file, withe various other metrics
    """

#SNR existance array : s0, s1,s2,s3,s4,s5,s6,s7,s8.  fields 0,3,4 are always false
    snrE = np.array([False, True, True,False,False,True,True,True,True],dtype = bool)
    f = np.genfromtxt(obsfile,comments='%')
    print('reading from a snr file ',obsfile)
    r,c = f.shape
    print('Number of rows:', r, ' Number of columns:',c)
#   store into new variable f
    sat = f[:,0]
    ele = f[:,1]
    azi = f[:,2]
    t =  f[:,3]
#   this is sometimes all zeros
    edot =  f[:,4]
#   looking for bad edot (since different rinex translators behave differently)
    median_edot = np.median(np.absolute(edot))
    s1 = f[:,6]
    s2 = f[:,7]
#   typically there is a zero in this row, but older files may have something
#   something that should not be used 
    s6 = f[:,5]


    s1 = np.power(10,(s1/20))  
    s2 = np.power(10,(s2/20))  
#
    s6 = s6/20
    s6 = np.power(10,s6)  
#
#   sometimes these records exist, sometimes not
#   depends on when the file was made, which version was used

    s5 = []
    s7 = []
    s8 = []
    if c > 8:
        s5 = f[:,8]
        if (sum(s5) > 0):
            s5 = s5/20; s5 = np.power(10,s5)  

    if c > 9:
        s7 = f[:,9]
        if (sum(s7) > 0):
            s7 = np.power(10,(s7/20))  
        else:
            s7 = []

    if c > 10:
        s8 = f[:,10]
        if (sum(s8) > 0):
            s8 = np.power(10,(s8/20))  
        else:
            s8 = []


    if (np.sum(s5) == 0):
        snrE[5] = False
#        print('no s5 data')
    if (np.sum(s6) == 0):
#        print('no s6 data')
        snrE[6] = False
    if (np.sum(s7) == 0):
#        print('no s7 data')
        snrE[7] = False
    if (np.sum(s8) == 0):
        snrE[8] = False
#        print('no s8 data')

#   now returned existence logical, snrE
    return sat, ele, azi, t, edot, s1, s2, s5, s6, s7, s8, median_edot, snrE

def find_satlist(f,snrExist):
    """
    inputs: frequency and boolean numpy array that tells you 
    if a signal is (potentially) legal
    outputs: list of satellites to use 
    author: kristine m. larson
    """
# set list of GPS satellites for now
# 
#   these are the only L2C satellites as of 18oct10
    l2c_sat = [1, 3, 5, 6, 7, 8, 9, 10, 12, 15, 17, 24, 25, 26, 27, 29, 30, 31, 32]
#   only L5 satellites thus far
    l5_sat = [1, 3,  6,  8, 9, 10, 24, 25, 26, 27, 30,  32]
#   assume l1 and l2 can be up to 32
    l1_sat = np.arange(1,33,1)
    satlist = []
    if f == 1:
        satlist = l1_sat
    if f == 20:
        satlist = l2c_sat
    if f == 2:
        satlist = l1_sat
    if f == 5:
        satlist = l5_sat
#   i do not think they have 26 - but ....
#   glonass L1
    if (f == 101) or (f==102):
# only have 24 frequencies defined
        satlist = np.arange(101,125,1)
#   galileo - 40 max?
#   use this to check for existence, mostly driven by whether there are 
#   extra columns (or if they are non-zero)
    gfs = int(f-200)

    if (f >  200) and (f < 210) and (snrExist[gfs]):
        satlist = np.arange(201,241,1)
#   galileo has no L2 frequency, so set that always to zero
    if f == 202:
        satlist = []
#   pretend there are 32 satellitesfor now
    if (f > 300):
        satlist = np.arange(301,333,1)

    if len(satlist) == 0:
        print('     illegal frequency: no sat list being returned')
    return satlist

def glonass_channels(f,prn):
    """
    inputs frequency and prn number
    returns wavelength for glonass satellite in meters
    logic from Simon Williams, matlab
    converted to function by KL, 2017 November
    converted to python by KL 2018 September
    """
#   we define glonass by adding 100.  remove this for definition of the wavelength
    if (prn > 100):
        prn = prn - 100
    lightSpeed = 299792458
    slot = [14,15,10,20,19,13,12,1,6,5,22,23,24,16,4,8,3,7,2,18,21,9,17,11]
    channel = [-7,0,-7,2,3,-2,-1,1,-4,1,-3,3,2,-1,6,6,5,5,-4,-3,4,-2,4,0]
    slot = np.matrix(slot)
    channel = np.matrix(channel)
#   main frequencies
    L1 = 1602e6
    L2 = 1246e6
#   deltas
    dL1 = 0.5625e6
    dL2 = 0.4375e6

    ch = channel[(slot == prn)]
    ch = int(ch)
#   wavelengths in meters
#    print(prn,ch,f)
    l = 0.0
    if (f == 101):
        l = lightSpeed/(L1 + ch*dL1)
    if (f == 102):
        l = lightSpeed/(L2 + ch*dL2)
    return l
def open_outputfile(station,year,doy):
    """
    inputs: station name, year, doy, and station name
    opens output file in REFL_CODE/year/results/station directory
    return fileID
    author kristine m. Larson
    if the results directory does not exist, it tries to make it. i think
    """
    fout = 0
#   primary reflector height output goes to this directory
    xdir = str(os.environ['REFL_CODE'])
    cdoy = '{:03d}'.format(doy)
#   extra file with rejected arcs
    frej=open('reject.txt','w+')
#   put a header in the file
    frej.write("%year, doy, maxF,sat,UTCtime, Azim, Amp,  eminO, emaxO,  Nv,freq,rise,Edot, PkNoise \n")
    filedir = xdir + '/' + str(year)  + '/results/' + station 
    filepath1 =  filedir + '/' + cdoy  + '.txt'
    print('output will go to:', filepath1)
    try:
        fout=open(filepath1,'w+')
#       put a header in the output file
        fout.write("%year, doy, maxF,sat,UTCtime, Azim, Amp,  eminO, emaxO,  Nv,freq,rise,EdotF, PkNoise  DelT     MJD \n")
        fout.write("% (1)  (2)   (3) (4)  (5)     (6)   (7)    (8)    (9)   (10) (11) (12) (13)  (14)     (15)     (16)\n")
        fout.write("%            m         hrs    deg   v/v    deg    deg                  hrs            min      \n")
    except:
        print('problem on first attempt - so try making results directory')
        cm = 'mkdir ' + xdir + '/' + str(year) + '/results/'
        os.system(cm)
        cm = 'mkdir ' + xdir + '/' + str(year) + '/results/' + station
        os.system(cm)
        try:
            fout=open(filepath1,'w+')
            print('successful open')
        except:
            print('problems opening the file')
            sys.exit()
    return fout, frej

def removeDC(dat,satNu, sat,ele, pele, azi,az1,az2,edot,seconds):
    """
#   remove direct signal using given elevation angle (pele) and azimuth 
    (az1,az2) constraints, return x,y as primary used data and windowed
    azimuth, time, and edot
#   removed zero points, which 10^0 have value 1.  used 5 to be sure?
    """
    p1 = pele[0]; p2 = pele[1]
#   look for data within these azimuth and elevation angle constraints
    x = ele[(sat == satNu) & (ele > p1) & (ele < p2) & (azi > az1) & (azi < az2) & (dat > 5)]
    y = dat[(sat == satNu) & (ele > p1) & (ele < p2) & (azi > az1) & (azi < az2) & (dat > 5)]
    edot = edot[(sat == satNu) & (ele > p1) & (ele < p2) & (azi > az1) & (azi < az2) & (dat > 5)]
    seconds = seconds[(sat == satNu) & (ele > p1) & (ele < p2) & (azi > az1) & (azi < az2) & (dat > 5)]
    azi = azi[(sat == satNu) & (ele > p1) & (ele < p2) & (azi > az1) & (azi < az2) & (dat > 5)]

    return x,y,sat,azi,seconds,edot

def quick_plot(plt_screen, gj,station,pltname):
    """
    inputs plt_screen variable (1 means go ahead) and integer variable gj
    which if > 0 there is something to plot
    also station name for the title
    pltname is png filename, if requested
    """
    if (plt_screen == 1  and gj > 0):
        plt.subplot(212)
        plt.xlabel('reflector height(m)')
        plt.ylabel('Spectral Amplitude')
        plt.subplot(211)
        plt.title(station)
        plt.ylabel('volts/volts')
        plt.xlabel('elev. angle(degrees)')
        if pltname != 'None':
            plt.savefig(pltname)
        else:
            print('plot file not saved ')
        plt.show()

def print_file_stats(ele,sat,s1,s2,s5,s6,s7,s8,e1,e2):
    """
    inputs 
    """
    gps = ele[(sat > 0) & (sat < 33) & (ele < e2) ]
    glonass = ele[(sat > 100) & (sat < 125) & (ele < e2) ]
    beidou = ele[(sat > 300) & (sat < 340) & (ele < e2) ]
    galileo = ele[(sat > 200) & (sat < 240) & (ele < e2) ]
    print('GPS     obs ', len(gps) )
    print('Glonass obs ', len(glonass))
    print('Galileo obs ', len(galileo))
    print('Beidou  obs ', len(beidou))

    return



def diffraction_correction(el_deg, temp=20.0, press=1013.25):
    """ Computes and return the elevation correction for refraction in the atmosphere.

    Computes and return the elevation correction for refraction in the atmosphere such that the elevation of the
    satellite plus the correction is the observed angle of incidence.

    Based on an empirical model by G.G. Bennet.
    This code was provided by Chalmers Group, Joakim Strandberg and Thomas Hobiger

    Parameters
    ----------
    el_deg : array_like
        A vector of true satellite elevations in degrees for which the correction is calculated.

    temp : float, optional
        Air temperature at ground level in degrees celsius, default 20 C.

    press : float, optional
        Air pressure at ground level in hPa, default 1013.25 hPa.

    Returns
    -------
    corr_el_deg : 1d-array
        The elevation correction in degrees.

    References:
    ----------
        Bennett, G. G. 'The calculation of astronomical refraction in marine navigation.'
        Journal of Navigation 35.02 (1982): 255-259.
    """
    el_deg = np.array(el_deg)

    corr_el_arc_min = 510/(9/5*temp + 492) * press/1010.16 * 1/np.tan(np.deg2rad(el_deg + 7.31/(el_deg + 4.4)))

    corr_el_deg = corr_el_arc_min/60

    return corr_el_deg
def mjd(y,m,d,hour,minute,second):
    """
    inputs: year, month, day, hour, minute,second
    output: modified julian day
    using information from http://infohost.nmt.edu/~shipman/soft/sidereal/ims/web/MJD-fromDatetime.html
    coded by kristine m. larson
    """
    if  (m <= 2):
        y, m = y-1, m+12
    if ((y, m, d) >= (1582, 10, 15)):
        A = int(y / 100)
        B = 2 - A + int(A / 4)
    else:
        B = 0
    C = int(365.25 * y)
    D = int(30.6001 *(m + 1))
    mjd = B + C + D + d - 679006
#   calculate seconds
    s = hour*3600 + minute*60 + second
    fracDay = s/86400
    return mjd, fracDay


def doy2ymd(year, doy):
    """
    inputs: year and day of year (doy)
    returns: some kind of datetime construct which can be used to get MM and DD
    """

    d = datetime.datetime(year, 1, 1) + datetime.timedelta(days=(doy-1))
    print('ymd',d)
    return d 

def getMJD(year,month,day,fract_hour):
    """
    inputs are year, month, day and fractional hour
    return is modified julian day (real8)
    """
#   convert fract_hour to HH MM SS
#   ignore fractional seconds for now
    hours = math.floor(fract_hour) 
    leftover = fract_hour - hours
    minutes = math.floor(leftover*60)
    seconds = math.floor(leftover*3600 - minutes*60)
#    print(fract_hour, hours, minutes, seconds)
    MJD, fracS = mjd(year,month,day,hours,minutes,seconds)
    MJD = MJD + fracS
    return MJD
def update_plot(plt_screen,x,y,px,pz):
    """
    input plt_screen integer value from gnssIR_lomb.
    (value of one means update the SNR and LSP plot)
    and values of the SNR data (x,y) and LSP (px,pz)
    """
    if (plt_screen == 1):
        plt.subplot(211)  
        plt.plot(x,y)
        plt.subplot(212)  
        plt.plot(px,pz)
def open_plot(plt_screen):
    """
    simple code to open a figure, called by gnssIR_lomb
    """
    if (plt_screen == 1):
        plt.figure()

def quick_rinex_snr(year, doy, station, option, orbtype):
    """
    inputs: year and day of year (integers) and station name
    option is for the snr creation
    orbtype can be nav or sp3.  if the former, then gpsSNR is used.
    if the later, then gnssSNR
    this assumes you follow my definitions for where things go,
    i.e. REFL_CODE and ORBITS
    """
    # define directory for the conversion executables
    exedir = os.environ['EXE']
    # FIRST, check to see if the SNR file already exists
    snrname_full = define_filename(station,year,doy,option)
    if (os.path.isfile(snrname_full) == True):
        print('snrfile already exists:', snrname_full)
    else:
    # SECOND MAKE SURE YOU HAVE THE ORBITS YOU NEED
        d = doy2ymd(year,doy); 
        month = d.month; day = d.day
        if orbtype == 'mgex':
            # this means you are using multi-GNSS and GFZ
            f,orbdir=getsp3file_mgex(year,month,day,'gbm')
            snrexe = exedir  + '/gnssSNR.e' 
        if orbtype == 'sp3':
            # this uses default IGS orbits, so only GPS
            f,orbdir=getsp3file(year,month,day)
            snrexe = exedir + '/gnssSNR.e' 
        if orbtype == 'gbm':
            # this uses GFZ multi-GNSS 
            f,orbdir=getsp3file_mgex(year,month,day,'gbm')
            snrexe = exedir + '/gnssSNR.e' 
        if orbtype == 'nav':
            f,orbdir=getnavfile(year, month, day) 
            snrexe = exedir  + '/gpsSNR.e' 
    # NOW MAKE SURE YOU HAVE THE RINEX FILE
        rinexfile,rinexfiled = rinex_name(station, year, month, day)
        print(rinexfile)
        if (os.path.isfile(rinexfile) == False):
            print('go get the rinex file')
            # new version
            rinex_unavco_obs(station, year, month, day) 
    # check to see if you found the rinex file
    # should check that the orbit really exists too
        oexist = os.path.isfile(orbdir + '/' + f) == True
        rexist = os.path.isfile(rinexfile) == True
        if (oexist and rexist):
            #convert to SNR file
            snrname = snr_name(station, year,month,day,option)
            orbfile = orbdir + '/' + f
            cmd = snrexe + ' ' + rinexfile + ' ' + snrname + ' ' + orbfile + ' ' + str(option)
            print(cmd); os.system(cmd)
            print('remove the rinexfile')
            os.system('rm -f ' + rinexfile)
#       move the snr file to its proper place
            if (os.stat(snrname).st_size == 0):
                print('you created a zero file size which could mean a lot of things')
                print('bad exe, bad snr option, do not really have the orbit file')
                os.system('rm -f ' + snrname)
            else:
                store_snrfile(snrname,year,station) 
        else:
            print('rinex file or orbit file does not exist, so there is nothing to convert')

def store_orbitfile(filename,year,orbtype):
    """
    simple code to move an orbit file to the right place 
    inputs are the filename, the year, and the kind of orbit
    (sp3 or nav)
    """
    xdir = str(os.environ['ORBITS']) + '/' + str(year)
    # check that directories exist
    if not os.path.isdir(xdir): #if year folder doesn't exist, make it
        os.makedirs(xdir)
    xdir = str(os.environ['ORBITS']) + '/' + str(year) + '/' + orbtype
    if not os.path.isdir(xdir): #if year folder doesn't exist, make it
        os.makedirs(xdir)
    if (os.path.isfile(filename) == True):
        cmd = 'mv ' + filename + ' ' + xdir 
        print('moving ', filename, ' to ', xdir)
        os.system(cmd)
    else:
        print('file did not exist, so it was not stored')
    return xdir


def store_snrfile(filename,year,station):
    """
    simple code to move an snr file to the right place 
    inputs are the filename, the year, and the station name
    """
    xdir = str(os.environ['REFL_CODE']) + '/' + str(year)
    # check that directories exist
    if not os.path.isdir(xdir): #if year folder doesn't exist, make it
        os.makedirs(xdir)
    xdir = xdir + '/snr'
    if not os.path.isdir(xdir): #if year folder doesn't exist, make it
        os.makedirs(xdir)
    xdir = xdir + '/' + station 
    if not os.path.isdir(xdir): #if year folder doesn't exist, make it
        os.makedirs(xdir)
    if (os.path.isfile(filename) == True):
        cmd = 'mv ' + filename + ' ' + xdir 
        print(cmd)
        os.system(cmd)
    else:
        print('file does not exist, so nothing was moved')

def rinex_name(station, year, month, day):
    """
    author: kristine larson
    given station (4 char), year, month, day, return rinexfile name
    and the hatanaka equivalent
    """
    doy,cdoy,cyyyy,cyy = ymd2doy(year,month,day)

    fnameo = station + cdoy + '0.' + cyy + 'o'
    fnamed = station + cdoy + '0.' + cyy + 'd'
    return fnameo, fnamed

def snr_name(station, year, month, day,option):
    """
    author: kristine larson
    given station (4 char), year, month, day, and snr option,
    return snr filename (and directory) using my system
    """
    doy,cdoy,cyyy,cyy = ymd2doy(year,month,day)

    fname = station + cdoy + '0.' + cyy + '.snr' + str(option)
    return fname

def nav_name(year, month, day):
    """
    kristine m. larson
    inputs are year month and day
    returns nav file name and directory
    """
    doy,cdoy,cyyyy,cyy = ymd2doy(year,month,day)
    navfilename = 'auto'  + cdoy + '0.' + cyy  +  'n'
    navfiledir = str(os.environ['ORBITS']) + '/' + cyyyy + '/nav'
    return navfilename,navfiledir

def sp3_name(year,month,day,pCtr):
    """
    kristine m. larson
    inputs are year month and day and processing center
    returns sp3 file name and directory
    """
    name,clkn=igsname(year,month,day)
    gps_week = name[3:7]
    sp3name = pCtr + name[3:8] + '.sp3'
    sp3dir = str(os.environ['ORBITS']) + '/' + str(year) + '/sp3'
    return sp3name, sp3dir

def rinex_unavco_obs(station, year, month, day):
    """
    author: kristine larson
    picks up a RINEX file from unavco.  
    normal observation file - not Hatanaka
    new version i was testing out
    """
    doy,cdoy,cyyyy,cyy = ymd2doy(year,month,day)
    rinexfile,rinexfiled = rinex_name(station, year, month, day) 
    unavco= 'ftp://data-out.unavco.org' 
    filename = rinexfile + '.Z'
    url = unavco+ '/pub/rinex/obs/' + cyyyy + '/' + cdoy + '/' + filename
    print(url)
    try:
        wget.download(url,filename)
        cmd = 'uncompress ' + filename
        os.system(cmd) 
    except:
        print('some kind of problem with download',rinexfile)


