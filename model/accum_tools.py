import numpy as np 
import datetime
from scipy.interpolate import interp1d
import math
import gps as g



def snr_RH(yr, doy, st, step=1, **kwargs):
    """
    author: Andrew Hoffman
    reads snr file and calls writeRH with nan value or median reflector height.
    """

    d = g.doy2ymd(year,doy)
    obsfile = g.snr_name(st,yr,month,day)
    try:      
        f = open(obsfile)
        f.close()
    	sat, ele, azi, t, edot, s1, s2, s5, s6, s7, s8, median_edot, snrE=read_snr_file(obsfile) 
    	if step==1:
    		RH = median_edot
    	else:
            #sort(t)
    		RH = np.median(edot)
    	writeRH(RH,st,yr,doy,step)
    except:
    	print('sorry - the snr file does not exist')
    	RH ='nan'
        writeRH(RH,st,yr,doy,step)
    return RH



def writeRH(RH,station,year,doy,step):
	"""
	author: Andrew Hoffman
	writes the RH to file that will serve as input for accumulation algorithm
	"""
	fout = 0
	xdir = str(os.enviorn['REFL_CODE'])
	cdoy = '{:03d}'.format(doy)
#   extra file with rejected arcs
    frej=open('reject.txt','w+')
#   put a header in the file
    frej.write("%year, doy, medianRH")
    filedir = xdir + '/' + str(year)  + '/results/' + station 
    filepath1 =  filedir + '/' + cdoy  + '.txt'
    print('output will go to:', filepath1)
    try:
        fout=open(filepath1,'w+')
#       put a header in the output file
        fout.write("%year, doy RH \n")
        fout.write("% (1)  (2)   (3)\n")
        fout.write("%  yr         doy  m  \n")
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






def find_breaks(yr1,doy1,yr2,doy2,st,**kwargs):
    """
    author: Andrew Hoffman
    finds breaks in the RH record
    """
    



def RH_accum(st,breaks_file,**kwargs):
    """
    author: Andrew Hoffman
    reads snr file and calls writeRH with nan value or median reflector height.
    """

    d = g.doy2ymd(year,doy)
    obsfile = g.snr_name(st,yr,month,day)
    try:      
        f = open(obsfile)
        f.close()
        sat, ele, azi, t, edot, s1, s2, s5, s6, s7, s8, median_edot, snrE=read_snr_file(obsfile) 
        if step==1:
            RH = median_edot
        else:
            #sort(t)
            RH = np.median(edot)
        writeRH(RH,st,yr,doy,step)
    except:
        print('sorry - the snr file does not exist')
        RH ='nan'
        writeRH(RH,st,yr,doy,step)
    return RH

