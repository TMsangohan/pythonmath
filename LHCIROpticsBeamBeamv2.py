# Author 	: Tom Mertens
# Date   	: 23/09/2015
# Version	: 1.0
# Version	: 1.1 (22/10/2015) - updated to use OrderedDict from collections
# Partial implementation of LHCIROpticsBeamBeam function in mathamatica developed by John Jowett CERN

import cern_pymad_domain_tfs as dom
import cern_pymad_io_tfs as io
import numpy as nm
import matplotlib.pyplot as plt
import math
from scipy.interpolate import interp1d
from collections import OrderedDict


# DEFINING CONVERSION CONSTANTS
r0m 	 = 2.5669694e-45
cns	 = 0.2999792458000000004


# converting python float to mathematica notation
def float_to_mathematica(x):
    return ("%e"%x).replace('e','*10^')

# getting the fractional part of a python float
def get_frac(n):
        frac, whole = math.modf(n)
        return frac

# python version of the mathematica Less function
def Less(x,y):
        if x < y:
                return True
        else:
                return False

# converting a pyhton dictionary to an mfs object to import in mathematica

def TableToMfs(table,tfssummary):
	mfs = 'mfs[{'
	for k in tfssummary:
		if isinstance(tfssummary[k],str):
			addval = '"' + str(tfssummary[k]).upper() + '"'
		elif isinstance(tfssummary[k],float):
			addval = float_to_mathematica(tfssummary[k])
		else:
			addval = str(tfssummary[k]).upper()
		mfs += '{"'+ str(k).upper() + '",'+  addval  + '},'
	mfs += '},{'
	for k in table:
		if str(k)=='x':
			mfs += '"XC",'
		elif str(k)=='y':
			mfs += '"YC",'
		elif str(k)=='px':
			mfs += '"PXC",'
		elif str(k)=='py':
			mfs += '"PYC",'
		else:
			mfs += '"'+ str(k).upper() + '",'
	mfs = mfs[:-1]
	mfs += '},{'
	for i in range(len(table['s'])):
		mfs += '{'
		for k in table:
			if isinstance(table.get(k)[i],str):
				addvals = '"' + str(table.get(k)[i]).upper() + '"'
			elif isinstance(table.get(k)[i],float):
				addvals = float_to_mathematica(table.get(k)[i])
			else:
				addvals = str(table.get(k)[i]).upper()
			mfs += addvals + ','
		mfs = mfs[:-1]
		mfs += '},'
	mfs = mfs[:-1]
	mfs += '}]'
	return mfs


# Reading in the tfs file to python library 

datab1 	= io.tfsDict("LHCB1-Optics.tfs")
datab2  = io.tfsDict("LHCB2-Optics.tfs")  
rb1 	= datab1[0]
rb2 	= datab2[0]
cols1   = ['name','keyword','y','px','py','betx','bety','alfx','alfy','mux','muy','dx','dy','dpx','dpy','s','x']
cols2   = ['y','betx','bety','dx','dy','s','x']
gamma1  = datab1[1]['gamma']
ex1	= datab1[1]['ex']
ey1 	= datab1[1]['ey']
sige1	= datab1[1]['sige']
n1	= datab1[1]['npart']
q1	= datab1[1]['charge']
m1  	= datab1[1]['mass']
s1	= nm.array(rb1.get('s'))
xc1	= nm.array(rb1.get('x'))
yc1	= nm.array(rb1.get('y'))
bx1	= nm.array(rb1.get('betx'))
by1	= nm.array(rb1.get('bety'))
dx1	= nm.array(rb1.get('dx'))
dy1	= nm.array(rb1.get('dy'))
ex2	= datab2[1]['ex']
ey2	= datab2[1]['ey']
sige2	= datab2[1]['sige']
n2	= datab2[1]['npart']
q2	= datab2[1]['charge']
s2	= nm.array(rb2.get('s'))
sc	= s1 / cns

dpFactor = 2 * q1 * (q2 * n2) *  r0m * 5.60959e26 / (m1 * gamma1)
dqFactor = dpFactor / (4 * math.pi) 

# interpolations
xc2 	= interp1d(s2,nm.array(rb2.get('x')))(s1)
yc2	= interp1d(s2,nm.array(rb2.get('y')))(s1)
bx2	= interp1d(s2,nm.array(rb2.get('betx')))(s1)
by2	= interp1d(s2,nm.array(rb2.get('bety')))(s1)
dx2	= interp1d(s2,nm.array(rb2.get('dx')))(s1)
dy2	= interp1d(s2,nm.array(rb2.get('dy')))(s1)

sx1	= nm.sqrt(bx1 * ex1 + dx1**2 * sige1**2)
sx2	= nm.sqrt(bx2 * ex2 + dx2**2 * sige2**2)
sy1	= nm.sqrt(by1 * ey1 + dy1**2 * sige1**2)
sy2	= nm.sqrt(by2 * ey2 + dy2**2 * sige2**2)
sepx	= xc1 - xc2
sepy	= yc1 - yc2
sepr	= nm.sqrt(sepx**2 + sepy**2)

sepxsx1 = map(math.fabs,sepx/sx1)
sepxsx2 = map(math.fabs,sepx/sx2)
sepxsy1 = map(math.fabs,sepx/sy1)
sepxsy2 = map(math.fabs,sepx/sy2)
sepysx1 = map(math.fabs,sepy/sx1)
sepysx2 = map(math.fabs,sepy/sx2)
sepysy1 = map(math.fabs,sepy/sy1)
sepysy2 = map(math.fabs,sepy/sy2)

seprxy1 = sepr / map(max,nm.vstack((sx1,sy1)).T)
seprxy2 = sepr / map(max,nm.vstack((sx2,sy2)).T)

bbdpx	= dpFactor * sepx / sepr**2
bbdpy	= dpFactor * sepy / sepr**2
bbdqx	= dqFactor * ( bx1 * (sepx**2 - sepy**2) ) / sepr**4
bbdqy	= dqFactor * ( -by1 * (sepx**2 - sepy**2) ) / sepr**4

# GENERATING THE DICTIONARY
rb1['cS']		= sc
rb1['SIGX']		= sx1
rb1['SEPX']     	= sepx
rb1['SEPXSIGX1']	= sepxsx1
rb1['SEPXSIGX2']	= sepxsx2
rb1['SEPXSIGY1']	= sepxsy1
rb1['SEPXSIGY2']	= sepxsy2
rb1['SIGY']     	= sy1
rb1['SEPY']     	= sepy
rb1['SEPYSIGX1']	= sepysx1
rb1['SEPYSIGX2']	= sepysx2
rb1['SEPYSIGY1']	= sepysy1
rb1['SEPYSIGY2']	= sepysy2
rb1['SEPR']     	= sepr
rb1['SEPRSIGXY1']     	= seprxy1
rb1['SEPRSIGXY2']     	= seprxy2
rb1['SIGX2']    	= sx2
rb1['SIGY2']     	= sy2
rb1['BBDPX']     	= bbdpx
rb1['BBDPY']     	= bbdpy
rb1['BBDQX']     	= bbdqx
rb1['BBDQY']     	= bbdqy

# UPDATING THE KEY VALUES
bunchch = n2 * q2

# WRITING OUTPUT TO FILE
datab1[1]['BUNCHCHARGE2'] = str(bunchch)
file = open("BeamBeamB1.m", "w")
file.write(TableToMfs(rb1,datab1[1]))
file.close()
