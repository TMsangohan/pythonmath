# Author 	: Tom Mertens
# Date   	: 23/09/2015
# Version	: 1.0
# Partial implementation of LHCIROpticsBeamBeam function in mathamatica developed by John Jowett CERN

import cern_pymad_domain_tfs as dom
import cern_pymad_io_tfs as io
import numpy as nm
import matplotlib.pyplot as plt
import math
from scipy.interpolate import interp1d

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

def TableToMfs(table,tfssummary,brho,orderedlist):
	mfs = 'mfs[{'
	for k in tfssummary:
		if isinstance(tfssummary[k],str):
			addval = '"' + str(tfssummary[k]).upper() + '"'
		elif isinstance(tfssummary[k],float):
			addval = float_to_mathematica(tfssummary[k])
		else:
			addval = str(tfssummary[k]).upper()
		mfs += '{"'+ str(k).upper() + '",'+  addval  + '},'
	mfs += '{"brho",'+ str(brho) + '}'
	mfs += '},{'
	for k in orderedlist:
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
		for k in orderedlist:
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

print m1,q1,q2,n2

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

print sepx, sepr

# INITIALIZATION OF VARIABLES FOR ADDING CORRECTOR STRENGTHS

L		= []
KW		= []
LRAD		= []
MUX		= []
MUY		= []
HKICK		= []
VKICK		= []
KMAX		= []
KMIN		= []
POLARITY	= []
CALIB		= []
secondtab	= []
kws 		= r.get('keyword')
Newr		= {}
orderedlist 	= [
        'name','keyword','parent','s','l','lrad','kick','hkick',
        'vkick','angle','k0l','k1l','k2l','k3l','x','y','px','py',
        'betx','bety','alfx','alfy','mux','muy','dx','dy','dpx','dpy',
        'kmax','kmin','calib','polarity','apertype','aper_1','n1',
        'MUXfrac','MUYfrac','HKICK','VKICK','KICKMIN','KICKMAX',
        'KICKPERCENTMIN','KICKPERCENTMAX','KICKCHECK','B/T','I/A'
]


# END INITIALIZATION OF VARIABLES


for i in range(len(r['s'])):
	if kws[i]=='HKICKER' or kws[i]=='VKICKER' or kws[i]=='MARKER' :
		L.append(r.get('l')[i])
		KW.append(r.get('keyword')[i])
		LRAD.append(r.get('lrad')[i])
		MUX.append(r.get('mux')[i])
		MUY.append(r.get('muy')[i])
		HKICK.append(r.get('hkick')[i])
		VKICK.append(r.get('vkick')[i])
		KMAX.append(r.get('kmax')[i])
		KMIN.append(r.get('kmin')[i])
		POLARITY.append(r.get('polarity')[i])
		CALIB.append(r.get('calib')[i])
		secondtab.append(i)

for kk in r:
	temp = [r[kk][i] for i in secondtab]
	Newr.update({kk:temp})

# TURNING LISTS INTO ARRAYS
L 	= nm.array(L)
KW	= nm.array(KW)
LRAD	= nm.array(LRAD)
MUX	= nm.array(MUX)
MUY	= nm.array(MUY)
HKICK	= nm.array(HKICK)
VKICK	= nm.array(VKICK)
KMAX	= nm.array(KMAX)
KMIN	= nm.array(KMIN)
POLARITY= nm.array(POLARITY)
CALIB	= nm.array(CALIB)

brho = 3.33564 * data[1]['pc']/data[1]['charge']

lcol         = map(max,nm.vstack((L,LRAD)).T)
muxfrac      = map(get_frac,MUX)
muyfrac      = map(get_frac,MUY)
calibcol     = CALIB
polcol       = POLARITY
kickmincol0  = (KMIN * lcol)/ brho
kickmaxcol0  = (KMAX * lcol)/ brho
kickmincol   = 0.5 * (polcol + 1) * kickmincol0 + 0.5 * (polcol - 1) * kickmaxcol0
kickmaxcol   = 0.5 * (polcol + 1) * kickmaxcol0 + 0.5 * (polcol - 1) * kickmincol0
kickhcol     = HKICK
kickvcol     = VKICK

kickcheckcol = []
for (x,y) in nm.vstack((kickmincol,kickmaxcol)).T:
	kickcheckcol.append(Less(x,y))
	
Bcol = []
for i in range(len(KW)):
	if KW[i]=='HKICKER':
		Bcol.append(kickhcol[i] / lcol[i])
	elif KW[i]=='VKICKER':
		Bcol.append(kickvcol[i] / lcol[i])
	else:
		Bcol.append(0.0)
Bcol 	= brho * nm.array(Bcol)

Icol = []
for i in range(len(KW)):
	if calibcol[i] < 1.0e-10:
		Icol.append(0.0)
	else:
		Icol.append(Bcol[i] / calibcol[i])

#print kickvcol[300]
#print lcol[300]
#print Bcol[300]
#print KMAX[300]

percentmaxcol = []
percentmincol = []

for i in range(len(KW)):
	if math.fabs(KMAX[i]) < 1.0e-10 :
		percentmaxcol.append(0.0)
	else:
		percentmaxcol.append(100. * Bcol[i] / KMAX[i])
	if math.fabs(KMIN[i]) < 1.0e-10:
		percentmincol.append(0.0)
	else:
		percentmincol.append(100. * Bcol[i] / KMIN[i])

AddCorrColumnKeys = ['MUXfrac','MUYfrac','HKICK','VKICK','KICKMIN','KICKMAX','KICKPERCENTMIN','KICKPERCENTMAX','KICKCHECK','B/T','I/A']
AddCorrColumnValues = [muxfrac,muyfrac,kickhcol,kickvcol,kickmincol,kickmaxcol,percentmincol,percentmaxcol,kickcheckcol,Bcol,Icol]

# CREATING A NEW DICTIONARY CONTAINING THE NEW COLUMNS
Dictcorr = zip(AddCorrColumnKeys,AddCorrColumnValues)

# MERGING THE OLD DICTIONARY AND THE ONE CONTAINING THE NEW COLUMNS
Newr.update(Dictcorr)

#plt.plot(data[0]['s'], data[0]['x'], 'b-')
#plt.axis([11000, 15000, -0.002, 0.002])
#plt.show()

#a = nm.array(data[0]['betx'])
#b= data[1]['ex']
#print b * a
#aa=nm.array(data[0]['x'])+nm.sqrt(b* a + nm.array(data[0]['dx'])**2 * data[1]['sige']**2)
#aaa=nm.array(data[0]['x'])-nm.sqrt(b* a + nm.array(data[0]['dx'])**2 * data[1]['sige']**2)

#plt.plot(data[0]['s'], data[0]['x'], 'b-',data[0]['s'], aa, 'r-',data[0]['s'], aaa, 'r-')

#print Newr['KICKPERCENTMAX']
#plt.plot(Newr['s'], Newr['KICKPERCENTMAX'],'b')
#plt.show()

#print TableToMfs(Newr,data[1],brho,orderedlist)

# SAVING THE NEW MFS OBJECT TO FILE

file = open("LCHAddCorrectorStrengths.m", "w")
file.write(TableToMfs(Newr,data[1],brho,orderedlist))
file.close()
