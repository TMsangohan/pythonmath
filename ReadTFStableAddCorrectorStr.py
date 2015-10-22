# Author 	: Tom Mertens
# Date   	: 23/09/2015
# Version	: 1.0
# Partial implementation of AddCorrectorStrenght function in mathamatica developed by John Jowett CERN

import cern_pymad_domain_tfs as dom
import cern_pymad_io_tfs as io
import numpy as nm
import matplotlib.pyplot as plt
import math

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

data 	= io.tfsDict("LHCB1-Optics.tfs") 
r 	= data[0]

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
