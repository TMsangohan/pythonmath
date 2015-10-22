# Author 	: Tom Mertens
# Date   	: 23/09/2015
# Version	: 1.0
# Partial implementation of AddQuadrupoleStrenght function in mathamatica developed by John Jowett CERN

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

K1L		= []
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
        'I1in/A','dBdxMIN','dBdx','dBdxMAX','K1MIN','K1','K1MAX',
        'K1PERCENTMIN','K1PERCENTMAX','K1CHECK'
]


# END INITIALIZATION OF VARIABLES


for i in range(len(r['s'])):
	if kws[i]=='QUADRUPOLE':
		K1L.append(r.get('k1l')[i])
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
K1L	= nm.array(K1L)
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

calibcol     = CALIB
polcol       = POLARITY
k1mincol0    = KMIN / brho
k1maxcol0    = KMAX / brho
k1mincol     = 0.5 * (polcol + 1) * k1mincol0 + 0.5 * (polcol - 1) * k1maxcol0
k1maxcol     = 0.5 * (polcol + 1) * k1maxcol0 + 0.5 * (polcol - 1) * k1mincol0
k1col	     = K1L / L

kickcheckcol = []
for (x,y) in nm.vstack((k1mincol,k1maxcol)).T:
	kickcheckcol.append(Less(x,y))
	
Icol = polcol * brho * k1col / calibcol

#print r.get('name')[secondtab[1]]
#print r.get('k1l')[secondtab[1]]
#print K1L[1]
#print KMAX[1]
#print k1col[1]
#print k1maxcol0[1]
#print brho
#print k1maxcol[1]
#print len(k1col), len(KMAX)
#print 100. * k1col[1]
#print 100. * k1col[1]/ KMAX[1]
#print (100. * k1col[1]) / KMAX[1]

percentmaxcol	= 100. * k1col / k1maxcol
percentmincol 	= 100. * k1col / k1mincol

#print percentmaxcol[1]
#print percentmincol[1]

dBdxcol 	= brho * k1col
dBdxmincol 	= brho * k1mincol
dBdxmaxcol	= brho * k1maxcol

AddQuadColumnKeys = [
	'I1in/A','dBdxMIN','dBdx','dBdxMAX','K1MIN','K1','K1MAX',
	'K1PERCENTMIN','K1PERCENTMAX','K1CHECK'
]
AddQuadColumnValues = [
	Icol,dBdxmincol,dBdxcol,dBdxmaxcol,
	k1mincol,k1col,k1maxcol,percentmincol,percentmaxcol,kickcheckcol]

# CREATING A NEW DICTIONARY CONTAINING THE NEW COLUMNS
DictQuad = zip(AddQuadColumnKeys,AddQuadColumnValues)

# MERGING THE OLD DICTIONARY AND THE ONE CONTAINING THE NEW COLUMNS
Newr.update(DictQuad)

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

file = open("LCHAddQuadrupoleStrengths.m", "w")
file.write(TableToMfs(Newr,data[1],brho,orderedlist))
file.close()
