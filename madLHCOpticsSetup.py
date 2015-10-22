empty 		= ()
base 		= '/afs/cern.ch/eng/lhc/optics/V6.503/'
scenario 	= 'runII/2015/'

# BEAM 1 DEFAULT SETTINGS

beamDefaultLib = {'SEQUENCE' : "LHCB1", 'BV':1,'ENERGY':7000.,'PARTICLE':"PROTON",'KBUNCH':2808,\
			'NPART':1.15E11,'EX':5.02645721478119E-10,'EY':5.02645721478119E-10,\
			'ET':8.54E-6,'SIGE':0.000111,'SIGT':0.077
			}

# CONSTRUCTING THE MADBEAM EXPRESSIONS

def MADBeam(beam="Beam1",inputDict=beamDefaultLib):
	output = beam + ': '
	for k, v in inputDict.iteritems():
		output += str(k) + ' = ' + str(v) + ', '
	output = output[:-2]+ ";"
	return output


print MADBeam()

beam1 = MADBeam()
beam2Lib = beamDefaultLib
beam2Lib['SEQUENCE']="LHCB2"
beam2 = MADBeam("Beam2",beam2Lib)

print MADBeam("Beam2",beam2Lib)

madLHCOpticsSetupDict = {'madLHCOpticsDatabaseStem' : base,\
			'madLHCOpticsScenario'     :  scenario,\
			'madLHCOpticsProlog'       : empty,\
			'madLHCOpticsEpilog'       : empty,\
			'madLHCSequenceFiles'	   : base + scenario + 'lhc_as-built.seq',\
			'madLHCSettingsFiles'      : base + scenario + 'opt_800_10000_800_3000.madx',\
			'madLHCSettingsRules'      : {},\
			'madLHCBeams'              : {"Beam1":beam1,"Beam2":beam2}
 			}

def madLHCOpticsSetup(madxinputinputDict=madLHCOpticsSetupDict):
	
	stem = inputDict['madLHCOpticsDatabaseStem'] 
	scenario = inputDict['madLHCOpticsScenario']
	
