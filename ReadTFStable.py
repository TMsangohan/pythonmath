import cern_pymad_domain_tfs as dom
import cern_pymad_io_tfs as io
import numpy as nm
import matplotlib.pyplot as plt
import math

data = io.tfsDict("LHCB1-Optics.tfs")

keys = []

for key, val in data[0].items():
    keys.append(key)

print keys
print data[1]

#plt.plot(data[0]['s'], data[0]['x'], 'b-')
plt.axis([-500, 500, -0.01, 0.01])
#plt.show()

a = nm.array(data[0]['betx'])
b= data[1]['ex']
print b * a
aa=nm.array(data[0]['x'])+nm.sqrt(b* a + nm.array(data[0]['dx'])**2 * data[1]['sige']**2)
aaa=nm.array(data[0]['x'])-nm.sqrt(b* a + nm.array(data[0]['dx'])**2 * data[1]['sige']**2)

plt.plot(data[0]['s'], data[0]['x'], 'b-',data[0]['s'], aa, 'r-',data[0]['s'], aaa, 'r-')

plt.show()