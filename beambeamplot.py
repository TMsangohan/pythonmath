import numpy as nm
import matplotlib.pyplot as plt
import scipy.constants as cons 


n 	= 10
betarel = 1
sigma 	= 5

r = nm.arange(-10.,10,0.01)

bbforce = - n * cons.e**2 / (2 * cons.pi * cons.epsilon_0 * r ) * (1 - nm.exp(-r**2/sigma))

plt.plot(r,- n * cons.e**2 / (2 * cons.pi * cons.epsilon_0 * r )  * (1-nm.exp(-r**2/5)),'r',r,- n * cons.e**2 / (2 * cons.pi * cons.epsilon_0 * r )  * (1-nm.exp(-r**2/10)),'b',r,- n * cons.e**2 / (2 * cons.pi * cons.epsilon_0 * r )  * (1-nm.exp(-r**2/2)),'g')
plt.show()
