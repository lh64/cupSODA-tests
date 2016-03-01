import pylab as plt
import numpy as np

cupsoda_data = np.genfromtxt('earm_new_simulate', delimiter=',', dtype=None, names=True)
cupsoda_data2 = np.genfromtxt('earm_slim', delimiter=',', dtype=None, names=True)

plt.plot(cupsoda_data['nsims'],cupsoda_data['cupsodatime'],'r-x',label='new simulate')
plt.plot(cupsoda_data['nsims'],cupsoda_data['pythontime'],'r-o')
plt.plot(cupsoda_data2['nsims'],cupsoda_data2['cupsodatime'],'b-x',label='nvme')
plt.plot(cupsoda_data2['nsims'],cupsoda_data2['pythontime'],'b-o')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc=0)
plt.savefig('earm_on_diablo.png')
plt.show()
quit()

plt.subplot(212)
plt.plot(cupsoda_data['nsims'],cupsoda_data['pythontime']-cupsoda_data['cupsodatime'])
#plt.xscale('log')
#plt.yscale('log')
plt.show()