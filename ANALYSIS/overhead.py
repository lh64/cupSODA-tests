import numpy as np
import pylab as plt
import os
from matplotlib.ticker import ScalarFormatter 

figs = os.path.join('..','FIGS')
if not os.path.exists(figs):
    os.makedirs(figs)

datafile = os.path.join('..', 'diablo-GPU_timing.csv')
cupsoda_data = np.genfromtxt(datafile, delimiter=',', dtype=None, names=True)

print cupsoda_data.dtype.names

colors = ['blue', 'red', 'green']

plt.figure('overhead')

xmax = 0
for i,model in enumerate(['tyson', 'ras', 'earm']):
    
    xdata = []
    for x in [d['nsims'] for d in cupsoda_data if d['model'] == model]:
        if x not in xdata:
            xdata.append(x)

    ydata = []
    for x in xdata:
        cs_times = [d['cupsodatime'] for d in cupsoda_data 
                      if  d['model'] == model
                      and d['nsims'] == x]
        py_times = [d['pythontime'] for d in cupsoda_data 
                      if  d['model'] == model
                      and d['nsims'] == x]
        print cs_times
        min_index = 0
        for j in range(1,len(cs_times)):
            if cs_times[j] < cs_times[j-1]:
                min_index = j
        ydata.append((py_times[min_index]-cs_times[min_index])/cs_times[min_index]*100)
        
    plt.plot(xdata, ydata, 'o-', ms=12, lw=3, mec=colors[i], color=colors[i], label=model) #, 'o', ms=12, mew=1.5, mfc='none')

    plt.legend(loc=0)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('# of simulations')
    plt.ylabel('(py_time-cs_time)/cs_time*100')
    if xmax < max(xdata)+0.25*max(xdata):
        xmax = max(xdata)+0.25*max(xdata)
        plt.xlim(xmax=xmax)
    plt.subplot(111).yaxis.set_major_formatter(ScalarFormatter())
    
    plt.savefig(os.path.join(figs,'overhead.pdf'))

plt.show()

