import numpy as np
import pylab as plt
import os

figs = os.path.join('..','FIGS')
if not os.path.exists(figs):
    os.makedirs(figs)

datafile = os.path.join('..', 'cupsoda_timings_all.csv')
cupsoda_data = np.genfromtxt(datafile, delimiter=',', dtype=None, names=True)

datafile = os.path.join('..', 'scipy_timings_all.csv')
scipy_data = np.genfromtxt(datafile, delimiter=',', dtype=None, names=True)

print cupsoda_data.dtype.names
print scipy_data.dtype.names

for model in ['tyson', 'ras', 'earm']:

    plt.figure(model)
    
    # SciPy
    xdata = [d['nsims'] for d in scipy_data if d['model'] == model]
    ydata = [d['scipytime'] for d in scipy_data if d['model'] == model]
    plt.plot(xdata, ydata, 'b', lw=3, label='SciPy (lsoda)')

    # cupSODA
    xdata = []
    for x in [d['nsims'] for d in cupsoda_data if d['model'] == model]:
        if x not in xdata:
            xdata.append(x)
    
    fmt = ['x','o','-']
    lmem = ['global', 'shared', 'hybrid']
    colors = ['red', 'green']
    labels=['PySB', 'cupSODA']
    for i,type in enumerate(['pythontime', 'cupsodatime']):
        for mem in [0,1,2]:
            ydata = []
            for x in xdata:
                ydata.append(min([d[type] for d in cupsoda_data 
                              if  d['model'] == model
                              and d['mem'] == mem
                              and d['nsims'] == x]))
            plt.plot(xdata, ydata, fmt[mem], ms=12, lw=3, mew=2, mfc='none', mec=colors[i], color=colors[i], label='%s-%s' % (labels[i],lmem[mem]))
    
    plt.legend(loc=0)
    plt.xlabel('# of simulations')
    plt.ylabel('time (s)')
    plt.xscale('log')
    plt.yscale('log')
    
    plt.savefig(os.path.join(figs,'%s_timings.pdf' % model))

plt.show()

