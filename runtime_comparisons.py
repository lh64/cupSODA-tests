import numpy as np
import pylab as plt
import os
import os.path

if not os.path.exists('FIGS'):
    os.makedirs('FIGS')

cupsoda_data = np.genfromtxt('cupsoda_timings_all.csv', delimiter=',', dtype=None, names=True)
scipy_data = np.genfromtxt('scipy_timings_all.csv', delimiter=',', dtype=None, names=True)

print cupsoda_data.dtype.names
print scipy_data.dtype.names

models = ['tyson', 'ras', 'earm']

for model in models:

    plt.figure(model)
    
    # SciPy
    xdata = [d['nsims'] for d in scipy_data if d['model'] == model]
    ydata = [d['scipytime'] for d in scipy_data if d['model'] == model]
    plt.plot(xdata, ydata, 'b', lw=3, label='SciPy')

    # cupSODA
    xdata = []
    for x in [d['nsims'] for d in cupsoda_data if d['model'] == model and d['mem'] == 0]:
        if x not in xdata:
            xdata.append(x)
    
    fmt = ['x','o','-']
    lmem = ['Global', 'Shared', 'Hybrid']
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
    
    plt.savefig(os.path.join('FIGS','%s_timings.pdf' % model))

plt.show()

