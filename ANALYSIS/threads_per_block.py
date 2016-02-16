import numpy as np
import pylab as plt
import os
from matplotlib.ticker import ScalarFormatter 

figs = os.path.join('..','FIGS')
if not os.path.exists(figs):
    os.makedirs(figs)

datafile = os.path.join('..', 'all-GPU_timing.csv')
cupsoda_data = np.genfromtxt(datafile, delimiter=',', dtype=None, names=True)

print cupsoda_data.dtype.names

colors = ['red', 'blue', 'green', 'magenta']
fmt = ['x','o','-']

for model in ['tyson', 'ras', 'earm']:
    
    plt.figure(model)
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
    plt.annotate(r"Global: x", (1.05,0.6), xycoords='axes fraction', fontsize=18)
    plt.annotate(r"Shared: o", (1.05,0.53), xycoords='axes fraction', fontsize=18)
    plt.annotate(r"Hybrid: --", (1.05,0.46), xycoords='axes fraction', fontsize=18)
    
    nsims = []
    for x in [d['nsims'] for d in cupsoda_data if d['model'] == model]:
        if x not in nsims and x > 10: # exclude nsims = 10
            nsims.append(x)
    for i,n in enumerate(nsims):
        for mem in [0,1,2]:
            xdata = [d['tpb'] for d in cupsoda_data
                     if  d['model'] == model
                     and d['nsims'] == n
                     and d['mem'] == mem]
            ydata = [d['cupsodatime'] for d in cupsoda_data
                     if  d['model'] == model
                     and d['nsims'] == n
                     and d['mem'] == mem]
            if mem < 2:
                ax.plot(xdata, ydata, fmt[mem], ms=12, lw=3, mew=2, mfc='none', mec=colors[i], color=colors[i])
            else:
                ax.plot(xdata, ydata, fmt[mem], ms=12, lw=3, mew=2, mfc='none', mec=colors[i], color=colors[i], label = 'nsims: %.0e' % float(n))
    ax.legend(loc='upper left', bbox_to_anchor=(1.,1.))
    ax.set_xscale('log', basex=2)
    ax.set_yscale('log')
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_ticks(xdata)
    ax.set_xlabel('threads/block')
    ax.set_ylabel('time (s)')
    
    plt.savefig(os.path.join(figs,'%s_tpb.pdf' % model))

plt.show()

     