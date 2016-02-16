import numpy as np
import pylab as plt
import os

figs = os.path.join('..','FIGS')
if not os.path.exists(figs):
    os.makedirs(figs)

datafile = os.path.join('..', 'all-GPU_timing.csv')
cupsoda_data = np.genfromtxt(datafile, delimiter=',', dtype=None, names=True)

datafile = os.path.join('..', 'scipy_timings_all.csv')
scipy_data = np.genfromtxt(datafile, delimiter=',', dtype=None, names=True)

#print cupsoda_data.dtype.names
#print cupsoda_data['card']
#print scipy_data.dtype.names


for model in ['tyson', 'ras', 'earm']:
    #if model == 'ras':
    #    break
    plt.figure(model)

    # SciPy
    xdata = [d['nsims'] for d in scipy_data if d['model'] == model and d['num_cpu'] == 1]
    ydata = [d['scipytime'] for d in scipy_data if d['model'] == model and d['num_cpu'] == 1]

    plt.plot(xdata, ydata, label='SciPy (lsoda)')
    # cupSODA
    xdata = []
    for x in [d['nsims'] for d in cupsoda_data if d['model'] == model]:
        if x not in xdata:
            xdata.append(x)

    fmt = ['x','o','-']
    lmem = ['global', 'shared', 'hybrid']
    cards = ['gtx980-mule','gtx980-diablo']
    colors = ['red', 'green','blue']
    labels=['PySB', 'cupSODA']

    for i,type in enumerate(cards):
        ydata = []
        Xdata = []
        for x in xdata:
            print x
            data = cupsoda_data[cupsoda_data['model']==model]
            print "=model",data
            data = data[data['card'] == cards[i]]
            print "=card",data
            data = data[data['nsims'] == x]
            print data
            print data['pythontime']
            if len(data) == 0:
                xdata.remove(x)
                continue
            ydata.append(float(data['pythontime']))
            Xdata.append(x)

        plt.plot(Xdata, ydata, '-o', ms=12, lw=3, mew=2, mfc='none', mec=colors[i], color=colors[i], label='%s' % (cards[i]))
    plt.legend(loc=0)
    plt.xlabel('# of simulations')
    plt.ylabel('time (s)')
    plt.xscale('log')
    plt.yscale('log')

    plt.savefig(os.path.join(figs,'%s_compare_gpu.pdf' % model))

plt.show()

