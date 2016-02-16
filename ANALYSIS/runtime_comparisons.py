import numpy as np
#import matplotlib
#matplotlib.use('AGG')
import matplotlib.pyplot as plt
#import pylab as plt
import os
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

figs = os.path.join('..','FIGS')
if not os.path.exists(figs):
    os.makedirs(figs)

datafile = os.path.join('..', '2016_all3.csv')
cupsoda_data = np.genfromtxt(datafile, delimiter=',', dtype=None, names=True)

datafile = os.path.join('..', 'scipy_timings_all.csv')
scipy_data = np.genfromtxt(datafile, delimiter=',', dtype=None, names=True)

fig = plt.figure(figsize=(12,6))
numbers = [1,3,5]
indies= [(0,0),(1,0),(2,0)]
for count, model in enumerate(['tyson', 'ras', 'earm']):

    #plt.figure(model)
    sharey=False
    if model == 'tyson':
        ax1 = plt.subplot2grid((3,2),indies[count])
    if model == 'ras':
        ax2 = plt.subplot2grid((3,2),indies[count],sharex=ax1)
    if model == 'earm':
        ax3 = plt.subplot2grid((3,2),indies[count],sharex=ax1)
    # SciPy
    xdata = [d['nsims'] for d in scipy_data if d['model'] == model and d['num_cpu'] == 1]
    ydata = [d['scipytime'] for d in scipy_data if d['model'] == model and d['num_cpu'] == 1]
    if model == 'tyson':
        ax1.plot(xdata, ydata, 'b-o',  label='SciPy (lsoda)',ms=12, lw=3, mew=0,)# mfc='none',)
    if model == 'ras':
        ax2.plot(xdata, ydata, 'b-o',  label='SciPy (lsoda)',ms=12, lw=3, mew=0,)# mfc='none',)
    if model == 'earm':
        ax3.plot(xdata, ydata, 'b-o',  label='SciPy (lsoda)',ms=12, lw=3, mew=0, )#mfc='none',)
    # cupSODA
    xdata = []
    for x in [d['nsims'] for d in cupsoda_data if d['model'] == model]:
        if x not in xdata:
            xdata.append(x)
    
    fmt = ['x','o','-v']
    lmem = ['global', 'shared', 'hybrid']
    colors = ['red', 'green']
    labels=['PySB/cupSODA', 'cupSODA']
    for i,type in enumerate(['pythontime', 'cupsodatime']):
        if type == 'cupsodatime':
            continue
        for mem in [2]:
            ydata = []
            for x in xdata:
                ydata.append([d[type] for d in cupsoda_data
                                  if  d['model'] == model
                                  and d['mem'] == mem
                                  and d['nsims'] == x
                                  ])
            if model == 'tyson':
                ax1.plot(xdata, ydata, fmt[mem] ,ms=12, lw=3, mew=2,  mec=colors[i], color=colors[i], label=labels[i])
            if model == 'ras':
                ax2.plot(xdata, ydata, fmt[mem],ms=12, lw=3, mew=2,  mec=colors[i], color=colors[i], label=labels[i])
            if model == 'earm':
                ax3.plot(xdata, ydata, fmt[mem], ms=12, lw=3, mew=2,  mec=colors[i], color=colors[i], label=labels[i])
            #plt.plot(xdata, ydata, fmt[mem], ms=12, lw=3, mew=2, mfc='none', mec=colors[i], color=colors[i], label=labels[i])
    plt.xscale('log')
    plt.yscale('log')
from matplotlib.ticker import MaxNLocator

ax1.set_xlim(0,10000)
ax1.set_ylim(0,100)

ax2.set_xlim(0,10000)
ax2.set_ylim(0,1000)

ax3.set_xlim(0,10000)
ax3.set_ylim(0,1000)


plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)

ax1.set_yticks(ax1.get_yticks()[3:])
ax2.set_yticks(ax2.get_yticks()[3:])

ax1.yaxis.set_tick_params(labelsize=14)
ax2.yaxis.set_tick_params(labelsize=14)
ax3.yaxis.set_tick_params(labelsize=14)
ax1.legend( fontsize=14, bbox_to_anchor=(.5, 1.15), fancybox=True)
ax2.set_ylabel('time (s)',fontsize=14)
ax3.set_xlabel("Number of simulations",fontsize=14)
distance = (-60,10)
ax1.annotate('A', xy=(0, 1), xycoords='axes fraction', fontsize=20,
                xytext=(-60,18), textcoords='offset points',
                ha='left', va='top')
ax2.annotate('B', xy=(0, 1), xycoords='axes fraction', fontsize=20,
                xytext=distance, textcoords='offset points',
                ha='left', va='top')
ax3.annotate('C', xy=(0, 1), xycoords='axes fraction', fontsize=20,
                xytext=(-60,10), textcoords='offset points',
                ha='left', va='top')


from earm.lopez_embedded import model
proteins_of_interest = []
for i in model.initial_conditions:
   proteins_of_interest.append(i[1].name)
vals = np.hstack((np.linspace(.7, .9, 5), np.logspace(0, .3, 5)))
image = np.loadtxt('../parameters_486_image_matrix.csv')
vmax = max(np.abs(image.min()), image.max())
vmin = -1 * vmax


#ax4 = plt.subplot(1,2,2,sharey =ax1)
ax4 = plt.subplot2grid((3,2),(0,1),rowspan=3)
all_runs = []
for i in range(0, len(image), len(vals)):
    tmp = image[:, i:i + len(vals)].flatten()
    tmp = tmp[tmp != 0]
    all_runs.append(tmp)

ax4.boxplot(all_runs, vert=False, labels=None, showfliers=False)
ax4.set_xlabel('Change of time of death (%)',fontsize=14)
xtickNames = plt.setp(ax4, yticklabels=proteins_of_interest)
ax4.yaxis.tick_right()
#plt.setp(xtickNames, rotation=45, fontsize=8)

ax4.annotate('D', xy=(0, 1), xycoords='axes fraction', fontsize=20,
                xytext=(-15, 18), textcoords='offset points',
                ha='left', va='top')
plt.tight_layout()
fig.subplots_adjust(hspace=.02,wspace=.1)

plt.savefig(os.path.join(figs,'all_timing_and_boxplot.png'),bbox_tight='True')

plt.show()