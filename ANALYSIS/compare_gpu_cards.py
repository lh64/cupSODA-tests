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
indies = [(0, 0), (1, 0), (2, 0)]
plt.figure()
for count, model in enumerate(['tyson', 'ras', 'earm']):
    if model == 'tyson':
        ax1 = plt.subplot2grid((3, 1), indies[count])
    if model == 'ras':
        ax2 = plt.subplot2grid((3, 1), indies[count], sharex=ax1)
    if model == 'earm':
        ax3 = plt.subplot2grid((3, 1), indies[count], sharex=ax1)
    xdata = []
    for x in [d['nsims'] for d in cupsoda_data if d['model'] == model]:
        if x not in xdata:
            xdata.append(x)

    fmt = ['^-', 's-', '*-']
    cards = ['gtx980-mule','gtx760-lolab','gtx980-diablo']
    colors = ['c', 'magenta', 'green', ]
    labels=[ 'gtx970','gtx760','gtx980-TI']

    for i,type in enumerate(cards):
        ydata = []
        Xdata = []
        for x in xdata:
            print x
            data = cupsoda_data[cupsoda_data['model']==model]
            data = data[data['card'] == cards[i]]
            data = data[data['nsims'] == x]
            print data
            if len(data) == 0:
                #xdata.remove(x)
                continue
            else:
                ydata.append(float(data['cupsodatime']))
                Xdata.append(x)

        if model == 'tyson':
            ax1.plot(Xdata, ydata,fmt[i], ms=10, lw=3, mew=2, mec=colors[i], color=colors[i],
                     label=labels[i])
        if model == 'ras':
            ax2.plot(Xdata, ydata,fmt[i], ms=10, lw=3, mew=2, mec=colors[i], color=colors[i],
                     label=labels[i])
        if model == 'earm':
            ax3.plot(Xdata, ydata,fmt[i], ms=10, lw=3, mew=2, mec=colors[i], color=colors[i],
                     label=labels[i])
        plt.xscale('log')
        plt.yscale('log')


plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)

ax1.yaxis.set_tick_params(labelsize=14)
ax2.yaxis.set_tick_params(labelsize=14)
ax3.yaxis.set_tick_params(labelsize=14)



ax1.legend(fontsize=14, bbox_to_anchor=(.75, 1.0), fancybox=True)
ax2.set_ylabel('time (s)', fontsize=14)
ax3.set_xlabel("Number of simulations", fontsize=14)

ax1.annotate('A', xy=(0, 1), xycoords='axes fraction', fontsize=20,
             xytext=(-60, 10), textcoords='offset points',
             ha='left', va='top')
ax2.annotate('B', xy=(0, 1), xycoords='axes fraction', fontsize=20,
             xytext=(-60, 10), textcoords='offset points',
             ha='left', va='top')
ax3.annotate('C', xy=(0, 1), xycoords='axes fraction', fontsize=20,
             xytext=(-60, 10), textcoords='offset points',
             ha='left', va='top')

ax1.annotate('Tyson', xy=(0, 1), xycoords='axes fraction', fontsize=20,
             xytext=(5, -5), textcoords='offset points',
             ha='left', va='top')
ax2.annotate('Ras/cAMP/PKA', xy=(0, 1), xycoords='axes fraction', fontsize=20,
             xytext=(5, -5), textcoords='offset points',
             ha='left', va='top')
ax3.annotate('EARM', xy=(0, 1), xycoords='axes fraction', fontsize=20,
             xytext=(5, -5), textcoords='offset points',
             ha='left', va='top')



y_lim = [.1,500]
x_lim = [-10,10000]
ax1.set_xlim(x_lim)
ax1.set_ylim(y_lim)
ax2.set_xlim(x_lim)
ax2.set_ylim([1,500])
ax3.set_xlim(x_lim)
ax3.set_ylim([1,500])


plt.tight_layout()
plt.subplots_adjust(hspace=0.0)

plt.savefig(os.path.join(figs, 'compare_gpu.eps' ), bbox_tight='True')
plt.savefig(os.path.join(figs, 'compare_gpu.png' ), bbox_tight='True')
plt.show()

