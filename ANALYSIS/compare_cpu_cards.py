import os
import numpy as np
import pylab as plt

figs = os.path.join('..', 'FIGS')
if not os.path.exists(figs):
    os.makedirs(figs)

datafile = os.path.join('..', 'scipy_times_new.csv')
cupsoda_data = np.genfromtxt(datafile, delimiter=',', dtype=None, names=True)

fig = plt.figure(figsize=(8, 6))
numbers = [1, 3, 5]
indies = [(0, 0), (1, 0), (2, 0)]
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
    xdata = np.unique(xdata)
    fmt = ['x', 'o', '-']
    lmem = ['global', 'shared', 'hybrid']
    cards = ['diablo', 'lolab-760','gtx-970-mule']
    colors = ['red', 'green', 'blue']
    labels = ['E5-2667 v3 @ 3.20GHz', 'i7-4820K CPU @ 3.70GHz', 'i7-5930K CPU @ 3.50GHz']

    for i, type in enumerate(cards):
        ydata = []
        Xdata = []
        for x in xdata:
            data = cupsoda_data[cupsoda_data['model'] == model]
            data = data[data['cpu'] == cards[i]]
            data = data[data['nsims'] == x]
            if len(data) == 0:
                continue
            ydata.append(float(data['scipytime']))
            Xdata.append(x)
        if model == 'tyson':
            ax1.plot(Xdata, ydata,'o-', ms=10, lw=3, mew=2, mec=colors[i], color=colors[i],
                     label=labels[i])
        if model == 'ras':
            ax2.plot(Xdata, ydata,'o-', ms=10, lw=3, mew=2, mec=colors[i], color=colors[i],
                     label=labels[i])
        if model == 'earm':
            ax3.plot(Xdata, ydata,'o-', ms=10, lw=3, mew=2, mec=colors[i], color=colors[i],
                     label=labels[i])

        plt.xscale('log')
        plt.yscale('log')

ax1.set_xlim(0, 10000)
ax1.set_ylim(0, 100)

ax2.set_xlim(0, 10000)
ax2.set_ylim(0, 1000)

ax3.set_xlim(0, 10000)
ax3.set_ylim(0, 1000)

plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)

ax1.set_yticks(ax1.get_yticks()[3:])
ax2.set_yticks(ax2.get_yticks()[3:])

ax1.yaxis.set_tick_params(labelsize=14)
ax2.yaxis.set_tick_params(labelsize=14)
ax3.yaxis.set_tick_params(labelsize=14)
ax1.legend(fontsize=14, bbox_to_anchor=(.85, 1.4), fancybox=True)
ax2.set_ylabel('time (s)', fontsize=14)
ax3.set_xlabel("Number of simulations", fontsize=14)
distance = (-60, 10)
ax1.annotate('A', xy=(0, 1), xycoords='axes fraction', fontsize=20,
             xytext=(-60, 18), textcoords='offset points',
             ha='left', va='top')
ax2.annotate('B', xy=(0, 1), xycoords='axes fraction', fontsize=20,
             xytext=distance, textcoords='offset points',
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

plt.savefig(os.path.join(figs, 'compare_cpu.eps' ), bbox_tight='True')
plt.savefig(os.path.join(figs, 'compare_cpu.png' ), bbox_tight='True')
plt.show()
