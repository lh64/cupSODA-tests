import os
import numpy as np
import pylab as plt

figs = os.path.join('..', 'FIGS')
if not os.path.exists(figs):
    os.makedirs(figs)

datafile = os.path.join('..', 'scipy_times_new.csv')
cpu_data = np.genfromtxt(datafile, delimiter=',', dtype=None, names=True)

datafile = os.path.join('..', 'all-GPU_timing.csv')
cupsoda_data = np.genfromtxt(datafile, delimiter=',', dtype=None, names=True)
fig = plt.figure(figsize=(12, 6))
numbers = [1, 3, 5]
indies = [(0, 0), (1, 0), (2, 0)]
for count, model in enumerate(['tyson', 'ras', 'earm']):

    if model == 'tyson':
        ax1 = plt.subplot2grid((3, 2), indies[count])
    if model == 'ras':
        ax2 = plt.subplot2grid((3, 2), indies[count], sharex=ax1)
    if model == 'earm':
        ax3 = plt.subplot2grid((3, 2), indies[count], sharex=ax1)

    xdata = []
    for x in [d['nsims'] for d in cpu_data if d['model'] == model]:
        if x not in xdata:
            xdata.append(x)
    xdata = np.unique(xdata)
    fmt = ['x', 'o', '-']
    lmem = ['global', 'shared', 'hybrid']
    cards = ['diablo', 'lolab-760','gtx-970-mule']
    colors = ['red', 'green', 'blue']
    labels = ['E5-2667 v3', 'i7-4820K', 'i7-5930K']

    for i, type in enumerate(cards):
        ydata = []
        Xdata = []
        for x in xdata:
            data = cpu_data[cpu_data['model'] == model]
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


indies = [(0, 1), (1, 1), (2, 1)]
for count, model in enumerate(['tyson', 'ras', 'earm']):
    if model == 'tyson':
        ax4 = plt.subplot2grid((3, 2), indies[count],sharey=ax1)
    if model == 'ras':
        ax5 = plt.subplot2grid((3, 2), indies[count], sharex=ax4,sharey=ax2)
    if model == 'earm':
        ax6 = plt.subplot2grid((3, 2), indies[count], sharex=ax5,sharey=ax3)
    xdata = []
    for x in [d['nsims'] for d in cupsoda_data if d['model'] == model]:
        if x not in xdata:
            xdata.append(x)

    fmt = ['x','o','-']
    cards = ['gtx980-mule','gtx980-diablo','gtx760-lolab']
    colors = ['red', 'green','blue']
    labels=['gtx970', 'gtx980','gtx760']

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
                xdata.remove(x)
                continue
            else:
                ydata.append(float(data['pythontime']))
                Xdata.append(x)

        if model == 'tyson':
            ax4.plot(Xdata, ydata,'o-', ms=10, lw=3, mew=2, mec=colors[i], color=colors[i],
                     label=labels[i])
        if model == 'ras':
            ax5.plot(Xdata, ydata,'o-', ms=10, lw=3, mew=2, mec=colors[i], color=colors[i],
                     label=labels[i])
        if model == 'earm':
            ax6.plot(Xdata, ydata,'o-', ms=10, lw=3, mew=2, mec=colors[i], color=colors[i],
                     label=labels[i])
        plt.xscale('log')
        plt.yscale('log')


plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
ax1.set_xlim(-10, 10000)
ax1.set_ylim(-10, 1000)
ax2.set_xlim(-10, 10000)
ax2.set_ylim(-10, 1000)
ax3.set_xlim(-10, 10000)
ax3.set_ylim(-10, 1000)

#ax1.set_yticks(ax1.get_yticks()[1:])
#ax2.set_yticks(ax2.get_yticks()[1:])
#ax3.set_yticks(ax3.get_yticks()[1:])

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



plt.setp(ax4.get_xticklabels(), visible=False)
plt.setp(ax5.get_xticklabels(), visible=False)

plt.setp(ax4.get_yticklabels(), visible=False)
plt.setp(ax5.get_yticklabels(), visible=False)
plt.setp(ax6.get_yticklabels(), visible=False)

ax4.yaxis.set_tick_params(labelsize=14)
ax5.yaxis.set_tick_params(labelsize=14)
ax6.yaxis.set_tick_params(labelsize=14)


ax4.legend(fontsize=14, bbox_to_anchor=(.85, 1.4), fancybox=True)
#ax5.set_ylabel('time (s)', fontsize=14)
ax6.set_xlabel("Number of simulations", fontsize=14)
y_lim = [-10,1000]
x_lim = [-10,10000]
ax1.set_xlim(x_lim)
ax1.set_ylim(y_lim)
ax2.set_xlim(x_lim)
ax2.set_ylim(y_lim)
ax3.set_xlim(x_lim)
ax3.set_ylim(y_lim)
ax4.set_xlim(x_lim)
# ax4.set_ylim(y_lim)
ax5.set_xlim(x_lim)
# ax5.set_ylim(y_lim)
ax6.set_xlim(x_lim)
# ax6.set_ylim(y_lim)

ax4.annotate('D', xy=(0, 1), xycoords='axes fraction', fontsize=20,
             xytext=(-20, 18), textcoords='offset points',
             ha='left', va='top')
ax5.annotate('E', xy=(0, 1), xycoords='axes fraction', fontsize=20,
             xytext=(-20, 10), textcoords='offset points',
             ha='left', va='top')
ax6.annotate('F', xy=(0, 1), xycoords='axes fraction', fontsize=20,
             xytext=(-20, 10), textcoords='offset points',
             ha='left', va='top')

plt.subplots_adjust(wspace=0.1)
plt.savefig(os.path.join(figs, 'compare_both_cpu.eps' ), bbox_tight='True')
plt.savefig(os.path.join(figs, 'compare_both_cpu.png' ), bbox_tight='True')
plt.show()

