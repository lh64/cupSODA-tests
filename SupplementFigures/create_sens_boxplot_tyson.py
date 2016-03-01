import numpy as np
# import matplotlib.pyplot as plt
import pylab as plt
import mpl_toolkits.axes_grid1 as axgrid
import os
from pysb.examples.tyson_oscillator import model

figs = os.path.join('..', 'FIGS')
if not os.path.exists(figs):
    os.makedirs(figs)

proteins_of_interest = []
for i in model.initial_conditions:
    proteins_of_interest.append(i[1].name)

colors = 'RdGy'
colors = 'PiYG'
colors = 'seismic'




vals = np.logspace(-.1, .1, 20)
image1 = np.loadtxt('sens_tyson_matrix.csv')
all_runs_1 = []

for i in range(0, len(image1), len(vals)):
    tmp2 = []
    for j in range(0,len(vals)):
        print i,j
        tmp = np.sum(image1[i:i+len(vals), j])
        print tmp
        #tmp = tmp[tmp != 0]
        tmp2.append(tmp)
    all_runs_1.append(tmp2)
print all_runs_1

vmax = max(np.abs(image1.min()), image1.max())
vmin = -1 * vmax
fig = plt.figure(figsize=(12, 7))
ax1 = fig.add_subplot(1, 2, 1)


n = len(image1)
im = ax1.imshow(image1, interpolation='nearest',
                origin='lower', cmap=plt.get_cmap(colors),
                vmin=vmin, vmax=vmax,
                extent=[0, n, 0, n]
                )
shape_label = np.arange(len(vals) / 2, len(image1), len(vals))
plt.xticks(shape_label, proteins_of_interest, rotation='vertical', fontsize=14)
plt.yticks(shape_label, proteins_of_interest, fontsize=14)
xticks = ([i for i in range(0, len(image1), len(vals))])
ax1.set_xticks(xticks, minor=True)
ax1.set_yticks(xticks, minor=True)
plt.grid(True, which='minor', axis='both', linestyle='--')
divider = axgrid.make_axes_locatable(ax1)
cax = divider.append_axes("top", size="5%", pad=0.3)
plt.colorbar(im, cax=cax, orientation='horizontal')


plt.grid(True, which='minor', axis='both', linestyle='--')

ax1.annotate('A', xy=(0, 1), xycoords='axes fraction', fontsize=20,
             xytext=(0, 75), textcoords='offset points',
             ha='left', va='top')

plt.tight_layout()


ax2 = plt.subplot(1, 2, 2)


ax2.boxplot(all_runs_1, vert=False, labels=None, showfliers=False)
ax2.set_xlabel('Change in maximum accumulation of cAMP (%)', fontsize=14)
xtickNames = plt.setp(ax2, yticklabels=proteins_of_interest)
ax2.yaxis.tick_left()



ax2.annotate('B', xy=(0, 1), xycoords='axes fraction', fontsize=20,
             xytext=(0, 25), textcoords='offset points',
             ha='left', va='top')

plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig('tyson_sensitivity_boxplot.eps', bbox_tight='True')
plt.savefig('tyson_sensitivity_boxplot.png', bbox_tight='True',dpi=300)
plt.show()

