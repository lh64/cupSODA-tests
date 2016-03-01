import numpy as np
# import matplotlib.pyplot as plt
import pylab as plt
import mpl_toolkits.axes_grid1 as axgrid
import os
from earm.lopez_embedded import model

figs = os.path.join('..', 'FIGS')
if not os.path.exists(figs):
    os.makedirs(figs)

proteins_of_interest = []
for i in model.initial_conditions:
    proteins_of_interest.append(i[1].name)

colors = 'RdGy'
colors = 'PiYG'
colors = 'seismic'




vals = np.hstack((np.linspace(.7, .9, 5), np.logspace(0, .3, 5)))
image1 = np.loadtxt('../parameters_486_gpu_new_image_matrix.csv')
image2 = np.loadtxt('../parameters_911_gpu_new_image_matrix.csv')
all_runs_1 = []
all_runs_2 = []
for i in range(0, len(image1), len(vals)):
    tmp = image1[:, i:i + len(vals)].flatten()
    tmp = tmp[tmp != 0]
    all_runs_1.append(tmp)
    tmp = image2[:, i:i + len(vals)].flatten()
    tmp = tmp[tmp != 0]
    all_runs_2.append(tmp)

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


# vmax = max(np.abs(image.min()), image.max())
# vmin = -1 * vmax

ax2 = fig.add_subplot(1, 2, 2)
im = ax2.imshow(image2, interpolation='nearest', origin='lower', cmap=plt.get_cmap(colors),
                vmin=vmin, vmax=vmax,extent=[0, n, 0, n])

plt.xticks(shape_label, proteins_of_interest, rotation='vertical', fontsize=14)
plt.yticks(shape_label, proteins_of_interest, fontsize=14)
xticks = np.arange(0, n)
xticks = ([i for i in range(0, len(image2), len(vals))])
ax2.set_xticks(xticks, minor=True)
ax2.set_yticks(xticks, minor=True)
ax2.yaxis.tick_right()
plt.grid(True, which='minor', axis='both', linestyle='--')

divider = axgrid.make_axes_locatable(ax2)
cax = divider.append_axes("top", size="5%", pad=0.3)

plt.colorbar(im, cax=cax, orientation='horizontal')
ax1.annotate('A', xy=(0, 1), xycoords='axes fraction', fontsize=20,
             xytext=(0, 75), textcoords='offset points',
             ha='left', va='top')
ax2.annotate('B', xy=(0, 1), xycoords='axes fraction', fontsize=20,
             xytext=(0, 75), textcoords='offset points',
             ha='left', va='top')

plt.tight_layout()
plt.savefig(os.path.join(figs, 'earm_heatplot1.eps'), bbox_tight='True')
plt.savefig(os.path.join(figs, 'earm_heatplot1.png'), bbox_tight='True')

# plt.close()

# ax4 = plt.subplot2grid((3,2),(0,1),rowspan=3)

plt.figure(figsize=(12,6))

ax4 = plt.subplot(1, 2, 1)
ax5 = plt.subplot(1, 2, 2)

ax4.boxplot(all_runs_1, vert=False, labels=None, showfliers=False)
ax4.set_xlabel('Change of time of death (%)', fontsize=14)
xtickNames = plt.setp(ax4, yticklabels=proteins_of_interest)
ax4.yaxis.tick_left()
# ax4.annotate('D', xy=(0, 1), xycoords='axes fraction', fontsize=20,
#              xytext=(-15, 18), textcoords='offset points',
#              ha='left', va='top')


ax5.boxplot(all_runs_2, vert=False, labels=None, showfliers=False)
ax5.set_xlabel('Percent change in time-to-death (%)', fontsize=14)
xtickNames = plt.setp(ax5, yticklabels=proteins_of_interest)
ax5.yaxis.tick_right()

ax4.annotate('C', xy=(0, 1), xycoords='axes fraction', fontsize=20,
             xytext=(0, 25), textcoords='offset points',
             ha='left', va='top')
ax5.annotate('D', xy=(0, 1), xycoords='axes fraction', fontsize=20,
             xytext=(0, 25), textcoords='offset points',
             ha='left', va='top')
plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.savefig(os.path.join(figs, 'earm_boxplot.eps'), bbox_tight='True')
plt.savefig(os.path.join(figs, 'earm_boxplot.png'), bbox_tight='True',dpi=300)
#plt.savefig(os.path.join(figs, 'both_earm_boxplot.eps'), bbox_tight='True')
#plt.savefig(os.path.join(figs, 'both_earm_boxplot.png'), bbox_tight='True')
plt.show()
plt.close()

image1 = image1[:10,-10:]
n = len(image1)
plt.figure()
ax7 = plt.subplot(1, 1, 1)
im = ax7.imshow(image1, interpolation='nearest',
                origin='lower', cmap=plt.get_cmap(colors),
                vmin=vmin, vmax=vmax,
                extent=[0, n, 0, n]
                )
vals = np.logspace(-.1, .1, 10)
vals = vals.round(2)
spacing = np.linspace(.5,9.5,10)
print spacing
plt.xticks(spacing, vals, rotation='vertical', fontsize=14)
plt.yticks(spacing, vals, fontsize=14)
plt.xlabel('L_0',fontsize=18)
plt.ylabel('Smac_0',fontsize=18)
ax7.xaxis.labelpad = 30
ax7.yaxis.labelpad = 30
divider = axgrid.make_axes_locatable(ax7)
cax = divider.append_axes("right", size="5%", pad=0.3)
plt.colorbar(im, cax=cax)

plt.tight_layout()

plt.savefig('right_corner.eps', bbox_tight='True')
plt.savefig('right_corner.png', bbox_tight='True')
