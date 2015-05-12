import pylab as plt
import numpy as np
import pysb.integrate
from pysb.integrate import odesolve
import matplotlib.cm as cm
import os

figs = os.path.join('..','FIGS')
if not os.path.exists(figs):
    os.makedirs(figs)

models = ['tyson', 'ras', 'earm']

for m in models:
    
    plt.figure(m)
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    if m == 'tyson':
        from pysb.examples.tyson_oscillator import model
        tspan = np.linspace(0, 100, 100)
    elif m == 'ras':
        from ras_amp_pka import model
        tspan = np.linspace(0,1500,100)
    elif m == 'earm':
        from earm.lopez_embedded import model
        tspan = np.linspace(0, 20000,100)
        
    colors = cm.rainbow(np.linspace(0, 1, len(model.observables)))
    
    x = odesolve(model, tspan, atol=1e-6, rtol=1e-6, nsteps=20000, verbose=True)

    for i,obs in enumerate(model.observables):
        ax.plot(tspan, x[obs.name]/np.nanmax(x[obs.name]), lw=3, label=obs.name.lstrip('obs_'), c=colors[i])
        
    if m == 'ras':
        ax.legend(loc='center left', bbox_to_anchor=(1.05,0.5), fontsize=8)
    else:
        ax.legend(loc='center left', bbox_to_anchor=(1.05,0.5))
        
    plt.xlabel('time (s)')
    plt.ylabel('normalized population')
    
    plt.savefig(os.path.join(figs,'%s_traj.pdf' % m))
    
plt.show()

