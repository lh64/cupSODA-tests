# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 22:00:07 2015

@author: James C. Pino
"""
import time   
import pysb.integrate
import numpy as np
import pylab as plt
import os
from pysb.tools.cupsoda import *
from pysb.bng import generate_equations
import sys

if sys.argv[1] == 'ras':
    from ras_amp_pka import model
    tspan = np.linspace(0,1500,1000)
if sys.argv[1] == 'earm':
    from earm.lopez_embedded import model
    tspan = np.linspace(0, 20000,1000)
if sys.argv[1] == 'tyson':
    from pysb.examples.tyson_oscillator import model
    tspan = np.linspace(0,100,100)

generate_equations(model)

#tspan = np.linspace(0, 20000,1000)
#tspan = np.linspace(0, 100,100)
#tspan = np.linspace(0,1500,1000)
set_cupSODA_path("/home/pinojc/CUPSODA")


params_names = [p.name for p in model.parameters]
init_name = [p[1].name for p in model.initial_conditions]
par_names = []
for parm in params_names:
    if parm in init_name:
        continue
    else:
        par_names.append(parm)
rate_params = model.parameters_rules()
rate_mask = np.array([p in rate_params for p in model.parameters])
nominal_values = np.array([p.value for p in model.parameters])
xnominal = np.log10(nominal_values[rate_mask])
par_dict = {par_names[i] : i for i in range(len(par_names))}
par_vals = np.array([model.parameters[nm].value for nm in par_names])


run ="cupSODA"


def main(number_particles,n_blocks):
    n_blocks = n_blocks
    num_particles = int(number_particles)
    if run == "cupSODA":
        c_matrix = np.zeros((num_particles, len(model.reactions)))
        rate_args = []
        for rxn in model.reactions:
            #print rxn['rate']
            rate_args.append([arg for arg in rxn['rate'].args if not re.match("_*s",str(arg))])
        for j in range(len(model.reactions)):
            rate = 1.0
            for r in rate_args[j]:
                x = str(r)
                if x in par_dict.keys():
                    rate *= par_vals[par_dict[x]] # model.parameters[x].value
                else:
                    rate *= float(x)
            c_matrix[:,j] = rate        
        
        #print par_vals - c_matrix
        # Initial concentrations
        MX_0 = np.zeros((num_particles,len(model.species)))
        for i in xrange(len(model.initial_conditions)):
            for j in xrange(len(model.species)):
                if str(model.initial_conditions[i][0]) == str(model.species[j]): # The ComplexPattern objects are not the same, even though they refer to the same species (ask about this)
                    x = model.initial_conditions[i][1]
                    MX_0[:,j] = [x.value for each in xrange(num_particles)]
                    break  
        solver = CupSODASolver(model, tspan, atol=1e-6, rtol=1e-6, verbose=False)
        Start = time.time()
        solver.run(c_matrix, MX_0 , n_blocks = n_blocks, \
                   outdir=os.path.join('.','CUPSODA_%s') % model.name, gpu=2, load_conc_data=False) 
        print 'sim = %s ,theads/block= %s, time = %s sec' % (num_particles,num_particles/n_blocks,time.time() - Start)
        for i in solver.y[0]:
            print i,
    if run =="scipy":
        solver = pysb.integrate.Solver(model, tspan, integrator='lsoda', \
                                       rtol=1e-6, atol=1e-6, nsteps=10000)
        Start = time.time()
        for i in xrange(num_particles):
            solver.run()
        print 'sim = %s , time = %s sec' % (num_particles,time.time() - Start)
if run =='scipy':                
    #main(1e1,1)
    #main(1e2,1)
    main(1e3,1)
    main(1e4,1)
    main(1e5,1)
if run == 'cupSODA':
    #for i in (32,64,128,256,512,1024,None):
    for i in (8,10,20,25,32,40,50):
        #main(1,i)
        #main(1e1,i)
        #main(1e2,i)
        #main(1e3,np.round(1e3/i))
        main(1e4,np.round(1e4/i))
        main(1e5,np.round(1e5//i))
        main(1e6,np.round(1e6/i))    



