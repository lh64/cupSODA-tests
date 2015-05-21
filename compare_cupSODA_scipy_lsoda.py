# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 22:00:07 2015

@author: James C. Pino
"""
import time
import pysb.integrate
import pysb
import numpy as np
import pylab as plt
import os
from pysb.tools.cupsoda import *
from pysb.bng import generate_equations
import sys

name = sys.argv[1]
#name = 'tyson'
if name == 'ras':
    from ras_amp_pka import model
    tspan = np.linspace(0,1500,100)
elif name == 'earm':
    from earm.lopez_embedded import model
    tspan = np.linspace(0, 20000,100)
elif name == 'tyson':
    from pysb.examples.tyson_oscillator import model
    tspan = np.linspace(0,100,100)

run = 'scipy-mp'
#run = 'scipy'
#run ="cupSODA"

generate_equations(model)

ATOL = 1e-6
RTOL = 1e-6
mxstep=20000
det = 1
vol = 0
card = 'K20c'
#card = 'gtx760'



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


global output
output = "model,nsims,tpb,mem,pythontime,rtol,atol,mxsteps,t_end,n_steps,deterministic,vol,card\n"
def main(number_particles):
    global output
    num_particles = int(number_particles)
    if run == "cupSODA":
        c_matrix = np.zeros((num_particles, len(model.reactions)))
        rate_args = []
        for rxn in model.reactions:
            rate_args.append([arg for arg in rxn['rate'].args if not re.match("_*s",str(arg))])
        for j in range(len(model.reactions)):
            rate = 1.0
            for r in rate_args[j]:
                x = str(r)
                if x in par_dict.keys():
                    rate *= par_vals[par_dict[x]]
                else:
                    rate *= float(x)
            c_matrix[:,j] = rate

        MX_0 = np.zeros((num_particles,len(model.species)))
        for i in xrange(len(model.initial_conditions)):
            for j in xrange(len(model.species)):
                if str(model.initial_conditions[i][0]) == str(model.species[j]): # The ComplexPattern objects are not the same, even though they refer to the same species (ask about this)
                    x = model.initial_conditions[i][1]
                    MX_0[:,j] = [x.value for each in xrange(num_particles)]
                    break
        for mem in [0,1,2]:
            for i in  (1,2,4,8,16,32,64):
                if num_particles < i :
                    continue
                if i == 64 and name == 'earm':
                    continue
                solver = CupSODASolver(model, tspan, atol=ATOL, rtol=RTOL, verbose=True)
                Start = time.time()
                solver.run(c_matrix, MX_0 , n_blocks = np.int(num_particles/i), \
                outdir=os.path.join('.','CUPSODA_%s')%model.name, gpu=2,max_steps=mxstep,load_conc_data=False,memory_usage=mem)
                output +='%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n'% \
                (name,num_particles,str(i),mem,time.time()-Start,RTOL,ATOL,len(tspan),mxstep,np.max(tspan),det,vol,card)
                print name,num_particles,str(i),mem,time.time()-Start,RTOL,ATOL,len(tspan),np.max(tspan),det,vol,card
                print 'out==',solver.yobs[0][0],solver.yobs[0][-1],'==out'
                os.system('rm -r %s'%os.path.join('.','CUPSODA_%s')%model.name)
                print 'removed directory'

    if run =="scipy":
        solver = pysb.integrate.Solver(model, tspan,rtol=RTOL, atol=ATOL, integrator='lsoda', nsteps=mxstep)
        Start = time.time()
        for i in xrange(number_particles):
            Time[i] = solver.run()
        print 'sim = %s , time = %s sec' % (number_particles,time.time() - Start)
    
    if run == 'scipy-mp':
        import multiprocessing 
        global solver
        solver = pysb.integrate.Solver(model, tspan,rtol=RTOL, atol=ATOL, integrator='lsoda', nsteps=mxstep)
        c_matrix = np.zeros((num_particles, len(nominal_values)))
        c_matrix[:,:] = nominal_values
        pool = multiprocessing.Pool()
        Start = time.time()
        pool.map(RUN,c_matrix)
        print 'sim = %s , time = %s sec' % (number_particles,time.time() - Start)
def RUN(params):
    solver.run(params)

if run =='scipy':
    for j in (10,100,1000,10000,100000):
        main(j)
if run =='scipy-mp':
    for j in (10,100,1000,10000,100000):
        main(j)
if run == 'cupSODA':
    set_cupSODA_path("/home/pinojc/CUPSODA")
    for j in (10,100,1000,10000,100000):
        main(j)
    outFile = open(sys.argv[2],'w')
    outFile.write(output)
    outFile.close()
    quit()




