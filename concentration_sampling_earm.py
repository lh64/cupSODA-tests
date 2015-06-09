# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 17:27:41 2014

@author: pinojc
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
import scipy.interpolate
import multiprocessing

#name = sys.argv[1]
name = 'earm'
if name == 'ras':
    from ras_amp_pka import model
    tspan = np.linspace(0,1500,100)
    simulations = [10,100,1000,10000,100000]
elif name == 'tyson':
    from pysb.examples.tyson_oscillator import model
    tspan = np.linspace(0,25,50)
    simulations = [10,100,1000,10000,100000]
elif name == 'earm':
    from earm.lopez_embedded import model
    tspan = np.linspace(0, 20000,100)
    simulations = [10,100,1000,10000]
from pysb.integrate import odesolve

ref = odesolve(model, tspan, integrator='lsoda', verbose=False)

run = 'cupSODA'
multi = False
ATOL = 1e-6
RTOL = 1e-6
mxstep=20000
det = 1
vol = 0
card = 'K20c'
#puma
CPU = 'Intel-Xeon-E5-2687W-v2'
GHZ = '3.40'


generate_equations(model)
solver = pysb.integrate.Solver(model, tspan, integrator='vode', rtol=1e-6, atol=1e-6,)
proteins_of_interest = []
for i in model.initial_conditions:
    proteins_of_interest.append(i[1].name)
key_list = []
initial_tmp = dict()
for species, keys in  model.initial_conditions:
    for j in proteins_of_interest:
        if  keys.name == j:
            initial_tmp[keys.name] = keys.value

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
            
# momp_obs_total = model.parameters['Smac_0'].value
# momp_data = np.array([9810.0, 180.0, momp_obs_total])
# momp_var = np.array([7245000.0, 3600.0, 1e4])
            
def likelihood1(initial_tmp):
    
    solver.run(initial_changes=initial_tmp)
    ysim_momp = solver.yobs['aSmac']
    if np.nanmax(ysim_momp) == 0:
        print 'No aSmac!'
        ysim_momp_norm = ysim_momp
        t10 = 0
        t90 = 0
    else:
        ysim_momp_norm = ysim_momp / np.nanmax(ysim_momp)
        st, sc, sk = scipy.interpolate.splrep(tspan, ysim_momp_norm)
        try:
            t10 = scipy.interpolate.sproot((st, sc-0.10, sk))[0]
            t90 = scipy.interpolate.sproot((st, sc-0.90, sk))[0]
        except IndexError:
            t10 = 0
            t90 = 0
    td = (t10 + t90) / 2
    return   td 
#tod = likelihood1(initial_tmp)
#print tod

def likelihood(ysim_momp):

    if np.nanmax(ysim_momp) == 0:
        #print 'No aSmac!'
        ysim_momp_norm = ysim_momp
        t10 = 0
        t90 = 0
        return -1
    else:
        ysim_momp_norm = ysim_momp / np.nanmax(ysim_momp)
        st, sc, sk = scipy.interpolate.splrep(tspan, ysim_momp_norm)
        try:
            t10 = scipy.interpolate.sproot((st, sc-0.10, sk))[0]
            t90 = scipy.interpolate.sproot((st, sc-0.90, sk))[0]
        except IndexError:
            t10 = 0
            t90 = 0
    td = (t10 + t90) / 2
    return  (td-tod)/tod 
def like(data):
    return np.sum((data - ref['YT'][-1]))

vals = [.75,1.,1.25]
#vals = np.hstack((np.linspace(.7,.9,5),np.logspace(0,1,5)))

#size_of_matrix = (len(proteins_of_interest)*len(proteins_of_interest)-len(proteins_of_interest))*len(vals)*len(vals)/2
size_of_matrix = (len(proteins_of_interest)*len(proteins_of_interest))*len(vals)*len(vals) 
for each in model.initial_conditions:
    print each
for each in model.species:
    print each

c_matrix = np.zeros((size_of_matrix, len(model.reactions)))
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

MX_0 = np.zeros((size_of_matrix,len(model.species)))
for i in xrange(len(model.initial_conditions)):
    for j in xrange(len(model.species)):
        if str(model.initial_conditions[i][0]) == str(model.species[j]):
            x = model.initial_conditions[i][1].value
            MX_0[:,j] = x
initi = MX_0.copy()
image = np.zeros((len(proteins_of_interest)*len(vals),len(proteins_of_interest)*len(vals)))
done = []
counter = 0
for i,one  in enumerate(proteins_of_interest):
    for j,two in enumerate(proteins_of_interest):
#         if i == j:
#             continue
        #if one+two in done:
        #    continue
        #done.append(one+two)
        #done.append(two+one)
        for a,c in enumerate(vals):
            for b,d in enumerate(vals):
                MX_0[counter,i] *= c
                if i != j:
                    MX_0[counter,j] *= d
                counter+=1
print size_of_matrix
print counter
#plt.plot(MX_0 - initi)
#plt.show()
#quit()
if run =='scipy':
    global output
    output = "model,nsims,scipytime,rtol,atol,mxsteps,t_end,n_steps,cpu,GHz,num_cpu\n"
    solver = pysb.integrate.Solver(model, tspan,rtol=RTOL, atol=ATOL, integrator='lsoda', nsteps=mxstep)
    def RUN(params):
        solver.run(params)
if run == 'cupSODA':
    set_cupSODA_path("/home/pinojc/CUPSODA")
    output = "model,nsims,tpb,mem,pythontime,rtol,atol,mxsteps,t_end,n_steps,deterministic,vol,card\n"

def main():
    
    global output
    if run == "cupSODA":
        global c_matrix , MX_0
        num_particles = len(MX_0)
        mem = 2
        i = 8
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
        
        for i in range(num_particles):
            np.savetxt('output_%s.txt'%str(i),solver.yobs[i])
    if run == 'scipy':

        c_matrix = np.zeros((num_particles, len(nominal_values)))
        c_matrix[:,:] = nominal_values
        pool = multiprocessing.Pool(processes=num_processes)
        Start = time.time()
        pool.map(RUN,c_matrix)
        TIME = time.time() - Start
        print 'sim = %s , time = %s sec' % (number_particles,TIME)
        output +='%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n'% \
                (name,num_particles,TIME,RTOL,ATOL,mxstep,np.max(tspan),len(tspan),CPU,GHZ,str(num_processes))

if multi == True:
    for i in [1,2,4,8,16,32]:
        print i
        num_processes = i
        for j in simulations:
            main(j)
else:
    main()
#print output
#outFile = open(sys.argv[2],'w')
#outFile.write(output)
#outFile.close()
image = np.zeros((len(proteins_of_interest)*len(vals),len(proteins_of_interest)*len(vals)))
X = 1
Y = 1
done = []
y = 0
x = 0
counter = 0
# for i in range(len(proteins_of_interest)):
#     y = i*len(vals)
#     for j in range(len(proteins_of_interest)):
#         x = j *len(vals)
#         if x >= y:
#             continue
#         for a in range(len(vals)):
#             for b in range(len(vals)):
#                 if   b >= a :
#                     continue
#                 elif  b < a :
#                     tmp = np.loadtxt('output_%s.txt'%str(counter))
#                     counter += 1
#                     tmp = likelihood(tmp[:,2])
#                     image[x+a,y+b] = tmp
counter = 0
started = False
for i in range(len(proteins_of_interest)):
    y = i*len(vals)
    for j in range(len(proteins_of_interest)):
        x = j *len(vals)
        for a in range(len(vals)):
            for b in range(len(vals)):
                tmp = np.loadtxt('output_%s.txt'%str(counter))
                counter += 1      
                if started == False:
                    tbid = tmp[:,0]
                    smac = tmp[:,1]
                    cparp = tmp[:,2]
                    started = True
                    plt.plot(tmp)
                    plt.show()
                else:
                    tbid = np.column_stack((tbid,tmp[:,0]))
                    smac = np.column_stack((smac,tmp[:,1]))
                    cparp = np.column_stack((cparp,tmp[:,2]))
  
plt.plot(tbid)
plt.show()
blank = np.zeros(np.shape(smac)[1])  
for i in range(np.shape(smac)[1]):
    #blank[i] = likelihood(smac[:,i])
    blank[i] = like(tbid[-1,i])

counter=0
for i in range(len(proteins_of_interest)):
    y = i*len(vals)
    for j in range(len(proteins_of_interest)):
        x = j *len(vals)
        if x == y:
            for a in range(len(vals)):
                for b in range(len(vals)):
                    counter+=1
            continue
        for a in range(len(vals)):
            for b in range(len(vals)):
                tmp = likelihood(smac[:,i])
                image[x+a,y+b] = tmp
                counter += 1
plt.imshow(image,interpolation='nearest',origin='lower',cmap=plt.get_cmap('bwr'))
plt.xticks(np.arange(len(vals)/2,len(image),len(vals)),proteins_of_interest,rotation='vertical')
plt.yticks(np.arange(len(vals)/2,len(image),len(vals)),proteins_of_interest)
plt.grid(True,which='minor')
plt.colorbar()
plt.savefig('earm_sensitivity.png')
plt.show()   





