import time
import pysb.integrate
import pysb
import numpy as np
import pylab as plt
import os
from pysb_cupsoda import set_cupsoda_path,CupsodaSolver
from pysb.bng import generate_equations
import scipy.interpolate
import multiprocessing
from pysb.integrate import odesolve
from earm.lopez_embedded import model
from pysb.util import update_param_vals,load_params
import re

tspan = np.linspace(0, 20000,100)

new_params = load_params('Params_486/pars_embedded_486.txt')
#new_params = load_params('Params_911/pars_embedded_911.txt')
savename = 'parameters_486'

update_param_vals(model, new_params)

set_cupsoda_path("/home/pinojc/CUPSODA")
run = 'cupSODA'
#run = 'scipy'

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
#proteins_of_interest = ['L_0','Mcl1_0','Bak_0','Bax_0']
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

            
def likelihood1(initial_tmp):
    solver = pysb.integrate.Solver(model, tspan,rtol=RTOL, atol=ATOL, integrator='lsoda', mxstep=mxstep)
    solver.run(param_values=initial_tmp)
    ysim_momp = solver.yobs['cPARP']
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
tod = likelihood1(initial_tmp)
print tod

def likelihood(ysim_momp):

    if np.nanmax(ysim_momp) == 0:
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


#vals = [.5,1.,10.]
vals = np.hstack((np.linspace(.7,.9,3),np.logspace(0,.3,5)))
print vals

#size_of_matrix = (len(proteins_of_interest)*len(proteins_of_interest))*len(vals)*len(vals) 
size_of_matrix = (len(proteins_of_interest)*(len(proteins_of_interest)-1))*len(vals)*len(vals)/2 

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
index_of_species_of_interest = {}
MX_0 = np.zeros((size_of_matrix,len(model.species)))
for i in xrange(len(model.initial_conditions)):
    for j in xrange(len(model.species)):
        if str(model.initial_conditions[i][0]) == str(model.species[j]):
            x = model.initial_conditions[i][1].value
            MX_0[:,j] = x
            if model.initial_conditions[i][1].name in proteins_of_interest: 
                index_of_species_of_interest[model.initial_conditions[i][1].name] = j


counter = 0
done = []
for i  in proteins_of_interest:
    for j in proteins_of_interest:
        if j in done:
            continue
        if i == j:
            continue
        for a,c in enumerate(vals):
            for b,d in enumerate(vals):
                x = index_of_species_of_interest[i]
                y = index_of_species_of_interest[j]
                MX_0[counter,x] *= c
                MX_0[counter,y] *= d
                counter+=1
    done.append(i)    
            
print("Number of simulations to run = %s" % counter)


def main():
    
    if run == "cupSODA":
        
        global c_matrix , MX_0
        num_particles = len(MX_0)
        mem = 2
        i = 8
        solver = CupsodaSolver(model, tspan, atol=ATOL, rtol=RTOL, verbose=False)
        Start = time.time()
        quit()
        solver.run(c_matrix, 
                   MX_0 ,
                   n_blocks = np.int(num_particles/i), 
                   outdir=os.path.join('.','CUPSODA_%s')%model.name, 
                   gpu = 0,
                   max_steps = mxstep,
                   load_conc_data = False,
                   memory_usage = mem)
        print 'out==',solver.yobs[0][0],solver.yobs[0][-1],'==out'
        os.system('rm -r %s'%os.path.join('.','CUPSODA_%s')%model.name)
        print 'removed directory'
        
        for i in range(num_particles):
            np.savetxt('output_%s.txt'%str(i),solver.yobs[i])
    if run == 'scipy':
        solver = pysb.integrate.Solver(model, tspan,rtol=RTOL, atol=ATOL, integrator='lsoda', mxstep=mxstep)
        print("started at %s"%time.time())
        Start = time.time()
        for i in range(size_of_matrix):
            solver.run(y0 = MX_0[i,:])
            np.savetxt('test-output_%s.txt'%str(i),solver.yobs)
        TIME = time.time() - Start
        print 'sim = %s , time = %s sec' % (size_of_matrix,TIME)
main()

tbid = np.zeros((len(tspan),size_of_matrix))
smac = np.zeros((len(tspan),size_of_matrix))
cparp = np.zeros((len(tspan),size_of_matrix))
counter = 0
for i in range(len(proteins_of_interest)):
    for j in range(i,len(proteins_of_interest)):
        if i == j:
            continue
        for a in range(len(vals)):
            for b in range(len(vals)):
                tmp = np.loadtxt('test-output_%s.txt'%str(counter))
                tbid[:,counter] = tmp[:,0]
                smac[:,counter] = tmp[:,1]
                cparp[:,counter] = tmp[:,2]
                counter += 1      

image = np.zeros((len(proteins_of_interest)*len(vals),len(proteins_of_interest)*len(vals)))
counter=0

for i in range(len(proteins_of_interest)):
    y = i*len(vals)
    for j in range(i,len(proteins_of_interest)):
        x = j *len(vals)
        if x == y:
            continue
        for a in range(len(vals)):
            for b in range(len(vals)):
                tmp = likelihood(cparp[:,counter])
                image[y+a,x+b] = tmp
                image[x+b,y+a] = tmp
                counter += 1
print counter
"""    #Good            
image = np.zeros((len(proteins_of_interest)*len(vals),len(proteins_of_interest)*len(vals)))
counter=0

for i in range(len(proteins_of_interest)):
    y = i*len(vals)
    #counter += len(vals)*len(vals)*i
    for j in range(i,len(proteins_of_interest)):
        x = j *len(vals)
        if x == y:
            for a in range(len(vals)):
                for b in range(len(vals)):
                    counter+=1
            continue
        for a in range(len(vals)):
            for b in range(len(vals)):
                tmp = likelihood(cparp[:,counter])
                image[x+b,y+a] = tmp
                image[y+a,x+b] = tmp
                counter += 1
print counter
"""
plt.imshow(image,interpolation='nearest',origin='lower',cmap=plt.get_cmap('bwr'))
plt.xticks(np.arange(len(vals)/2,len(image),len(vals)),proteins_of_interest,rotation='vertical')
plt.yticks(np.arange(len(vals)/2,len(image),len(vals)),proteins_of_interest)
plt.grid(True,which='minor')
plt.colorbar()
plt.savefig('earm_heatplot_%s.pdf'%savename)
plt.show()  
    
all = []
for i in range(0,len(image),len(vals)):
    tmp = image[:,i:i+len(vals)].flatten()
    print np.shape(tmp),
    tmp = tmp[tmp!=0]
    print np.shape(tmp)
    all.append(tmp)
#plt.figure(figsize=(3,3))
plt.boxplot(all,vert=False,labels=proteins_of_interest,showfliers=False)
plt.xlabel('Time of death change(%)')
plt.savefig('earm_boxplot_%s.pdf'%savename)
plt.show()