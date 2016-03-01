import time
import pysb.integrate as integrate
import pysb
import numpy as np
import re
import os
from pysb_cupsoda import set_cupsoda_path, CupsodaSolver
from pysb.bng import generate_equations
import sys
import multiprocessing

#name = sys.argv[1]
name = 'earm'
if name == 'ras':
    from ras_amp_pka import model

    tspan = np.linspace(0, 1500, 100)
    simulations = [10, 100, 1000, 10000, 100000]
    simulations = [2**8,2**9,2**10, 2**11,2**12,2**13,2**14,2**15,2**16,2**17]
elif name == 'tyson':
    from pysb.examples.tyson_oscillator import model
    tspan = np.linspace(0, 100, 100)
    simulations = [10, 100, 1000, 10000, 100000]
    simulations = [10, 100, 1000]
    simulations = [2**8,2**9,2**10, 2**11,2**12,2**13,2**14,2**15,2**16,2**17]
elif name == 'earm':
    from earm.lopez_embedded import model
    tspan = np.linspace(0, 20000, 100)
    #simulations = [ 1000, 5000, 10000, 15000]
    simulations = [10, 100, 1000, 10000]#,100000]
    simulations = [2**8,2**9,2**10, 2**11,2**12,2**13,2**14,2**15,2**16,2**17]

run = 'cupSODA'
multi = False
#run ="scipy"

generate_equations(model)

ATOL = 1e-6
RTOL = 1e-6
mxstep = 20000
det = 1
vol = 0
# card = 'K20c'
# card = 'gtx760'
card = 'gtx980-diablo'

# puma
CPU = 'Intel-Xeon-E5-2687W-v2'
GHZ = '3.40'

# macbook
# CPU = 'Intel-core-I5'
# GHZ = '2.5'


params_names = [p.name for p in model.parameters_rules()]
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
par_dict = {par_names[i]: i for i in range(len(par_names))}
print len(rate_params)
print len(model.reactions)
quit()
if run == 'scipy':
    global output
    output = "model,nsims,scipytime,rtol,atol,mxsteps,t_end,n_steps,cpu,GHz,num_cpu\n"
    solver = pysb.integrate.Solver(model, tspan, rtol=RTOL, atol=ATOL, integrator='lsoda', )


    def run_pool(params):
        solver.run(params)
if run == 'cupSODA':
    set_cupsoda_path("/home/pinojc/git/cupSODA")
    solver = CupsodaSolver(model, tspan, atol=ATOL, rtol=RTOL, verbose=False)
    output = "model,nsims,tpb,mem,cupsodatime,pythontime,rtol,atol,mxsteps,t_end,n_steps,deterministic,vol,card\n"


def main(number_particles):
    num_particles = int(number_particles)
    c_matrix = np.zeros((num_particles, len(nominal_values)))
    c_matrix[:, :] = nominal_values
    global output
    if run == "cupSODA":
        mx_0 = np.zeros((num_particles, len(model.species)))
        for i in xrange(len(model.initial_conditions)):
            for j in xrange(len(model.species)):
                if str(model.initial_conditions[i][0]) == str(model.species[j]):
                    x = model.initial_conditions[i][1]
                    mx_0[:, j] = x.value
                    break
        #for mem in [0, 1, 2]:
        #    for i in (1, 2, 4, 8, 16, 32, 64):
        #        if num_particles < i:
        #            continue
        #        if i == 64 and name == 'earm':
        #            continue
        mem = 2
        i = 32
        start_time = time.time()
        solver.run(c_matrix,
                   mx_0,
                   n_blocks=np.ceil(num_particles / i),
                   outdir=os.path.join('/tmp/ramdisk', 'CUPSODA_%s') % model.name,
                   gpu=0,
                   max_steps=mxstep,
                   load_conc_data=False,
                   memory_usage=mem)
        end_time = time.time()
        new_line = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % \
                  (name, num_particles, str(i), mem,solver.time, end_time - start_time, RTOL, ATOL, len(tspan), mxstep,
                   np.max(tspan), det, vol, card)
        output += new_line
        print new_line
        """
        t = tspan
        x = solver.yobs["YT"]

        import pylab as plt
        plt.plot(t, solver.yobs["YT"][0,:], lw=2, label='YT')
        plt.plot(t, solver.yobs["M"][0,:], lw=2, label='M')

        plt.legend(loc=0)
        plt.xlabel('time')
        plt.ylabel('population')

        plt.show()
        """
        print 'out==', solver.yobs[0][0], solver.yobs[0][-1], '==out'
        os.system('rm -r %s' % os.path.join('/tmp/ramdisk', 'CUPSODA_%s') % model.name)
        #print 'removed directory'

    if run == 'scipy':

        pool = multiprocessing.Pool(processes=num_processes)
        start_time = time.time()
        #pool.map(run_pool, c_matrix)
        for i in xrange(number_particles):
            solver.run()
        total_time = time.time() - start_time
        print 'sim = %s , time = %s sec' % (number_particles, total_time)
        output += '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % \
                  (name,
                   num_particles,
                   total_time,
                   RTOL,
                   ATOL,
                   mxstep,
                   np.max(tspan),
                   len(tspan),
                   CPU,
                   GHZ,
                   str(num_processes))


if multi:
    for i in [1, 2, 4, 8, 16, 32]:
        print i
        num_processes = i
        for j in simulations:
            main(j)
else:
    num_processes = 1
    main(2**13)
    #for j in simulations:
    #    main(j)
print output
#outFile = open(sys.argv[2], 'w')
#outFile.write(output)
#outFile.close()

