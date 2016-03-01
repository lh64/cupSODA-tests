import os
import multiprocessing
import time
import matplotlib
import numpy as np
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import pysb
from pysb.bng import generate_equations
from pysb.integrate import odesolve
from pysb_cupsoda import set_cupsoda_path, CupsodaSolver

option = 'tyson'
if option == 'ras':
    from ras_amp_pka import model
    tspan = np.linspace(0, 1000, 301)
    observable = 'obs_cAMP'
    savename = 'ras'

if option == 'tyson':
    from pysb.examples.tyson_oscillator import model
    tspan = np.linspace(0, 200, 1001)
    observable = 'Y3'
    savename = 'tyson'

set_cupsoda_path("/home/pinojc/git/cupSODA")
run = 'cupSODA'
#run = 'scipy'

ATOL = 1e-6
RTOL = 1e-6
mxstep = 20000
det = 1
vol = 0
card = 'K20c'
# puma
CPU = 'Intel-Xeon-E5-2687W-v2'
GHZ = '3.40'

generate_equations(model)
proteins_of_interest = []
for i in model.initial_conditions:
    proteins_of_interest.append(i[1].name)
if option == 'tyson':
    proteins_of_interest.remove('__source_0')

key_list = []
initial_tmp = dict()
for species, keys in model.initial_conditions:
    for j in proteins_of_interest:
        if keys.name == j:
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
par_dict = {par_names[i]: i for i in range(len(par_names))}
par_vals = np.array([model.parameters[nm].value for nm in par_names])


def likelihood_one(parameters):
    solver = pysb.integrate.Solver(model, tspan, rtol=RTOL, atol=ATOL,
                                   integrator='lsoda', mxstep=mxstep)
    solver.run(param_values=parameters)
    out = solver.yobs[observable]
    if option == 'ras':
        plt.figure()

        plt.plot(solver.tspan, solver.yobs[observable], linewidth=2, )
        max_y = np.where(out == out.max())

        plt.plot(solver.tspan[max_y], out.max(), 'o', color='red', markersize=14,mew=3,mfc='none', alpha=.75)
        plt.vlines(solver.tspan[max_y], 0,out.max()+1000000,linestyles='dashed')
        plt.ylabel('cAMP molecules',fontsize = 16)
        plt.xlabel('Time (s)',fontsize = 16)
        plt.ylim(0,125000)
        plt.savefig(observable + '.png')
        plt.savefig(observable + '.eps')
        plt.close()
        return out.max()
    elif option == 'tyson':
        plt.figure()
        plt.title(observable)
        plt.plot(solver.tspan, out)
        plt.savefig(observable + '.png')
        plt.close()

        timestep = solver.tspan[:-1]
        y = out[:-1] - out[1:]
        times = []
        prev = y[0]
        for n in range(1, len(y)):
            if y[n] > 0 > prev:
                times.append(timestep[n])
            prev = y[n]
        plt.plot(timestep, y)
        plt.savefig('slope.png')
        times = np.array(times)

        #print times[1] - times[0]
        return times[1] - times[0]


times = likelihood_one(initial_tmp)
print times



def likelihood(ysim_momp):
    out = ysim_momp
    if option == 'ras':
        #print out.max()
        return (out.max() - times) / times * 100.
    else:
        timestep = tspan[:-1]
        y = out[:-1] - out[1:]
        freq = 0
        local_times = []
        prev = y[0]
        for n in range(1, len(y)):
            if y[n] > 0 > prev:
                local_times.append(timestep[n])
                freq += 1
            prev = y[n]
        local_times = np.array(local_times)
        print local_times
        print np.average(local_times)/len(local_times)
        local_times = local_times[1] - local_times[0]
        print local_times
        return (times - local_times) / times * 100.


vals = np.logspace(-.1, .1, 20)

n_sam = len(vals)
n_proteins = len(proteins_of_interest)

size_of_matrix = (n_proteins * n_proteins - n_proteins) * (n_sam * n_sam) / 2

c_matrix = np.zeros((size_of_matrix, len(nominal_values)))
c_matrix[:, :] = nominal_values

MX_0 = np.zeros((size_of_matrix, len(model.species)))

index_of_species_of_interest = {}

for i in xrange(len(model.initial_conditions)):
    for j in xrange(len(model.species)):
        if str(model.initial_conditions[i][0]) == str(model.species[j]):
            x = model.initial_conditions[i][1].value
            MX_0[:, j] = x
            if model.initial_conditions[i][1].name in proteins_of_interest:
                index_of_species_of_interest[model.initial_conditions[i][1].name] = j

counter = 0
done = []
for i in proteins_of_interest:
    for j in proteins_of_interest:
        if j in done:
            continue
        if i == j:
            continue
        for a, c in enumerate(vals):
            for b, d in enumerate(vals):
                x = index_of_species_of_interest[i]
                y = index_of_species_of_interest[j]
                MX_0[counter, x] *= c
                MX_0[counter, y] *= d
                counter += 1
    done.append(i)

print("Number of simulations to run = %s" % counter)

if not os.path.exists('OUT'):
    os.mkdir('OUT')
solver = pysb.integrate.Solver(model, tspan, rtol=RTOL, atol=ATOL,
                               integrator='lsoda', mxstep=mxstep)

def run_pool(k):
    solver.run(y0=MX_0[k, :])
    np.savetxt('OUT/test-output_%s.txt' % str(k), solver.yobs[observable])

def main():
    if run == "cupSODA":
        global c_matrix, MX_0
        num_particles = len(MX_0)
        mem = 2
        threads_per_block = 8
        solver = CupsodaSolver(model, tspan, atol=ATOL, rtol=RTOL, verbose=False)
        start_time = time.time()
        solver.run(c_matrix,
                   MX_0,
                   #n_blocks=np.int(num_particles / threads_per_block),
                   outdir=os.path.join('.', 'CUPSODA_%s') % model.name,
                   gpu=0,
                   max_steps=mxstep,
                   load_conc_data=False,
                   memory_usage=mem)
        print 'out==', solver.yobs[0][0], solver.yobs[0][-1], '==out'
        #os.system('rm -r %s' % os.path.join('.', 'CUPSODA_%s') % model.name)
        plt.plot(solver.tspan,solver.yobs[:][observable].T)
        plt.savefig('test.png')
        for threads_per_block in range(num_particles):
            np.savetxt('OUT/test-output_%s.txt' % str(threads_per_block), solver.yobs[observable][threads_per_block])
    if run == 'scipy':


        pool = multiprocessing.Pool(processes=4)
        start_time = time.time()
        points = range(size_of_matrix)
        pool.map(run_pool,points )
        #for k in range(size_of_matrix):
            #solver.run(y0=MX_0[k, :])
            #np.savetxt('OUT/test-output_%s.txt' % str(k), solver.yobs[observable])
        time_taken = time.time() - start_time
        print 'sim = %s , time = %s sec' % (size_of_matrix, time_taken)


#main()


def load_results():
    camp = np.zeros((len(tspan), size_of_matrix))
    counter = 0
    for i in range(len(proteins_of_interest)):
        for j in range(i, len(proteins_of_interest)):
            if i == j:
                continue
            for a in range(len(vals)):
                for b in range(len(vals)):
                    tmp = np.loadtxt('OUT/test-output_%s.txt' % str(counter))
                    camp[:, counter] = tmp
                    counter += 1
    return camp


camp = load_results()
image = np.zeros((len(proteins_of_interest) * len(vals), len(proteins_of_interest) * len(vals)))
counter = 0
for i in range(len(proteins_of_interest)):
    y = i * len(vals)
    for j in range(i, len(proteins_of_interest)):
        x = j * len(vals)
        if x == y:
            continue
        for a in range(len(vals)):
            for b in range(len(vals)):
                tmp = likelihood(camp[:, counter])
                image[y + a, x + b] = tmp
                image[x + b, y + a] = tmp
                counter += 1
print("Number of simulations to ran = %s" % counter)

vmax = max(np.abs(image.min()), image.max())
vmin = -1 * vmax
np.savetxt('sens_%s_matrix.csv' % savename, image)
plt.figure()
print image.shape
plt.imshow(image[20:,:20], interpolation='nearest', origin='lower', cmap=plt.get_cmap('bwr'), vmin=vmin, vmax=vmax)
#plt.xticks(np.arange(len(vals) / 2, len(image), len(vals)), proteins_of_interest, rotation='vertical')
#plt.yticks(np.arange(len(vals) / 2, len(image), len(vals)), proteins_of_interest)
#plt.grid(True, which='minor')
plt.colorbar()
plt.savefig('sens_image_matrix_%s.png' % savename, bbox_inches='tight')
plt.show()

all_runs = []
for i in range(0, len(image), len(vals)):
    tmp = image[:, i:i + len(vals)].flatten()
    tmp = tmp[tmp != 0]
    all_runs.append(tmp)

plt.figure()
plt.boxplot(all_runs, vert=False, labels=proteins_of_interest, showfliers=False)
# plt.xlabel('Change in max cAMP (%)')
plt.tight_layout()
plt.savefig('sens_boxplot_%s.png' % savename, bbox_inches='tight')
plt.show()
