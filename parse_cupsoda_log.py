# -*- coding: utf-8 -*-
"""
Created on Sun May 10 20:56:38 2015

@author: James C. Pino
"""
import sys
import numpy as np

print "To parse we need a log file (verbose output from cupSODA) and results, created in pysb run"
print "python parse_cupsoda_log.py new_tyson_log.txt new_tyson_results.log TYSON.txt"
File = open(sys.argv[1],'r')
log = np.genfromtxt(sys.argv[2],delimiter=',',dtype='string')
print np.shape(log)
tmp = np.zeros((np.shape(log)[0],1))

count = 1
for line in File:
    print line
    if line.startswith("Running time"):
        tmp[count]=line.split()[2]
        count+=1

output = "model,nsims,tpb,mem,pythontime,rtol,atol,mxsteps,t_end,n_steps,deterministic,vol,card\n"
end = np.column_stack((log[:,:5],tmp))
end = np.column_stack((end,log[:,5:]))
end[0,5]='cupsodatime'


print end
np.savetxt(sys.argv[3],end,delimiter=",", fmt="%s")


File.close()
