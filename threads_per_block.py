import numpy as np
import pylab as plt



quit()

#earm = np.loadtxt('to_export_earm.csv',dtype='string',delimiter=',')
#earm_data = earm[1:,:].astype('float')

ras = np.loadtxt('threads_per_block_ras.csv',delimiter=',',skiprows=1)
ras[ras==0]=np.nan
def plot_ras():
    titles = [10,100,1000,10000,100000]
    for n,i in enumerate([0,4,8,12,16]):
        plt.plot(ras[:,i],ras[:,i+1],'o-',label='global')
        plt.plot(ras[:,i],ras[:,i+2],'o-',label='shared')
        plt.plot(ras[:,i],ras[:,i+3],'o-',label='both')
        plt.title('%s simulations'%titles[n])    
        plt.legend(loc=0)
        plt.ylabel('Time(s)')
        plt.xlabel('Threads/block')
        plt.savefig('ras_%s.png'%titles[n])
        plt.show()
    
plot_ras()

tyson = np.genfromtxt('threads_per_block_tyson.csv',delimiter=',',skiprows=1)
tyson[tyson==0]=np.nan

def plot_tyson():
    titles = [10,100,1000,10000,100000]
    for n,i in enumerate([0,4,8,12,16]):
        plt.plot(tyson[:,i],tyson[:,i+1],'o-',label='global')
        plt.plot(tyson[:,i],tyson[:,i+2],'o-',label='shared')
        plt.plot(tyson[:,i],tyson[:,i+3],'o-',label='both')
        plt.title('%s simulations'%titles[n])    
        plt.legend(loc=0)
        plt.ylabel('Time(s)')
        plt.xlabel('Threads/block')
        plt.savefig('tyson_%s.png'%titles[n])
        plt.show()
plot_tyson()
        
earm = np.genfromtxt('threads_per_block_earm.csv',delimiter=',',skiprows=1)
earm[earm==0]=np.nan        
def plot_earm():
    titles = [10,100,1000,10000]
    for n,i in enumerate([0,4,8,12]):
        plt.plot(earm[:,i],earm[:,i+1],'o-',label='global')
        plt.plot(earm[:,i],earm[:,i+2],'o-',label='shared')
        plt.plot(earm[:,i],earm[:,i+3],'o-',label='both')
        plt.title('%s simulations'%titles[n])    
        plt.legend(loc=0)
        plt.ylabel('Time(s)')
        plt.xlabel('Threads/block')
        plt.savefig('earm_%s.png'%titles[n])
        plt.show()
plot_earm()        