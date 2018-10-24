import numpy as np
import camb
import time

bad_chisq=1e100
param_list=['ombh2','omch2','H0','tau','As','ns']
np_power=2 #number of params to set via InitPower instead of set_cosmology

def update_model(params,cosmology):
    np=len(param_list)
    assert(np==len(cosmology))
    p2=params.copy()
    np_normal=np-np_power
    p2.set_cosmology(ombh2=cosmology[0],omch2=cosmology[1],H0=cosmology[2],tau=cosmology[3])
    p2.InitPower.set_params(As=cosmology[4],ns=cosmology[5])
    return p2

def wmap_chisq(cosmology,wmap,pars_in):
    t0=time.time()
    try:
        pars=update_model(pars_in,cosmology)
        results=camb.get_results(pars)
        power=results.get_cmb_power_spectra(pars,CMB_unit='muK')['total']
    except:
        return bad_chisq
    t1=time.time()
    inds=np.asarray(wmap[:,0],dtype='int')
    pred=power[inds,0]
    chisq=np.sum( (pred-wmap[:,1])**2/wmap[:,2]**2)
    t2=time.time()
    #print t1-t0,t2-t1
    return chisq

pars=camb.CAMBparams()
#mycosmo=[0.02, 0.10, 70,0.05,2e-9,0.98]
mycosmo=[2.28e-02,   0.11177,   70.68, 4.57e-02,   2.0138e-09,   9.763e-01]
mycosmo=np.asarray(mycosmo)

par2=update_model(pars,mycosmo)
t1=time.time()
results=camb.get_results(par2)
t2=time.time()

errs=np.asarray([2.4e-4, 7.34e-4,0.7,1.8e-3,7e-12,0.005])/2
wmap=np.loadtxt('wmap_tt_spectrum_9yr_v5.txt')
chisq=wmap_chisq(mycosmo,wmap,pars)

ombh=np.linspace(0.022,0.024,50)
cosmo_use=mycosmo.copy()


nsamp=1000
chains=np.zeros([nsamp,1+len(mycosmo)])
for iter in range(nsamp):
    new_cosmo=mycosmo+np.random.randn(len(errs))*errs
    new_chisq=wmap_chisq(new_cosmo,wmap,pars)
    accept=False
    if new_chisq<bad_chisq:
        thresh=np.exp(-0.5*(new_chisq-chisq))
        if np.random.rand()<thresh:
            accept=True
    print iter, ' new_chisq is',new_chisq, ' and accept ',accept,new_cosmo
    if accept:
        chisq=new_chisq
        mycosmo=new_cosmo
    chains[iter,0]=chisq
    chains[iter,1:]=mycosmo
