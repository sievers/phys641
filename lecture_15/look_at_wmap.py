import numpy as np
from matplotlib import pyplot as plt
import camb

def get_chisq(wmap,ps):
    ell=np.asarray(wmap[:,0],dtype='int')
    delt=wmap[:,1]-ps[ell]
    chisq=np.sum( (delt/wmap[:,2])**2)
    return chisq

def update_model(pars,cosmology):
    p2=pars.copy()
    p2.set_cosmology(ombh2=cosmology[0],omch2=cosmology[1],H0=cosmology[3],tau=cosmology[5])
    p2.InitPower.set_params(As=cosmology[2],ns=cosmology[4])
    return p2

def get_power_spectrum(cosmology,pars):
    pars2=update_model(pars,cosmology)
    results=camb.get_results(pars2)
    powspec=results.get_cmb_power_spectra(pars2,CMB_unit='muK')['total']
    return powspec
    #chisq=get_chisq(wmap,powspec)
    #return chisq

def get_wmap_chisq(cosmology, wmap,pars):
    powspec=get_power_spectrum(cosmology,pars)
    powspec=powspec[:,0]
    return get_chisq(wmap,powspec)

    #order omega_b, omega_cdm, As, h0, n_s, tau


def get_step(step_size):
    npar=len(step_size)
    return np.random.randn(npar)*step_size

def run_mcmc(cosmology,wmap,pars,step_size,nstep=100):
    chisq=get_wmap_chisq(cosmology,wmap,pars)
    chain=np.zeros([nstep,len(cosmology)+1])
    for iter in range(nstep):
        #print iter
        new_cosmology=cosmology+get_step(step_size)        
        new_chisq=get_wmap_chisq(new_cosmology,wmap,pars)

        like=np.exp(-0.5*(new_chisq-chisq))
        accept=np.random.rand()<like
        print iter,chisq,new_chisq,accept,new_cosmology[2],like

        if accept:
            cosmology=new_cosmology
            chisq=new_chisq
        chain[iter,1:]=cosmology
        chain[iter,0]=chisq
    return chain
tt=np.loadtxt('wmap_tt_spectrum_9yr_v5.txt')
plt.ion()

pars=camb.CAMBparams()
cosmology=np.asarray([0.02,0.1157,2.125e-9,70,1,0.05])
chisq=get_wmap_chisq(cosmology,tt,pars)
step_size=np.asarray([0,0.0008,0,0,0,0])

chain=run_mcmc(cosmology,tt,pars,step_size,nstep=150)




#important parameters are:  matter density, baryon density, amplitude, distance to last scattering, n_s, tau

