import numpy as np
from matplotlib import pyplot as plt

chains=np.loadtxt('wmap9_pol_corr_chains_v4_hippo.txt')
for i in range(1,chains.shape[1]):
    dat=chains[:,i].copy()
    myval=np.mean(dat)
    mystd=np.std(dat)
    dat=dat-myval
    datft=np.fft.rfft(dat)
    mycorr=np.fft.irfft(datft*np.conj(datft))
    mylen=np.min(np.where(mycorr<0))
    nsamp=len(dat)/mylen
    print 'mean/error on parameter ',i-1,' are ',myval,mystd,' with roughly ',nsamp,' indepenent samples'

chains_norm=chains[:,1:].copy()
for i in range(chains_norm.shape[1]):
    chains_norm[:,i]=chains_norm[:,i]-chains_norm[:,i].mean()
    chains_norm[:,i]=chains_norm[:,i]/chains_norm[:,i].std()
mycorr=np.dot(chains_norm.transpose(),chains_norm)/chains_norm.shape[0]
print 'correlation matrix is: '
print mycorr

plt.ion()
plot_thresh=0.5
npar=mycorr.shape[0]
print 'npar is',npar
for i in range(npar):
    for j in range(i+1,npar):
        if np.abs(mycorr[i,j])>plot_thresh:
            plt.clf();
            plt.plot(chains[:,i+1],chains[:,j+1],'.')
            outname='wmap_corrs_'+repr(i)+'_'+repr(j)+'.png'
            plt.savefig(outname)
    
