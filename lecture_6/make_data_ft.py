import numpy as np

n=50
niter=1000
mycov=np.ones([n,n])+np.eye(n) #example on PS
mycorr=mycov[0,:] #kind of by definition, if stationary

mycorr_ft=np.fft.rfft(mycorr)
dat_white=np.random.randn(niter,n)
dat_white_ft=np.fft.rfft(dat_white,axis=1)
sqrt_mycorr_ft_mat=np.repeat([np.sqrt(np.abs(mycorr_ft))],niter,axis=0)
dat_corr_ft=dat_white_ft*sqrt_mycorr_ft_mat
dat_corr=np.fft.irfft(dat_corr_ft)


dat_cov=np.dot(dat_corr.transpose(),dat_corr)/niter
