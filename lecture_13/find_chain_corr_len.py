import numpy as np

#chains=np.loadtxt('chain_tt_kiyo_dmpow_wcov.txt_4')
chains=np.loadtxt('chain_tt_kiyo_dmpow_decorr_wcov.txt_4')
for i in range(3,chains.shape[1]):
    vec=chains[:,i].copy()
    vec=vec-vec.mean()
    vecft=np.fft.rfft(vec)
    mycorr=np.fft.irfft(vecft*np.conj(vecft))
    ii=np.min(np.where(mycorr<0))
    print 'for parameter ',i,' correlation goes to zero at ',ii
