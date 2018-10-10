import numpy as np
import simple_read_ligo as rl


def make_ft_vec(n):
    #make a routine to make a vector that goes from
    #0 to n/2, then -n/2 up to -1
    x=np.arange(n)
    x[n/2:]=x[n/2:]-n
    assert(x[-1]==-1)
    return x*1.0  #make sure we return it as a float, not an int

def smooth_dat_gauss(dat,fwhm):
    n=len(dat)
    x=make_ft_vec(n)
    sig=fwhm/np.sqrt(8*np.log(2))
    smooth_vec=np.exp(-0.5*x**2/sig**2)
    smooth_vec=smooth_vec/smooth_vec.sum()
    datft=np.fft.rfft(dat)
    vecft=np.fft.rfft(smooth_vec)
    dat_smoothft=datft*vecft
    dat_smooth=np.fft.irfft(dat_smoothft,n=n)
    return dat_smooth

fname='H-H1_LOSC_4_V2-1126259446-32.hdf5'
print 'reading file ',fname
strain,dt,utc=rl.read_file(fname)

template_name='GW150914_4_template.hdf5'
th,tl=rl.read_template(template_name)

#make a window
n=len(strain)
x=np.linspace(-1.0,1.0,n)*np.pi
win=0.5+0.5*np.cos(x)

strainft=np.fft.rfft(strain*win)
ps_raw=np.abs(strainft)**2
ps_smooth=smooth_dat_gauss(ps_raw,10)

dnu=1/(n*dt)
nu=np.arange(len(ps_smooth))*dnu

ninv=1.0/ps_smooth
ninv[nu<10]=0
ninv[nu>1600]=0 #because it seems to crash there
dat_whitened_ft=strainft*np.sqrt(ninv)



dat_whitened=np.fft.irfft(dat_whitened_ft,n)

thft=np.fft.rfft(th*win)
thft_whitened=thft*np.sqrt(ninv)
th_whitened=np.fft.irfft(thft_whitened)
top=np.fft.irfft(np.conj(thft_whitened)*dat_whitened_ft)

#ninv[nu>


