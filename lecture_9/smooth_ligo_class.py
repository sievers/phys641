import numpy as np
import simple_read_ligo as rl


def smooth_data(vec,fwhm):
    n=len(vec)
    x=np.arange(n)
    x[n/2:]=x[n/2:]-n
    print x[0],x[-1]
    sig=fwhm/np.sqrt(8*np.log(2))
    y=np.exp(-0.5*x**2/sig**2)
    y=y/np.sum(y)
    
    vecft=np.fft.rfft(vec)
    yft=np.fft.rfft(y)
    vec_smooth=np.fft.irfft(yft*vecft,n)
    return vec_smooth


fname='H-H1_LOSC_4_V2-1126259446-32.hdf5'
print 'reading file ',fname
strain,dt,utc=rl.read_file(fname)

template_name='GW150914_4_template.hdf5'
tl,th=rl.read_template(template_name)


n=len(strain)
#x=np.arange(n)
#x[n/2:]=x[n/2:]-n

x=np.linspace(-1,1,n)
win=0.5*np.cos(np.pi*x)+0.5

Fk=np.abs(np.fft.rfft(win*strain))**2


N=smooth_data(Fk,20)
Nmhalf=1/np.sqrt(N)
Nmhalf[:200]=0
Nmhalf[53000:]=0

strainft=np.fft.rfft(win*strain)
strainft_white=Nmhalf*strainft
strain_white=np.fft.irfft(strainft_white,len(strain))


th_ft=np.fft.rfft(th*win)
thft_white=th_ft*Nmhalf
th_white=np.fft.irfft(thft_white,len(th))

mf_ft=strainft_white*np.conj(thft_white)
mf=np.fft.irfft(mf_ft,n)
