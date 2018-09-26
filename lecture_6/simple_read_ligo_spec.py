import numpy
import numpy as np
from matplotlib import pyplot as plt
import h5py
import glob
plt.ion()

def read_template(filename):
    dataFile=h5py.File(filename,'r')
    template=dataFile['template']
    th=template[0]
    tl=template[1]
    return th,tl
def read_file(filename):
    dataFile=h5py.File(filename,'r')
    dqInfo = dataFile['quality']['simple']
    qmask=dqInfo['DQmask'][...]

    meta=dataFile['meta']
    gpsStart=meta['GPSstart'].value
    #print meta.keys()
    utc=meta['UTCstart'].value
    duration=meta['Duration'].value
    strain=dataFile['strain']['Strain'].value
    dt=(1.0*duration)/len(strain)

    dataFile.close()
    return strain,dt,utc



#fnames=glob.glob("[HL]-*.hdf5")
#fname=fnames[0]
fname='H-H1_LOSC_4_V2-1126259446-32.hdf5'
print 'reading file ',fname
strain,dt,utc=read_file(fname)

#th,tl=read_template('GW150914_4_template.hdf5')
template_name='GW150914_4_template.hdf5'
th,tl=read_template(template_name)


#spec,nu=measure_ps(strain,do_win=True,dt=dt,osamp=16)
#strain_white=noise_filter(strain,numpy.sqrt(spec),nu,nu_max=1600.,taper=5000)

#th_white=noise_filter(th,numpy.sqrt(spec),nu,nu_max=1600.,taper=5000)
#tl_white=noise_filter(tl,numpy.sqrt(spec),nu,nu_max=1600.,taper=5000)


#matched_filt_h=numpy.fft.irfft(numpy.fft.rfft(strain_white)*numpy.conj(numpy.fft.rfft(th_white)))
#matched_filt_l=numpy.fft.irfft(numpy.fft.rfft(strain_white)*numpy.conj(numpy.fft.rfft(tl_white)))

x=np.linspace(-1,1,len(strain))*np.pi
win=0.5+0.5*np.cos(x)
tot_len=dt*len(strain)
dnu=1.0/tot_len

datft=np.fft.rfft(win*strain)
nuvec=np.arange(len(datft))*dnu
plt.clf();
plt.loglog(nuvec[1:],np.abs(datft[1:])**2)
plt.savefig('ligo_spec_win_wfreq.png')
plt.axis([40,2000,1e-41,1e-31])
plt.savefig('ligo_spec_win_wfreq_zoom.png')

