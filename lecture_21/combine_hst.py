import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits

def apply_mat(map,beam_fts):
    nbeam=len(beams)
    maps_out=[None]*nbeam
    mapft=np.fft.fft2(map)
    for i in range(nbeam):
        maps_out[i]=np.real(np.fft.ifft2(mapft*beam_fts[i]))
    return maps_out
def apply_mat_transpose(maps,beam_fts):
    tot=0
    nmap=len(maps)
    for i in range(nmap):
        tot=tot+np.real(np.fft.ifft2(np.fft.fft2(maps[i])*beam_fts[i]))
    return tot
def apply_noise(maps,noises=None):
    if noises is None:
        return maps
    

def map2map(map,beamfts,noises=None):
    maps=apply_mat(map,beamfts)
    maps=apply_noise(maps,noises)
    map_out=apply_mat_transpose(maps,beamfts)
    return map_out

def run_pcg(x,b,beam_fts,noises=None,niter=40):
    Ax=map2map(x,beam_fts,noises)
    r=b-Ax
    rtr=np.sum(r**2)
    p=r.copy()
    print 'starting residual is ',rtr
    for i in range(niter):
        Ap=map2map(p,beam_fts,noises)
        pAp=np.sum(p*Ap)
        alpha=rtr/pAp
        x_next=x+alpha*p
        r_next=r-alpha*Ap
        rtr_next=np.sum(r_next**2)
        beta=rtr_next/rtr
        p_next=r_next+beta*p
        
        rtr=rtr_next
        r=r_next
        p=p_next
        x=x_next
        print 'at end of iteration ',i,' residual is ',rtr
    return x

    


plt.ion()
#you can get e.g. HST Deep field images from https://archive.stsci.edu/prepds/xdf/
hdu=fits.open('/Users/sievers/courses/phys_641/hst/hlsp_xdf_hst_acswfc-60mas_hudf_f435w_v1_sci.fits')
image_raw=hdu[0].data
hdu.close()

image_raw=image_raw[1536:3584,1536:3584] #for speed, let's take a smaller patch

noise_est=np.median(np.abs(image_raw-np.median(image_raw)))/0.67 #0.67 is factor from median to std deviation

beam_shapes=[[0.5,5.0],[5.0,0.5]]
nbeam=len(beam_shapes)

maps=[None]*nbeam
beams_ft=[None]*nbeam
beam_var=[None]*nbeam

nx=image_raw.shape[0]
ny=image_raw.shape[1]
xvec=np.arange(image_raw.shape[0])
yvec=np.arange(image_raw.shape[1])

xvec[nx/2:]=xvec[nx/2:]-nx
yvec[ny/2:]=yvec[ny/2:]-ny

xmat=np.repeat([xvec],ny,axis=0).transpose()
ymat=np.repeat([yvec],nx,axis=0)
mapft=np.fft.fft2(image_raw)
for i in range(nbeam):
    rmat=(xmat/beam_shapes[i][0])**2+(ymat/beam_shapes[i][1])**2
    beam=np.exp(-0.5*rmat)
    beam=beam/beam.sum()
    beam_var[i]=np.sum(beam**2)
    beams_ft[i]=np.fft.fft2(beam)
    maps[i]=np.real(np.fft.ifft2(mapft*beams_ft[i]))
    maps[i]=maps[i]+np.random.randn(nx,ny)*noise_est

noises=None

ninv_maps=apply_noise(maps,noises)
rhs=apply_mat_transpose(ninv_maps,beams_ft)
my_map=run_pcg(0*image_raw,rhs,beams_ft,noises,niter=15)

plt.clf();
plt.subplot(2,2,1);plt.imshow(np.sqrt(image_raw),vmax=0.1);plt.colorbar();plt.axis([400,600,500,700])
plt.subplot(2,2,2);plt.imshow(np.sqrt(maps[0]),vmax=0.1);plt.colorbar();plt.axis([400,600,500,700])
plt.subplot(2,2,3);plt.imshow(np.sqrt(maps[1]),vmax=0.1);plt.colorbar();plt.axis([400,600,500,700])
plt.subplot(2,2,4);plt.imshow(np.sqrt(my_map),vmax=0.1);plt.colorbar();plt.axis([400,600,500,700])
plt.savefig('reconstructed_15iter.png')
chisq=0
for i in range(nbeam):
    pred=np.real(np.fft.ifft2(np.fft.fft2(my_map)*beams_ft[i]))
    chisq=chisq+np.sum( (pred-maps[i])**2)
print 'chisq after first iterations is ',chisq

my_map2=run_pcg(my_map,rhs,beams_ft,noises,niter=40)
plt.clf();
plt.subplot(2,2,1);plt.imshow(np.sqrt(image_raw),vmax=0.1);plt.colorbar();plt.axis([400,600,500,700])
plt.subplot(2,2,2);plt.imshow(np.sqrt(maps[0]),vmax=0.1);plt.colorbar();plt.axis([400,600,500,700])
plt.subplot(2,2,3);plt.imshow(np.sqrt(maps[1]),vmax=0.1);plt.colorbar();plt.axis([400,600,500,700])
plt.subplot(2,2,4);plt.imshow(np.sqrt(my_map2),vmax=0.1);plt.colorbar();plt.axis([400,600,500,700])

chisq=0
for i in range(nbeam):
    pred=np.real(np.fft.ifft2(np.fft.fft2(my_map2)*beams_ft[i]))
    chisq=chisq+np.sum( (pred-maps[i])**2)
print 'chisq after second iterations is ',chisq
plt.savefig('reconstructed_40iter.png')

#plt.clf();plt.imshow(np.sqrt(image_raw),vmax=0.1);plt.colorbar();plt.axis([400,600,500,700])
