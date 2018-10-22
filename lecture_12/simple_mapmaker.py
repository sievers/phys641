import numpy as np
from scipy import sparse
import pyfits
from matplotlib import pyplot as plt


def read_tod_from_fits(fname,hdu=1):
    f=pyfits.open(fname)
    raw=f[hdu].data
    pixid=raw['PIXID']
    dets=np.unique(pixid)
    ndet=len(dets)
    nsamp=len(pixid)/len(dets)
    print 'nsamp and ndet are ',ndet,nsamp,len(pixid),' on ',fname
    #print raw.names                                                                                                                                                       
    dat={}
    #this bit of odd gymnastics is because a straightforward reshape doesn't seem to leave the data in                                                                     
    #memory-contiguous order, which causes problems down the road                                                                                                          
    #also, float32 is a bit on the edge for pointing, so cast to float64                                                                                                   
    dx=raw['DX']
    #dat['dx']=np.zeros([ndet,nsamp],dtype=type(dx[0]))                                                                                                                 
    dat['dx']=np.zeros([ndet,nsamp],dtype='float64')
    dat['dx'][:]=np.reshape(dx,[ndet,nsamp])[:]
    dy=raw['DY']
    #dat['dy']=np.zeros([ndet,nsamp],dtype=type(dy[0]))                                                                                                                 
    dat['dy']=np.zeros([ndet,nsamp],dtype='float64')
    dat['dy'][:]=np.reshape(dy,[ndet,nsamp])[:]

    tt=np.reshape(raw['TIME'],[ndet,nsamp])
    tt=tt[0,:]
    dt=np.median(np.diff(tt))
    dat['dt']=dt
    pixid=np.reshape(pixid,[ndet,nsamp])
    pixid=pixid[:,0]
    dat['pixid']=pixid
    dat_calib=raw['FNU']
    #dat['dat_calib']=np.zeros([ndet,nsamp],dtype=type(dat_calib[0]))                                                                                                   
    dat['dat_calib']=np.zeros([ndet,nsamp],dtype='float64') #go to double because why not                                                                               
    dat_calib=np.reshape(dat_calib,[ndet,nsamp])
    dat['dat_calib'][:]=dat_calib[:]
    dat['fname']=fname
    f.close()
    return dat

def simple_pcg_squared(x,b,A,Minv=None,niter=20):
    AT=A.transpose()
    r=b-AT*(A*x)
    if not(Minv is None):
        z=r*Minv
    else:
        z=r.copy()
    p=z.copy()
    x=0
    zTr=np.dot(r,z)
    for iter in range(niter):
        Ap=AT*(A*p)
        pAp=np.dot(p,Ap)
        rTr=np.dot(r,r)
        print 'iter is ',iter,' with residual ',rTr
        alpha=zTr/pAp
        x=x+alpha*p
        r_new=r-alpha*Ap
        if not(Minv is None):
            z_new=r_new*Minv
        else:
            z_new=r_new.copy()
        zTr_new=np.dot(r_new,z_new)
        beta=zTr_new/zTr
        p=z_new+beta*p
        r=r_new
        z=z_new
        zTr=zTr_new
    return x
        
def simple_cg_squared(x,b,A,niter=20):
    AT=A.transpose()
    r=b-AT*(A*x)

    p=r.copy()
    x=0
    rTr=np.dot(r,r)
    for iter in range(niter):
        print 'iter is ',iter,' with residual squared ',rTr
        Ap=AT*(A*p)
        pAp=np.dot(p,Ap)
        alpha=np.dot(r,r)/pAp
        x=x+alpha*p
        r_new=r-alpha*Ap
        rTr_new=np.dot(r_new,r_new)
        beta=rTr_new/rTr
        p=r_new+beta*p
        r=r_new
        rTr=rTr_new
    return x


plt.ion()


fname='/Users/sievers/mustang/actcl_0059/data/pointing/Signal_TOD-AGBT18A_175_02-s5.fits'
dat=read_tod_from_fits(fname)

pixsize=2.0
pixsize=pixsize*np.pi/180/3600

xpix=np.asarray( (dat['dx']-dat['dx'].min())/pixsize,dtype='int')
ypix=np.asarray( (dat['dy']-dat['dy'].min())/pixsize,dtype='int')

nxpix=xpix.max()+1
nypix=ypix.max()+1
npix=nxpix*nypix
ipix=xpix+nxpix*ypix


dat_calib=dat['dat_calib'].copy()
ndet=dat_calib.shape[0]
for i in range(ndet):
    dat_calib[i,:]=dat_calib[i,:]-np.median(dat_calib[i,:])

cm=np.median(dat_calib,axis=0)
mat=np.zeros([len(cm),3])
mat[:,0]=1.0
mat[:,1]=np.linspace(-1,1,len(cm))
mat[:,2]=cm

lhs=np.dot(mat.transpose(),mat)
rhs=np.dot(mat.transpose(),dat_calib.transpose())
fitp=np.dot(np.linalg.inv(lhs),rhs)
pred=np.dot(mat,fitp).transpose()
dat_cm=dat_calib-pred




dat_vec=np.reshape(dat_cm,dat_cm.size)
ndat=len(dat_vec)
ipix_vec=np.reshape(ipix,ipix.size)
ind=np.arange(len(ipix_vec))
#A=sparse.csr_matrix((dat_vec,(ipix_vec,ind)),shape=[npix,ndat])
#A=sparse.csr_matrix((dat_vec,(ind,ipix_vec)),shape=[ndat,npix])
A=sparse.csr_matrix((np.ones(len(dat_vec)),(ind,ipix_vec)),shape=[ndat,npix])

map_dirty=A.transpose()*dat_vec
mm=np.reshape(map_dirty,[nypix,nxpix])

map_norm=simple_cg_squared(0*map_dirty,map_dirty,A)
mm_norm=np.reshape(map_norm,[nypix,nxpix])


hits=A.transpose()*np.ones(len(dat_vec))
hits_inv=0*hits
hits_inv[hits>0]=1.0/hits[hits>0]
hh=np.reshape(hits,[nypix,nxpix])

map_norm_pcg=simple_pcg_squared(0*map_dirty,map_dirty,A,hits_inv,niter=2)
mm_pcg=np.reshape(map_norm_pcg,[nypix,nxpix])
