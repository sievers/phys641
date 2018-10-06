import numpy as np
from matplotlib import pyplot as plt

plt.ion()

###Problem 1
#We are solving A^T N^-1 Am=A^T N^-1 d
#with QR, A=QR, so we have
# R^T Q^T N^-1 QRm=R^T Q^T N^-1 d
#R is invertible (is A is non-singular), and Q orthogonal.
#so, so can multiply on the left by R^T^-1, so:
# Q^T N^-1 QRm=Q^T N^-1 d
#note that QQ^T is not equal to I, just Q^TQ, so we can't zap that Q.
#This leaves us with (Q^T N^-1QR)m=(Q^T N^-1d)
#or m=(Q^T N^-1 Q R)^-1(Q^TN^-1 d), or since QR=A, we could also write (Q^T N^-1 A)^-1 (Q^T N^-1 d)
#where the second way is probably a bit more numerically stable than the first.
#also, if we don't have a noise matrix, then we *can* cancel out a Q^TQ to simplify a bit.
#that leaves us with (Q^T Q R)m=Q^Td, or Rm=Q^Td, m=R^-1 Q^T d
def linfit_qr(d,A,N=None):
    Q,R=np.linalg.qr(A)
    if N is None:
        Rinv=np.linalg.inv(R)
        QTd=np.dot(Q.transpose(),d)
        fitp=np.dot(Rinv,QTd)
        return fitp
    else:
        #if we do have a noise, note that if we have N^-1d and N^-1A, we don't need to do anything else
        #it's much faster to make them if N is diagonal, so have a separate code path for diagonal
        #vs. full correlation matrix noise.
        if len(N.shape)==1:  #check for a diagonal noise, which will go faster
            Ninv_d=d/N
            Ninv_A=np.zeros(A.shape)
            for i in range(A.shape[1]):
                Ninv_A[:,i]=A[:,i]/N
        else:
            Ninv=np.linalg.inv(N)
            Ninv_d=np.dot(Ninv,d)
            Ninv_A=np.dot(Ninv,A)
        #once we have Ninv_A and Ninv_d, the two paths look the same
        lhs=np.dot(Q.transpose(),Ninv_A)
        rhs=np.dot(Q.transpose(),Ninv_d)
        lhs_inv=np.linalg.inv(lhs)
        fitp=np.dot(lhs_inv,rhs)
        return fitp


npt=1000

x=np.linspace(-1,1,npt)
y_true=x**2+2.5
y=y_true+np.random.randn(npt)*0.1
ord=30
A=np.zeros([npt,ord+1])
A[:,0]=1.0
for i in range(ord):
    A[:,i+1]=x*A[:,i]

#lets check with an explicit noise matrix set to constant along diagonal to make sure code behaves correctly
N1=np.ones(npt)*3.0
N2=np.eye(npt)*3.0

fitp0=linfit_qr(y,A)
fitp1=linfit_qr(y,A,N1)
fitp2=linfit_qr(y,A,N2)

pred2=np.dot(A,fitp2)
scat2=np.std(y-pred2)

lhs_ref=np.dot(A.transpose(),A)
rhs_ref=np.dot(A.transpose(),y)
fitp_ref=np.dot(np.linalg.inv(lhs_ref),rhs_ref)
pred_ref=np.dot(A,fitp_ref)
scat_ref=np.std(y-pred_ref)

print 'reference fitp is ',fitp_ref
print 'noiseless QR fitp is ',fitp0
print 'diagonal noise QR fitp is ',fitp1
print 'matrix noise QR fitp is ',fitp2

print 'for ',ord,'th order fit, ratio of RMS residuals is ',scat_ref/scat2
