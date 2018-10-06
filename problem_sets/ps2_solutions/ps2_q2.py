import numpy as np
def cheb_mat(x,ord):
    mat=np.zeros([len(x),ord+1])
    mat[:,0]=1.0
    mat[:,1]=x
    for i in range(1,ord):
        mat[:,i+1]=2*x*mat[:,i]-mat[:,i-1]
    return mat


x=np.linspace(-1,1,10000)
y=np.exp(x)
ord=100
mat=cheb_mat(x,ord)
lhs=np.dot(mat.transpose(),mat)
rhs=np.dot(mat.transpose(),y)
lhs_inv=np.linalg.inv(lhs)
fitp=np.dot(lhs_inv,rhs)

pred=np.dot(mat,fitp)

#print RMS error - if this number is large, our fit did not do well
fit_rms=np.std(pred-y)
fit_max=np.max(np.abs(pred-y))
print 'for order ',ord,' fit, rms and max errors are  ',fit_rms,fit_max

trunc=6
fitp_trunc=fitp.copy()
fitp_trunc[trunc+1:]=0
max_predicted=np.sum(np.abs(fitp[trunc+1:])) #max error we could have gotten from the coefficients we dropped

pred_trunc=np.dot(mat,fitp_trunc)
trunc_rms=np.std(y-pred_trunc)
trunc_max=np.max(np.abs(y-pred_trunc))
print 'for order ',trunc,' truncation, rms and max errors are ',trunc_rms,trunc_max
print 'for order ',trunc,' truncation, predicted max error is ',max_predicted,' with ratio ',trunc_max/max_predicted #this number sould be no more than 1


mat_trunc=cheb_mat(x,trunc)
lhs=np.dot(mat_trunc.transpose(),mat_trunc)
rhs=np.dot(mat_trunc.transpose(),y)
fitp=np.dot(np.linalg.inv(lhs),rhs)
pred=np.dot(mat_trunc,fitp)

fit_rms=np.std(pred-y)
fit_max=np.max(np.abs(pred-y))
print 'for order ',trunc,' fit, rms and max errors are ',fit_rms,fit_max


