import numpy
from matplotlib import pyplot as plt

x=numpy.arange(-10,10,0.01)
template=0*x
template[x>0]=numpy.exp(-1.0*x[x>0])

n=1.0


ft1=numpy.fft.rfft(template)
ft2=numpy.fft.rfft(numpy.flipud(template))

ans1=numpy.fft.irfft(ft1*ft1)
ans2=numpy.fft.irfft(ft1*ft2)
print ans1.max()
print ans2.max()

plt.ion()
