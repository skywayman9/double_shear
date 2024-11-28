# from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
# from mpl_toolkits.mplot3d import Axes3D
# import scipy.interpolate as interpolate


# # For latex font, i guess so 
plt.rc('text', usetex=True)
plt.rc('font', family='arial')

plt.rcParams.update({'font.size': 12})

x,y,d,u,v,vort=np.loadtxt('Rslt0007.plt', delimiter=None, unpack=True,skiprows=3)


g=160
k=160


x = x.reshape(g,k)
y = y.reshape(g,k)
d = d.reshape(g,k)
vort = vort.reshape(g,k)
u = u.reshape(g,k)
v = v.reshape(g,k)


print(v.min())
print(v.max())

plt.contour(x,y,abs(vort),40,linewidths=0.5,colors=('r'))


plt.ylabel(r'\textbf{y}')
plt.xlabel(r'\textbf{x}')
fig1 = plt.gcf()
# plt.colorbar()
fig1.set_size_inches(w=5,h=5, forward=True)

# fig1.savefig('linear56-dsl.png', dpi=600,bbox_inches='tight', pad_inches = 0.1)
plt.show()