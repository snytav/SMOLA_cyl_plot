
import numpy as np
from sympy import *
r,th,z,omega,k = symbols('r th z omega,k')

k,m = symbols("k,m")
Sum(k**2, (k, 1, m))

from sympy import *

i, j = symbols("i j", integer=True)

arr = Array([1, 2])

summation(arr[i]*sin(i),(j, 0, i), (i, 0, len(arr)-1))

i = symbols('i', integer=True)
arr = Array([1, 0, 1,0,1])
y = Sum(arr[i]*sin(i),(i,0,4))

y.evalf()

beta = Array([0, -0.575, 0, -0.000799, 0, -0.00000156])

m = symbols("m", integer=True)

Sum(beta[m],(m,1,5))

print(beta[5])

from sympy.plotting import plot
# x,n, p=var('x,n, p')
from sympy.functions import besseli,besselj
#from spiral_field import kc,Bc
import numpy as np

r,z,kc,z0,Bc,m,B0 = symbols('r z kc z0 Bc m B0')
#Ath = Sum(-Bc/kc*beta[m]*besseli(m,kc*r)*sin(kc*(z-z0)),(m,1,5))
# Ath == 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Ath = 0.0

import matplotlib.pyplot as plt
import numpy as np
N = 100
rr = np.linspace(0,0.08,N)
y = np.zeros_like(rr)

#Aths = Ath.subs(z,0).subs(kc,0.349).subs(Bc,0.5).subs(z0,0.0).subs(kc,30*0.349).subs(z,0.0)
# for j,r_n in enumerate(rr):
#     y[j] = Aths.subs(r,r_n)

plt.plot(rr,y)

Az = Sum(beta[m],(m,1,5))
Az

y,x = symbols("y x")
y = x**3 + 5*x
diff(y,x)

Az = Sum(-B0*beta[m]*r*diff(besseli(m,m*k*r),r)*cos(m*(th-k*z)),(m,1,5))  #d/d besseli
k_n = 2*np.pi/0.18

Az

az = Az.subs(th,0.0).subs(r,rr[10]).subs(z,0.0).evalf()
az

az = [Az.subs(B0,0.01).subs(k,k_n).subs(th,0.0).subs(r,xr).subs(z,0.0).evalf() for xr in rr]
az

plt.plot(rr,az)

Az_sub = Az.subs(k,k_n).subs(z,0).subs(B0,0.01)

#plt.plot(rr,Az_sub)

# uniform Ath
Ath = B0*r/2

fi = omega/kc * (r*Ath+ Az/k)
# kc = k*c
# c = 3e10 m/sec
# omega = 1e6 rad/sec
# 0< r < 0.08 m
# z < 2.16 m
# B0 = 100 Gs = 0.01 Tl

Ar = Sum(B0*beta[m]*r*besseli(m,m*k*r)*sin(m*(th-k*z)),(m,1,5))

# NOT phi but fi_s accoding to A.Philippova's presentation p.4
fi_s = Sum(B0*beta[m]*besseli(m,m*k*r)*sin(m*(th-k*z)),(m,1,5))
fi_s = fi_s.subs(B0,0.01).subs(k,k_n)
fi_s

xr = rr[5]
f =fi_s.subs(th,0.01*np.pi).subs(r,xr).subs(z,0.0).evalf()
f

rr[-1]

ff = [fi_s.subs(th,0.01*np.pi).subs(r,xr).subs(z,0.0).evalf() for xr in rr]
print(ff)

#type(ff)
plt.plot(rr,ff)
plt.title('fi_s at z = 0, th = $\pi$/100 ')
plt.xlabel('r, m')
plt.ylabel('fi_s, Volts')

xm = 0.08
N = 3
xx = np.linspace(-xm,xm,N)
yy = np.linspace(-xm,xm,N)
X,Y = np.meshgrid(xx,yy)
Z = np.zeros_like(X)

for i,x in enumerate(xx):
    for j,y in enumerate(yy):
        print(i,j,x,y)
        xr = np.sqrt(x**2+y**2)
        xth = np.arctan2(y,x)
        if xr < xm:
           print(i,j)
           Z[i][j] = fi_s.subs(th,xth).subs(r,xr).subs(z,0.0).evalf()

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
#X, Y = np.meshgrid(x_space, y_space)
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
cmap='viridis',
linewidth=0, antialiased=False)

import matplotlib.pyplot as plt
import numpy as np



# Create the mesh in polar coordinates and compute corresponding Z.
rr = np.linspace(0, 0.1, N)
pp = np.linspace(0, 2*np.pi, N)
R, P = np.meshgrid(rr, pp)
cylZ = np.zeros_like(P)

#p



r,th,z = symbols("r th z")
# Z = ((R**2 - 1)**2)
for i,xr in enumerate(rr):
    for j,xp in enumerate(pp):
          #print(i,j,xr,xp,fi_s)
          t =  fi_s.subs(r,xr).subs(z,0.0).subs(th,xp).evalf()
          qq = 0
          #print(t)
          cylZ[i][j] = t

# Express the mesh in the cartesian system.
X, Y = R*np.cos(P), R*np.sin(P)
from cylindrical_surface_plot import polar_plot
polar_plot(pp,rr,cylZ)
# Z = transfrorm_cylindrical_matrix_to_cartesian(rr,pp,cylZ)
# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# # Plot the surface.
# ax.plot_surface(X, Y, Z, cmap='coolwarm')
#
# # Tweak the limits and add latex math labels.
# #ax.set_zlim(0, 1)
# ax.set_xlabel(r'$\phi_\mathrm{real}$')
# ax.set_ylabel(r'$\phi_\mathrm{im}$')
# ax.set_zlabel(r'$V(\phi)$')
#
# plt.show()



Z

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf=ax.plot_surface(X, Y, Z, rstride=1, cstride=1,
                cmap='coolwarm', edgecolor='none')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z');
plt.colorbar(surf)