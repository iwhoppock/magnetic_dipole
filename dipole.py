import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits import mplot3d
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


from matplotlib import rc
rc('text', usetex=True)



dt = 1e-2;
mass = 1.;
charge = 1.;
vAc = 3e-4;

duration = 400000;

v = np.array([0., 0.075, 0.05]);
x = np.array([10., 0., 0.]);

E = np.array([0., 0., 0.]);

X = np.zeros((duration,3)); 
V = np.zeros((duration,3)); 

def boris(x,v,charge,mass,vAc,dt,B,E):	
    t = charge / mass * B * 0.5 * dt;
    s = 2. * t / (1. + t*t);
    v_minus = v + charge / (mass * vAc) * E * 0.5 * dt;
    v_prime = v_minus + np.cross(v_minus,t);
    v_plus = v_minus + np.cross(v_prime,s);
    v = v_plus + charge / (mass * vAc) * E * 0.5 * dt;
    x += v * dt;
    return [x,v];

def dipole(x):
	B = np.array([0., 0., 0.]);
	position_mag = np.sqrt(np.dot(x,x));
	m = 1000.;
	B[0] = (3. * m * x[0] * x[2]) / (position_mag**5);
	B[1] = (3. * m * x[1] * x[2]) / (position_mag**5);
	B[2] = ((m) / (position_mag**3)) * (( (3 * x[2] * x[2]) / (position_mag**2) ) - 1.0);
	return B;

for time in range(duration):
    B = dipole(x);
    [x,v] = boris(x,v,charge,mass,vAc,dt,B,E);
    X[time,:] = x;
    V[time,:] = v; 

ax = fig.add_subplot(111, projection='3d')

ax.plot3D(X[:,0],X[:,1],X[:,2],'k',linewidth=2.0)    

ax.set_xlabel(r'$x/d_{\rm p}$',fontsize=16)
ax.set_ylabel(r'$y/d_{\rm p}$',fontsize=16)
ax.set_zlabel(r'$z/d_{\rm p}$',fontsize=16)

plt.show()







