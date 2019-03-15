# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from control import ss
from control.matlab import step, impulse, initial, lsim
import matplotlib.pyplot as plt

import numpy as np
from Cit_par_new import *

print("Assemblying system matrices ...")

def sysmat(C1,C2,C3):
    A = -np.linalg.inv(C1)*C2
    B = -np.linalg.inv(C1)*C3
    return A,B



C1_sym=np.matrix([[-2*muc*c/V0, 0, 0, 0],
                 [0, (CZadot-2*muc)*c/V0, 0, 0],
                 [0, 0, -c/V0, 0],
                 [0, Cmadot*c/V0, 0, -2*muc*KY2*c/V0]])

C2_sym=np.matrix([[CXu, CXa, CZ0, CXq],
                 [CZu, CZa, -CX0, CZq+2*muc],
                 [0, 0, 0, 1],
                 [Cmu, Cma, 0, Cmq]])

C3_sym=np.matrix([[CXde*V0],
                  [CZde],
                  [0],
                  [Cmde*V0/c]])

C1_asym = np.matrix([[(CYb-2*mub)*(b/V0),0,0,0],
             [0,-.5*b/V0,0,0],
             [0,0,-2*mub*KX2*(b**2)/(V0**2),2*mub*KXZ*(b**2)/(V0**2)],
             [Cnb*b/V0,0,2*mub*KXZ*(b**2)/(V0**2),-2*mub*KZ2*(b**2)/(V0**2)]])
    
    
C2_asym=np.matrix([[CYb, CL, CYp*b/(2*V0), (CYr-4*mub)*b/(2*V0)],
                    [0, 0, b/(2*V0), 0],
                    [Clb, 0, Clp*b/(2*V0), Clr*b/(2*V0)],
                    [Cnb, 0, Cnp*b/(2*V0), Cnr*b/(2*V0)]])
C3_asym=np.matrix([[CYda, CYdr],
                   [0, 0],
                   [Clda, Cldr],
                   [Cnda, Cndr]])

asm = sysmat(C1_asym,C2_asym,C3_asym)
ssm = sysmat(C1_sym,C2_sym,C3_sym)

A_asym = asm[0]
B_asym = asm[1]

A_sym = ssm[0]
B_sym = ssm[1]

C_asym=np.identity(4)
D_asym=np.zeros((4,2))
C_sym=np.identity(4)
D_sym=np.zeros((4,1))


sysAsym = ss(A_asym,B_asym,C_asym,D_asym)
sysSym=ss(A_sym, B_sym, C_sym, D_sym)



print("Symmetric eigenvalues = {}".format(np.linalg.eig(A_sym)[0]))
print("Asymmetric eigenvalues = {}".format(np.linalg.eig(A_asym)[0]))

#step_response=step()
#initial_response=initial()
#impulse_response=impulse()

#[y_asym,t_asym] = step(sysAsym);
#plt.plot(t_asym,y_asym)

T=np.arange(0, 300, 0.1)
[y_sym,t_sym,x_sym] = lsim(sysSym, np.deg2rad(0.3), T)

fig1, ax1 = plt.subplots(2, 2)
fig1.suptitle("Symmetric")

ax1[0,0].plot(t_sym, y_sym[:,0], label="u")
ax1[0,0].legend()

ax1[0,1].plot(t_sym, np.rad2deg(y_sym[:,1] + alpha0), label=r"$\alpha$")
ax1[0,1].legend()

ax1[1,0].plot(t_sym, np.rad2deg(y_sym[:,2] + th0), label=r"$\theta$")
ax1[1,0].legend()

ax1[1,1].plot(t_sym, y_sym[:,3], label="q")
ax1[1,1].legend()


T=np.arange(0, 2000, 0.5)
[y_sym,t_sym] = step(sysAsym,T)
plt.figure()
plt.title("ASymmetric")
plt.plot(t_sym,y_sym[:,0], label=r"$\beta$")
plt.plot(t_sym,y_sym[:,1], label=r"$\phi$")
plt.plot(t_sym,y_sym[:,2], label="p")
plt.plot(t_sym,y_sym[:,3], label="r")
plt.legend()

plt.show()

