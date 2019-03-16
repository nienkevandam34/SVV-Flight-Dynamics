# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from control import ss
from control.matlab import step, impulse, initial, lsim
import matplotlib.pyplot as plt

import grafiekjesmakenyay

import numpy as np
from Cit_par_appC import *

import read_mat_data

data, unit, description, keys = read_mat_data.read_mat("reference_data.mat")

print("Assemblying system matrices ...")

def sysmat(C1,C2,C3):
    A = -np.linalg.inv(C1)*C2
    B = -np.linalg.inv(C1)*C3
    return A,B



C1_sym=np.matrix([[-2*muc*c/(V0**2), 0, 0, 0],
                 [0, (CZadot-2*muc)*c/V0, 0, 0],
                 [0, 0, -c/V0, 0],
                 [0, Cmadot*c/V0, 0, -2*muc*KY2*c/(V0**2)]])

C2_sym=np.matrix([[CXu/V0, CXa, CZ0, CXq*c/V0],
                 [CZu/V0, CZa, -CX0, (CZq+2*muc)*c/V0],
                 [0, 0, 0, c/V0],
                 [Cmu/V0, Cma, 0, Cmq*c/V0]])

C3_sym=np.matrix([[CXde],
                  [CZde],
                  [0],
                  [Cmde]])

C1_asym = np.matrix([[(CYbdot-2*mub)*(b/V0),0,0,0],
             [0,-.5*b/V0,0,0],
             [0,0,-2*mub*KX2*(b**2)/(V0**2),2*mub*KXZ*(b**2)/(V0**2)],
             [Cnbdot*b/V0,0,2*mub*KXZ*(b**2)/(V0**2),-2*mub*KZ2*(b**2)/(V0**2)]])
    
    
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


sym_eig = np.linalg.eig(A_sym)[0]
asym_eig = np.linalg.eig(A_asym)[0]
print("\nSymmetric eigenvalues: \nshort period = {}\nphugoid = {}".format(sym_eig[0], sym_eig[2]))
print("\nAsymmetric eigenvalues: \nHighly damped aperiodic rolling = {}\nDutch roll = {}\nAperiodic spiral motion = {}\n".format(sym_eig[0], sym_eig[1], sym_eig[3]))

#step_response=step()
#initial_response=initial()
#impulse_response=impulse()

#[y_asym,t_asym] = step(sysAsym);
#plt.plot(t_asym,y_asym)

T=np.arange(0, 57*60 + 30 - 54*60 - 17 + 0.1, 0.1)
[y_sym,t_sym,x_sym] = lsim(sysSym, 0.3259, T)

#tlistu, delistu = grafiekjesmakenyay.plot_real_data((54,17), (57,30), "", data)
tlista, delista = grafiekjesmakenyay.plot_real_data((54,17), (57,30), "vane_AOA", data)
tlistt, delistt = grafiekjesmakenyay.plot_real_data((54,17), (57,30), "Ahrs1_Pitch", data)
tlistq, delistq = grafiekjesmakenyay.plot_real_data((54,17), (57,30), "Ahrs1_bPitchRate", data)

#ax.plot(t, np.rad2deg(y_sym[:,1] + alpha0), label="simulated")
#ax.set(xlabel="elapsed time [s]")
#ax.legend()

fig1, ax1 = plt.subplots(2, 2)
fig1.suptitle("Symmetric")

#ax1[0,0].plot(tlistu,delistu, label="test data")
ax1[0,0].plot(t_sym, y_sym[:,0], label="u", color="orange")
ax1[0,0].set(xlabel="elapsed time [s]", ylabel="u [m/s]", title="disturbance velocity")
ax1[0,0].legend()

ax1[0,1].plot(tlista-tlista[0], delista, label="test data")
ax1[0,1].plot(tlista-tlista[0], y_sym[:,1] + np.rad2deg(alpha0), label=r"$\alpha$")
ax1[0,1].set(xlabel="elapsed time [s]", ylabel=r"$\alpha$ [°]", title="angle of attack")
ax1[0,1].legend()

ax1[1,0].plot(tlistt-tlistt[0], delistt, label="test data")
ax1[1,0].plot(tlistt-tlistt[0], y_sym[:,2] + np.rad2deg(th0), label=r"$\theta$")
ax1[1,0].set(xlabel="elapsed time [s]", ylabel=r"$\theta$ [°]", title="pitch angle")
ax1[1,0].legend()

ax1[1,1].plot(tlistq-tlistq[0], delistq, label="test data")
ax1[1,1].plot(tlistq-tlistq[0], y_sym[:,3], label="q")
ax1[1,1].set(xlabel="elapsed time [s]", ylabel=r"q [°/s]", title="pitch rate")
ax1[1,1].legend()


T=np.arange(0, 2000, 0.5)
[y_asym,t_asym,x_sym] = lsim(sysAsym, np.ones((4000,2))*0.3, T)

fig2, ax2 = plt.subplots(2, 2)
plt.title("ASymmetric")
ax2[0,0].plot(t_asym, y_asym[:,0], label=r"$\beta$")
ax2[0,0].legend()

ax2[0,1].plot(t_asym, y_asym[:,1], label=r"$\phi$")
ax2[0,1].legend()

ax2[1,0].plot(t_asym, y_asym[:,2], label="p")
ax2[1,0].legend()

ax2[1,1].plot(t_asym, y_asym[:,3], label="r")
ax2[1,1].legend()

plt.legend()

plt.show()

