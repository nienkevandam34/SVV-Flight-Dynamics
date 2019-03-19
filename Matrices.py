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
print("\nAsymmetric eigenvalues: \nHighly damped aperiodic rolling = {}\nDutch roll = {}\nAperiodic spiral motion = {}\n".format(asym_eig[0], asym_eig[1], asym_eig[3]))



# Phugoid: 53:57 - 59:10
# Aperiodic Roll: 59:10 - 60:35
# Short Period: 60:35 - 61:57
# Dutch Roll: 61:57 - 62:47
# Dutch Roll Yaw Damping: 62:47 - 65:20
# Spiral: 65:20 - ...



#step_response=step()
#initial_response=initial()
#impulse_response=impulse()

#[y_asym,t_asym] = step(sysAsym);
#plt.plot(t_asym,y_asym)

name_sym_eigenm = ["Phugoid", "Short Period"]

nr_sym_eigenm = len(name_sym_eigenm)
nr_asym_eigenm = 4

T_st_sym = [(53, 40), (60,25)]
T_en_sym = [(57, 30), (60,31)]

inp_sym = [-0.3259, -1.74]

for i in range(nr_sym_eigenm):
    T = np.arange(0, T_en_sym[i][0]*60 + T_en_sym[i][1] - T_st_sym[i][0]*60 - T_st_sym[i][1] + 0.1, 0.1)
    [y_sym,t_sym,x_sym] = lsim(sysSym, inp_sym[i], T)
    
    #tlistu, delistu = grafiekjesmakenyay.plot_real_data(T_st_sym[i], T_en_sym[i], "Dadc1_tas", data)
    tlista, delista = grafiekjesmakenyay.plot_real_data(T_st_sym[i], T_en_sym[i], "vane_AOA", data)
    tlistt, delistt = grafiekjesmakenyay.plot_real_data(T_st_sym[i], T_en_sym[i], "Ahrs1_Pitch", data)
    tlistq, delistq = grafiekjesmakenyay.plot_real_data(T_st_sym[i], T_en_sym[i], "Ahrs1_bPitchRate", data)
    
    #ax.plot(t, np.rad2deg(y_sym[:,1] + alpha0), label="simulated")
    #ax.set(xlabel="elapsed time [s]")
    #ax.legend()
    
    fig1, ax1 = plt.subplots(2, 2)
    fig1.suptitle("Symmetric: " + name_sym_eigenm[i])
    
    #ax1[0,0].plot(tlistu-tlistu[0], delistu, label="test data")
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



tlistb, delistb = grafiekjesmakenyay.plot_real_data((61,57), (62,10), "Fms1_trueHeading", data)
tlisth, delisth = grafiekjesmakenyay.plot_real_data((61,57), (62,10), "Ahrs1_Roll", data)
tlistp, delistp = grafiekjesmakenyay.plot_real_data((61,57), (62,10), "Ahrs1_bYawRate", data)
tlistr, delistr = grafiekjesmakenyay.plot_real_data((61,57), (62,10), "Ahrs1_bRollRate", data)

T=np.arange(0, len(tlistb)/10, 0.1)
dummy, input_dutchroll = grafiekjesmakenyay.plot_real_data((61,57), (62,10), "delta_r", data)
[y_asym,t_asym,x_sym] = lsim(sysAsym, np.vstack((np.zeros(len(T)), input_dutchroll)).T, T)

figt, axt = plt.subplots(1, 1)
axt.plot(T, input_dutchroll)


fig2, ax2 = plt.subplots(2, 2)
fig2.suptitle("Asymmetric")

#ax2[0,0].plot(tlistb-tlistb[0], delistb, label="test data")
ax2[0,0].plot(tlistb-tlistb[0], y_asym[:,0], label=r"$\beta$", color="orange")
ax2[0,0].set(xlabel="elapsed time [s]", ylabel=r"$\beta$ [°]", title="sideslip angle")
ax2[0,0].legend()

ax2[0,1].plot(tlisth-tlisth[0], delisth, label="test data")
ax2[0,1].plot(tlisth-tlisth[0], y_asym[:,1], label=r"$\phi$")
ax2[0,1].set(xlabel="elapsed time [s]", ylabel=r"$\phi$ [°]", title="roll angle")
ax2[0,1].legend()

ax2[1,0].plot(tlistp-tlistp[0], delistp, label="test data")
ax2[1,0].plot(tlistp-tlistp[0], y_asym[:,2], label="p")
ax2[1,0].set(xlabel="elapsed time [s]", ylabel="p [°/s]", title="yaw rate")
ax2[1,0].legend()

ax2[1,1].plot(tlistr-tlistr[0], delistr, label="test data")
ax2[1,1].plot(tlistr-tlistr[0], y_asym[:,3], label="r")
ax2[1,1].set(xlabel="elapsed time [s]", ylabel="r [°/s]", title="roll rate")
ax2[1,1].legend()

plt.legend()

plt.show()

