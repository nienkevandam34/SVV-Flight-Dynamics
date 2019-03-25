# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from control import ss
from control.matlab import step, impulse, initial, lsim

import matplotlib.pyplot as plt

import numpy as np
import eigenvalues 

import os


import grafiekjesmakenyay

import Cit_par_new

import read_mat_data

# NOTE: the reference data and the real data have a different sign for the 
# aileron (and maybe also for the rudder and elevator??)
use_reference_data = False

if use_reference_data:
    data_name = "reference_data.mat"
else:
    data_name = os.path.join("flight data", "FTISxprt-20190318_124004.mat")

data, unit, description, keys = read_mat_data.read_mat(data_name)



# collect parameters derived from stationary measurements
# first boolean is to indicate whether or not to use the reference data, the
# second one is to print and plot the curves (elevator-trim etc.)
CLarad, CD0, e = Cit_par_new.stat_meas_1(use_reference_data, False)
Cmde, Cma = Cit_par_new.stat_meas_2(use_reference_data, False)

lambda_short = eigenvalues.findeigenvalues(CLarad, CD0, e, data, "Short Period")
lambda_phug = eigenvalues.findeigenvalues(CLarad, CD0, e, data, "Phugoid")
lambda_roll = eigenvalues.findeigenvalues(CLarad, CD0, e, data, "Highly Damped Aperiodic Roll")
lambda_dutch = eigenvalues.findeigenvalues(CLarad, CD0, e, data, "Dutch Roll")
lambda_spiral = eigenvalues.findeigenvalues(CLarad, CD0, e, data, "Aperiodic Spiral")

lambda_list_sym  = [lambda_short, lambda_phug]
lambda_list_asym = [lambda_roll, lambda_dutch, lambda_spiral]

def sysmat(C1,C2,C3):
    A = -np.linalg.inv(C1)*C2
    B = -np.linalg.inv(C1)*C3
    return A,B


def sym_matrices(c, muc, V0, KY2, CX0, CXa, CXde, CXu, CXq, CZ0, CZa, CZadot, 
                 CZde, CZu, CZq, Cmadot, Cmu, Cmq):
    
    print("Assemblying symmetric system matrices ...")
    
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
    
    ssm = sysmat(C1_sym,C2_sym,C3_sym)
    
    A_sym = ssm[0]
    B_sym = ssm[1]
    
    C_sym=np.identity(4)
    D_sym=np.zeros((4,1))
    
    return ss(A_sym, B_sym, C_sym, D_sym), A_sym


def asym_matrices(b, mub, V0, KX2, KZ2, KXZ, CYb, CYbdot, CYp, CYr, CYdr, CYda,
                  CL, Clda, Clb, Clp, Clr, Cldr, Cnda, Cnb, Cnbdot, Cnp, Cnr, 
                  Cndr):
    
    print("Assemblying asymmetric system matrices ...")
    
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
    
    
    A_asym = asm[0]
    B_asym = asm[1]
    
    
    
    C_asym=np.identity(4)
    D_asym=np.zeros((4,2))
    
    return ss(A_asym,B_asym,C_asym,D_asym), A_asym



name_sym_eigenm = ["Short Period", "Phugoid"]

if use_reference_data:
    # reference times
    T_st_sym = [(60,25), (53, 40)]
    T_en_sym = [(60,31), (57, 30)]

else:
    # our times
    T_st_sym = [(61,35), (58, 56)]
    T_en_sym = [(62,20), (61, 4)]

T_st_sym_s = [T_st_sym[0][0]*60 + T_st_sym[0][1], T_st_sym[1][0]*60 + T_st_sym[1][1]]
T_en_sym_s = [T_en_sym[0][0]*60 + T_en_sym[0][1], T_en_sym[1][0]*60 + T_en_sym[1][1]]

inp_sym = [-1.74, -0.3259]

for i in range(len(name_sym_eigenm)):
    
    print("\nSimulating "+name_sym_eigenm[i]+" ...")
    
    # determine length of simulation (add the 0.1 to the end time to include 
    # the end time)
    T = np.arange(0, T_en_sym[i][0]*60 + T_en_sym[i][1] - T_st_sym[i][0]*60 - T_st_sym[i][1] + 0.1, 0.1)
    
    # get elevator input from flight data (only elevator in symmetric 
    # condition) and use this as elevator input for the model
    dummy, elevator_input = grafiekjesmakenyay.plot_real_data(T_st_sym[i], T_en_sym[i], "delta_e", data)
    
    # build system matrices with stability constants
    (hp0, V0, alpha0, th0, m, e, CD0f, CLaf, W, muc, mub, KX2, KZ2, 
     KXZ, KY2, Cmac, CNwa, CNha, depsda, CL, CD, CX0, CXu, CXa, CXadot, 
     CXq, CXde, CZ0, CZu, CZa, CZadot, CZq, CZde, Cmu, Cmadot, Cmq, CYb, 
     CYbdot, CYp, CYr, CYda, CYdr, Clb, Clp, Clr, Clda, Cldr, Cnb, 
     Cnbdot, Cnp, Cnr, Cnda, Cndr, c, b) = Cit_par_new.stab_coef(T_st_sym_s[i], 
                                                                 T_en_sym_s[i], 
                                                                 CLarad, CD0, e,
                                                                 data)
    
    sysSym, A = sym_matrices(c, muc, V0, KY2, CX0, CXa, CXde, CXu, CXq, CZ0, 
                             CZa, CZadot, CZde, CZu, CZq, Cmadot, Cmu, Cmq)
    
    # determine eigenvalues
    eig = np.linalg.eig(A)[0]
    print("Eigenvalue = {}".format(eig[2*i]))
    print("Analytical eigenvalue = {}".format(lambda_list_sym[i][0])) # "percentage off =",lambda_list_sym[i]/eig[2*i]*100)
    print("Difference real = {}".format((np.absolute(np.real(lambda_list_sym[i][0])-np.real(eig[2*i])))/np.real(eig[2*i])*100))
    print("Difference imaginary = {}".format((np.absolute(np.imag(lambda_list_sym[i][0])-np.imag(eig[2*i])))/np.imag(eig[2*i])*100))
    
    # simulate system response
    [y_sym,t_sym,x_sym] = lsim(sysSym, elevator_input, T)#lsim(sysSym, inp_sym[i], T)
    
    # collect the actual response data (from flight test)
    #tlistu, delistu = grafiekjesmakenyay.plot_real_data(T_st_sym[i], T_en_sym[i], "Dadc1_tas", data)
    tlista, delista = grafiekjesmakenyay.plot_real_data(T_st_sym[i], T_en_sym[i], "vane_AOA", data)
    tlistt, delistt = grafiekjesmakenyay.plot_real_data(T_st_sym[i], T_en_sym[i], "Ahrs1_Pitch", data)
    tlistq, delistq = grafiekjesmakenyay.plot_real_data(T_st_sym[i], T_en_sym[i], "Ahrs1_bPitchRate", data)
    
    fig1, ax1 = plt.subplots(2, 2)
    fig1.suptitle("Symmetric: " + name_sym_eigenm[i])
    
    #ax1[0,0].plot(tlistu-tlistu[0], delistu, label="test data")
    ax1[0,0].plot(t_sym, y_sym[:,0], label="u", color="orange")
    ax1[0,0].set(xlabel="elapsed time [s]", ylabel="u [m/s]", title="disturbance velocity")
    ax1[0,0].legend()
    
    ax1[0,1].plot(tlista-tlista[0], delista, label="test data")
    ax1[0,1].plot(tlista-tlista[0], y_sym[:,1] + np.rad2deg(alpha0), label="simulated data$")
    ax1[0,1].set(xlabel="elapsed time [s]", ylabel=r"$\alpha$ [°]", title="angle of attack")
    ax1[0,1].legend()
    
    ax1[1,0].plot(tlistt-tlistt[0], delistt, label="test data")
    ax1[1,0].plot(tlistt-tlistt[0], y_sym[:,2] + np.rad2deg(th0), label="simulated data")
    ax1[1,0].set(xlabel="elapsed time [s]", ylabel=r"$\theta$ [°]", title="pitch angle")
    ax1[1,0].legend()
    
    q0 = delistq[0]
    ax1[1,1].plot(tlistq-tlistq[0], delistq, label="test data")
    ax1[1,1].plot(tlistq-tlistq[0], y_sym[:,3] + q0, label="simulated data")
    ax1[1,1].set(xlabel="elapsed time [s]", ylabel=r"q [°/s]", title="pitch rate")
    ax1[1,1].legend()





name_asym_eigenm = ["Aperiodic Roll", "Dutch Roll", "Spiral"]

if use_reference_data:
    # reference times
    T_st_asym = [(59,10), (61,50), (65,22)]
    T_en_asym = [(59,30), (62,35), (65,50)]

else:
    # our times
    T_st_asym = [(65,4), (63,3), (68,50)]
    T_en_asym = [(68,0), (63,23), (69,20)]

T_st_asym_s = [T_st_asym[0][0]*60 + T_st_asym[0][1], T_st_asym[1][0]*60 + T_st_asym[1][1], T_st_asym[2][0]*60 + T_st_asym[2][1]]
T_en_asym_s = [T_en_asym[0][0]*60 + T_en_asym[0][1], T_en_asym[1][0]*60 + T_en_asym[1][1], T_en_asym[2][0]*60 + T_en_asym[2][1]]

for i in range(len(name_asym_eigenm)):
    
    print("\nSimulating "+name_asym_eigenm[i]+" ...")
    
    # collect the actual response data (from flight test)
    tlistb, delistb = grafiekjesmakenyay.plot_real_data(T_st_asym[i], T_en_asym[i], "Fms1_trueHeading", data)
    tlisth, delisth = grafiekjesmakenyay.plot_real_data(T_st_asym[i], T_en_asym[i], "Ahrs1_Roll", data)
    tlistp, delistp = grafiekjesmakenyay.plot_real_data(T_st_asym[i], T_en_asym[i], "Ahrs1_bYawRate", data)
    tlistr, delistr = grafiekjesmakenyay.plot_real_data(T_st_asym[i], T_en_asym[i], "Ahrs1_bRollRate", data)
    
    # determine length of simulation (add the 0.1 to the end time to include 
    # the end time)
    T = np.arange(0, T_en_asym[i][0]*60 + T_en_asym[i][1] - T_st_asym[i][0]*60 - T_st_asym[i][1] + 0.1, 0.1)
    
    # get aileron and rudder input from flight data (for asymmetric condition) 
    # and use this as aileron and rudder input for the model
    dummy, aileron_input = grafiekjesmakenyay.plot_real_data(T_st_asym[i], T_en_asym[i], "delta_a", data)
    dummy, rudder_input  = grafiekjesmakenyay.plot_real_data(T_st_asym[i], T_en_asym[i], "delta_r", data)
    
    # build system matrices with stability constants
    (hp0, V0, alpha0, th0, m, e, CD0f, CLaf, W, muc, mub, KX2, KZ2, 
     KXZ, KY2, Cmac, CNwa, CNha, depsda, CL, CD, CX0, CXu, CXa, CXadot, 
     CXq, CXde, CZ0, CZu, CZa, CZadot, CZq, CZde, Cmu, Cmadot, Cmq, CYb, 
     CYbdot, CYp, CYr, CYda, CYdr, Clb, Clp, Clr, Clda, Cldr, Cnb, 
     Cnbdot, Cnp, Cnr, Cnda, Cndr, c, b) = Cit_par_new.stab_coef(T_st_asym_s[i], 
                                                                 T_en_asym_s[i],
                                                                 CLarad, CD0, e,
                                                                 data)
    
    sysAsym, A = asym_matrices(b, mub, V0, KX2, KZ2, KXZ, CYb, CYbdot, CYp, 
                               CYr, CYdr, CYda, CL, Clda, Clb, Clp, Clr, Cldr, 
                               Cnda, Cnb, Cnbdot, Cnp, Cnr, Cndr)
    
    # determine eigenvalues
    eig = np.linalg.eig(A)[0]
    
    
    if name_asym_eigenm[i] == "Spiral":
        print("Eigenvalue = {}".format(eig[-1]))
        print("Analytical eigenvalue = {}".format(lambda_list_asym[-1]))
        print("Difference Real = {}".format((np.absolute(np.real(lambda_list_asym[-1])-np.real(eig[-1])))/np.real(eig[-1])*100))
        print("Difference Imaginary = {}".format((np.absolute(np.imag(lambda_list_asym[-1])-np.imag(eig[-1])))/np.imag(eig[-1])*100))
        
    else:
        print("Eigenvalue = {}".format(eig[i]))
        print("Analytical eigenvalue = {}".format(lambda_list_asym[i]))
        print("Difference Real = {}".format((np.absolute(np.real(lambda_list_asym[i])-np.real(eig[i])))/np.real(eig[i])*100))
        print("Difference Imaginary = {}".format((np.absolute(np.imag(lambda_list_asym[i])-np.imag(eig[i])))/np.imag(eig[i])*100))
    
    
    # simulate system response
    [y_asym,t_asym,x_sym] = lsim(sysAsym, np.vstack((-1*aileron_input, -1*rudder_input)).T, T)
    
    
    fig2, ax2 = plt.subplots(2, 2)
    fig2.suptitle("Asymmetric: " + name_asym_eigenm[i])
    
    #ax2[0,0].plot(tlistb-tlistb[0], delistb, label="test data")
    ax2[0,0].plot(tlistb-tlistb[0], y_asym[:,0], label=r"$\beta$", color="orange")
    ax2[0,0].set(xlabel="elapsed time [s]", ylabel=r"$\beta$ [°]", title="sideslip angle")
    ax2[0,0].legend()
    
    phi0 = delisth[0]
    ax2[0,1].plot(tlisth-tlisth[0], delisth, label="test data")
    ax2[0,1].plot(tlisth-tlisth[0], y_asym[:,1] + phi0, label=r"$\phi$")
    ax2[0,1].set(xlabel="elapsed time [s]", ylabel=r"$\phi$ [°]", title="roll angle")
    ax2[0,1].legend()
    
    r0 = delistr[0]
    ax2[1,0].plot(tlistr-tlistr[0], delistr, label="test data")
    ax2[1,0].plot(tlistr-tlistr[0], y_asym[:,2] + r0, label="p")
    ax2[1,0].set(xlabel="elapsed time [s]", ylabel="r [°/s]", title="yaw rate")
    ax2[1,0].legend()
    
    p0 = delistp[0]
    ax2[1,1].plot(tlistp-tlistp[0], delistp, label="test data")
    ax2[1,1].plot(tlistp-tlistp[0], y_asym[:,3] + p0, label="r")
    ax2[1,1].set(xlabel="elapsed time [s]", ylabel="p [°/s]", title="roll rate")
    ax2[1,1].legend()
    
    plt.legend()

plt.show()

