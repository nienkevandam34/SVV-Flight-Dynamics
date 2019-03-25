# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 11:20:18 2019
@author: nienke
"""

import numpy as np
import Cit_par_new

def findeigenvalues(CLarad, CD0, e, data, motion):
    Cmde,Cma= Cit_par_new.stat_meas_2(show_plots=False)
    
    if motion=="Phugoid":
        tstart=61*60+35
        tend=62*60+20
        (hp0, V0, alpha0, th0, m, e, CD0f, CLaf, W, muc, mub, KX2, KZ2, KXZ, KY2, Cmac, CNwa, CNha, depsda, CL, CD, CX0, CXu, CXa, CXadot, CXq, CXde, CZ0, CZu, CZa, CZadot, CZq, CZde, Cmu, Cmadot, Cmq, CYb,  CYbdot, CYp, CYr, CYda, CYdr, Clb, Clp, Clr, Clda, Cldr, Cnb, Cnbdot, Cnp, Cnr, Cnda, Cndr, c, b)=Cit_par_new.stab_coef(tstart, tend, CLarad, CD0, e, data)
        
        # qdot = 0, a = 0
        A_3 = -4*muc**2
        B_3 = 2*muc*CXu
        C_3 = -CZu * CZ0
        coeff_3 = [A_3, B_3, C_3]

        # qdot = 0, adot = 0
        A_4 = 2*muc*(CZa*Cmq - 2*muc*Cma)
        B_4 = 2*muc*(CXu*Cma - Cmu*CXa) + Cmq*(CZu*CXa - CXu*CZa)
        C_4 = CZ0*(Cmu*CZa - CZu*Cma)
        coeff_4 = [A_4, B_4, C_4]
        lambda_ = [np.roots(coeff_3)*V0/c, np.roots(coeff_4)*V0/c]
        
    elif motion=="Short Period":
        tstart=58*60+56
        tend=61*60+4
        (hp0, V0, alpha0, th0, m, e, CD0f, CLaf, W, muc, mub, KX2, KZ2, KXZ, KY2, Cmac, CNwa, CNha, depsda, CL, CD, CX0, CXu, CXa, CXadot, CXq, CXde, CZ0, CZu, CZa, CZadot, CZq, CZde, Cmu, Cmadot, Cmq, CYb,  CYbdot, CYp, CYr, CYda, CYdr, Clb, Clp, Clr, Clda, Cldr, Cnb, Cnbdot, Cnp, Cnr, Cnda, Cndr, c, b)=Cit_par_new.stab_coef(tstart, tend, CLarad, CD0, e, data)
        
        # A*labda^2 + B*labda + C = 0
        
        # constant V
        A_1a = 2*muc*KY2*(2*muc-CZadot)
        B_1a = -2*muc*KY2*CZa - (2*muc+CZq)*Cmadot - (2*muc - CZadot)*Cmq
        C_1a = CZa*Cmq - (2*muc + CZq)*Cma
        coeff_1a = [A_1a, B_1a, C_1a]
        
        # constant V, Cz omitted
        A_1 = 4 * muc**2 * KY2
        B_1 = -2*muc*(KY2*CZa + Cmadot + Cmq)
        C_1 = CZa*Cmq - 2*muc*Cma
        coeff_1 = [A_1, B_1, C_1]
        
        # constant V, gamma constant
        A_2 = -2*muc*KY2
        B_2 = Cmadot + Cmq
        C_2 = Cma
        coeff_2 = [A_2, B_2, C_2]
        lambda_ = [np.roots(coeff_1a)*V0/c, np.roots(coeff_1)*V0/c, np.roots(coeff_2)*V0/c]
        
    elif motion=="Highly Damped Aperiodic Roll":
        tstart=65*60+4
        tend=68*60+0
        (hp0, V0, alpha0, th0, m, e, CD0f, CLaf, W, muc, mub, KX2, KZ2, KXZ, KY2, Cmac, CNwa, CNha, depsda, CL, CD, CX0, CXu, CXa, CXadot, CXq, CXde, CZ0, CZu, CZa, CZadot, CZq, CZde, Cmu, Cmadot, Cmq, CYb,  CYbdot, CYp, CYr, CYda, CYdr, Clb, Clp, Clr, Clda, Cldr, Cnb, Cnbdot, Cnp, Cnr, Cnda, Cndr, c, b)=Cit_par_new.stab_coef(tstart, tend, CLarad, CD0, e, data)
        
        lambda_ = Clp/(4*mub*KX2)*V0/b
        
    elif motion=="Dutch Roll":
        tstart=63*60+3
        tend=63*60+23
        (hp0, V0, alpha0, th0, m, e, CD0f, CLaf, W, muc, mub, KX2, KZ2, KXZ, KY2, Cmac, CNwa, CNha, depsda, CL, CD, CX0, CXu, CXa, CXadot, CXq, CXde, CZ0, CZu, CZa, CZadot, CZq, CZde, Cmu, Cmadot, Cmq, CYb,  CYbdot, CYp, CYr, CYda, CYdr, Clb, Clp, Clr, Clda, Cldr, Cnb, Cnbdot, Cnp, Cnr, Cnda, Cndr, c, b)=Cit_par_new.stab_coef(tstart, tend, CLarad, CD0, e, data)
        
        # pb/2V = 0
        A_6 = 8*mub**2*KZ2
        B_6 = -2*mub*(Cnr + 2*KZ2*CYb)
        C_6 = 4*mub*Cnb + CYb*Cnr
        coeff_6 = [A_6, B_6, C_6]
        
        # pb/2V = 0, only yaw rotation
        A_7 = -2*mub*KZ2
        B_7 = 1/2*Cnr
        C_7 = -Cnb
        coeff_7 = [A_7, B_7, C_7]
        lambda_ = [np.roots(coeff_6)*V0/b, np.roots(coeff_7)*V0/b]
        
    elif motion=="Aperiodic Spiral":
        tstart=68*60+50
        tend=70*60+50
        
        (hp0, V0, alpha0, th0, m, e, CD0f, CLaf, W, muc, mub, KX2, KZ2, KXZ, KY2, Cmac, CNwa, CNha, depsda, CL, CD, CX0, CXu, CXa, CXadot, CXq, CXde, CZ0, CZu, CZa, CZadot, CZq, CZde, Cmu, Cmadot, Cmq, CYb,  CYbdot, CYp, CYr, CYda, CYdr, Clb, Clp, Clr, Clda, Cldr, Cnb, Cnbdot, Cnp, Cnr, Cnda, Cndr, c, b)=Cit_par_new.stab_coef(tstart, tend, CLarad, CD0, e, data)
        
        # All linear and angular accelerations are neglected
        lambda_ = (2*CL*(Clb*Cnr - Cnb*Clr))/(Clp*(CYb*Cnr + 4*mub*Cnb) - Cnp*(CYb*Clr + 4*mub*Clb))*V0/b


    return lambda_



# Symmetric flight






#print("symmetric: ") 
#print("Short period oscillation: ", lambda_1a*V0/c, lambda_1*V0/c, lambda_2*V0/c)
#print("Phugoid: ",lambda_3*V0/c, lambda_4*V0/c)

# Asymmetric flight





#print("asymmetric: ")
#print("Highly damped aperiodic rolling: ",lambda_5*V0/b)
#print("Dutch roll: ", lambda_6*V0/b, lambda_7*V0/b)
#print("Aperiodic spiral motion: ", lambda_8*V0/b)