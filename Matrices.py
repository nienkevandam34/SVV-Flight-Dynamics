# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from control import ss
from control.matlab import step, impulse, initial
import matplotlib.pyplot as plt

import numpy as np
from Cit_par import *
def sysmat(C1,C2,C3):
    A = -np.linalg.inv(C1)*C2
    B = -np.linalg.inv(C1)*C3
    return A,B


4
C1_sym=np.matrix([[-2*muc*c/V0, 0, 0, 0],
                 [0, CZadot*c/V0, 0, 0],
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


#step_response=step()
#initial_response=initial()
#impulse_response=impulse()

#[y_asym,t_asym] = step(sysAsym);
#plt.plot(t_asym,y_asym)

T=np.arange(0,10000,0.5)
[y_sym,t_sym] = step(sysSym,T)
y_sym[:,1]=np.array(y_sym[:,1])+np.array(alpha0)

plt.plot(t_sym,y_sym[:,1])

plt.show()

