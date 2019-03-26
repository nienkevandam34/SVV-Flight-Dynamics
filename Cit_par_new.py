import numpy as np

import matplotlib.pyplot as plt

import subprocess


# specifiy the path to the thrust.exe file, if thrust(1).exe is in the same 
# folder as this Python file, no changes are required. Otherwise, specify the
# full path (and name)
path_to_thrust_file = "thrust(1).exe"



# begin pre-specified values --------------------------------------------------
# Constant values concerning atmosphere and gravity
rho0   = 1.2250          # air density at sea level [kg/m^3] 
Tgrad  = -0.0065         # temperature gradient in ISA [K/m]
Temp0  = 288.15          # temperature at sea level in ISA [K]
R      = 287.05          # specific gas constant [m^2/sec^2K]
g      = 9.81            # [m/sec^2] (gravity constant)
gamma  = 1.4


# Aircraft geometry
S      = 30.00	          # wing area [m^2]
Sh     = 0.2 * S          # stabiliser area [m^2]
Sh_S   = Sh / S	          # [ ]
lh     = 0.71 * 5.968     # tail length [m]
c      = 2.0569	          # mean aerodynamic cord [m]
lh_c   = lh / c	          # [ ]
b      = 15.911	          # wing span [m]
bh     = 5.791	          # stabilser span [m]
A      = b ** 2 / S       # wing aspect ratio [ ]
Ah     = bh ** 2 / Sh     # stabilser aspect ratio [ ]
Vh_V   = 1	              # [ ]
ih     = -2 * np.pi / 180 # stabiliser angle of incidence [rad]
# end pre-specified values ----------------------------------------------------
p0     = 101325           # sea-level pressure [Pa]  
ramp_mass = 6689.22       # kg



def flight_variables(hp, Vc, Tm):
    # use global constants (not a fan of this tough)
    global p0, Tgrad, Temp0, R, gamma, g, rho0
    
    # pressure
    p  = p0*(1 + Tgrad*hp/Temp0)**(-g/(Tgrad*R))
    
    # mach number, temperature, speed of sound, VTAS, rho, VEAS
    M    = np.sqrt(2/(gamma-1) * ( (1 + p0/p*( (1 + (gamma - 1)/(2*gamma) * \
                   (rho0/p0) * Vc*Vc)**(gamma/(gamma - 1)) - 1))**((gamma - 1) \
                                                                  /gamma) - 1))
    T    = Tm/(1 + M*M*(gamma-1)/2)
    a    = np.sqrt(gamma*R*T)
    VTAS = M*a
    rho  = p/(R*T)
    VEAS = VTAS*np.sqrt(rho/rho0)
    
    return p, M, T, a, VTAS, rho, VEAS



def stat_meas_1(use_reference_data=False, show_plots=False):
    
    print("Using stationary measurements 1 ...")
    
    if use_reference_data:
        # Data from REFERENCE_Post_Flight_Datasheet_Flight.xlsx Stationary 
        # measurements CL-CD series 1
        aoa       = np.array([1.7, 2.4, 3.6, 5.4, 8.7, 10.6])               # deg
        hp        = np.array([5010, 5020, 5020, 5030, 5020, 5110]) * 0.3048 # m
        Vc        = np.array([249, 221, 192, 163, 130, 118]) * 0.5144444    # m/s
        Tm        = np.array([12.5, 10.5, 8.8, 7.2, 6, 5.2]) + 273.15       # K
        le_ff     = np.array([798, 673, 561, 463, 443, 474]) * 0.000125998  # kg/s
        re_ff     = np.array([813, 682, 579, 484, 467, 499]) * 0.000125998  # kg/s
        tot_FU    = np.array([360, 412, 447, 478, 532, 570]) * 0.453592     # kg
    
    else:
        # Data from Post_Flight_Datasheet_18_03_V3.xlsx Stationary measurements 
        # CL-CD series 1
        aoa       = np.array([1.65, 2.433333333, 3.416666667, 5.5, 7.6, 10.65])                                   # deg
        hp        = np.array([7998.333333, 9088.333333, 9070, 9080, 8951.666667, 9206.666667]) * 0.3048           # m
        Vc        = np.array([242.8333333, 216.8333333, 190.8333333, 157.6666667, 138.3333333, 118]) * 0.5144444  # m/s
        Tm        = np.array([-0.65, -5.25, -6.5, -8, -8.783333333, -10.16666667]) + 273.15                       # K
        le_ff     = np.array([704, 593.3333333, 481.5, 441, 385.6666667, 423.6666667]) * 0.000125998              # kg/s
        re_ff     = np.array([765, 651.8333333, 538.5, 493.1666667, 414.3333333, 458.3333333]) * 0.000125998      # kg/s
        tot_FU    = np.array([375.5, 617.8333333, 640.1666667, 670.3333333, 704.8333333, 757.6666667]) * 0.453592 # kg
    
    p, M, T, a, VTAS, rho, VEAS = flight_variables(hp, Vc, Tm)
    
    # Reynolds number
    mu = (1.458e-6 * T**(3/2))/(T + 110.4)
    Re = rho*VTAS*c/mu
    
    # ramp mass, total fuel used, mass, weight
    m         = ramp_mass - tot_FU
    W         = m*g
    
    # CL
    CL = 2*W/(rho*S*VTAS*VTAS)
    
    
    # prepare for thrust calculations
    T_ISA    = Temp0 + Tgrad*(hp)                                       # K
    delta_T  = T - T_ISA                                                # -
    
    thrust_calc_in_file = open("matlab.dat", "w")
    
    for i in range(len(delta_T)):
        thrust_calc_in_file.write("{} {} {} {} {}\n".format(float(hp[i]), float(M[i]), float(delta_T[i]), float(le_ff[i]), float(re_ff[i])))
    
    thrust_calc_in_file.close()
    
    # call thrust.exe to calculate the thrust
    subprocess.call([path_to_thrust_file])
    
    # process calculated thrust
    thrust_le = []
    thrust_re = []
    
    with open("thrust.dat") as f:
        all_lines = f.readlines()
    f.close()
    
    for line in all_lines:
        thrust_le.append(float(line.strip().split()[0]))
        thrust_re.append(float(line.strip().split()[1]))
    
    thrust_le = np.array(thrust_le)
    thrust_re = np.array(thrust_re)
    
    thrust_imp = thrust_le + thrust_re
    
    # CD
    CD = 2*thrust_imp/(rho*S*VTAS*VTAS)
    
    
    # Results
    CLadeg, CLa0deg = np.polyfit(aoa, CL, 1) # lift curve gradient and CL at alpha = 0 per deg
    CLarad, CLa0rad = np.polyfit(np.deg2rad(aoa), CL, 1) # lift curve gradient and CL at alpha = 0 per rad
    a2, CD0 = np.polyfit(CL**2, CD, 1)
    e = (sum(CL**2)/len(CL**2))/(((sum(CD)/len(CD))-CD0)*(np.pi*A)) # use avg CL and CD
    
    if show_plots:        
        print("\nCLa per deg = {}".format(CLadeg))
        print("CLa per rad = {}".format(CLarad))
        
        print("CD0 = {}".format(CD0))
        print("e = {}".format(e))
        
        fig, ax = plt.subplots(1, 2)
        
        fig.suptitle(r"Clean configuration, gear and HLDs retracted, $M$ from {:0.2} to {:0.2} and $R_e$ from {:0.2e} to {:0.2e}".format(min(M), max(M), min(Re), max(Re)))
        
        ax[0].plot(aoa, CL, marker="o", linestyle="none", label="reference data")
        ax[0].plot(aoa, CLadeg*aoa+CLa0deg, label="least square fit")
        ax[0].set(xlabel=r"$\alpha$ [°]", ylabel=r"$C_L$", title=r"$C_L$ vs $\alpha$")
        ax[0].legend()
        
        ax[1].plot(CL, CD, marker="o", linestyle="none", label="reference data")
        ap, bp, cp = np.polyfit(CL, CD, 2)
        ax[1].plot(CL, ap*CL**2 + bp*CL + cp, label="parabolic fit")
        ax[1].set(xlabel=r"$C_L$", ylabel=r"$C_D$", title=r"$C_L$ vs $C_D$")
        ax[1].legend()
        
        plt.show()
    
    return CLarad, CD0, e



# Citation 550 - Linear simulation

# xcg = 0.25 * c

def stat_meas_2(use_reference_data=False, show_plots=False):
    
    print("Using stationary measurements 2 ...")
    
    if use_reference_data:
        # Data from REFERENCE_Post_Flight_Datasheet_Flight.xlsx 
        # shift in centre of gravity
        de_cgsh     = np.array([0, -0.5])              # deg  
        dde         = de_cgsh[-1] - de_cgsh[0]         # deg
        alpha_cgsh  = np.array([5.3, 5.3])             # deg
        hp_cgsh     = np.array([5730, 57909]) * 0.3048 # m
        Vc_cgsh     = np.array([161, 161]) * 0.5144444 # m/s
        Tm_cgsh     = np.array([5.0, 5.0]) + 273.15    # K
        tot_FU_cgsh = np.array([881, 910]) * 0.453592  # kg
        
        dxcg        = -0.05389589                      # thomas' calculation
        
        # stationary measurement Elevator Trim Curve
        alpha_et   = np.array([5.3, 6.3, 7.3, 8.5, 4.5, 4.1, 3.4])                 # deg
        de_et      = np.array([0, -0.4, -0.9, -1.5, 0.4, 0.6, 1])                  # deg
        hp_et      = np.array([6060, 6350, 6550, 6880, 6160, 5810, 5310]) * 0.3048 # m
        Vc_et      = np.array([161, 150, 140, 130, 173, 179, 192]) * 0.5144444     # m/s
        Tm_et      = np.array([5.5, 4.5, 3.5, 2.5, 5.0, 6.2, 8.2]) + 273.15        # K
        le_ff_et   = np.array([462, 458, 454, 449, 465, 472, 482]) * 0.000125998   # kg/s
        re_ff_et   = np.array([486, 482, 477, 473, 489, 496, 505]) * 0.000125998   # kg/s
        tot_FU_et  = np.array([664, 694, 730, 755, 798, 825, 846]) * 0.453592      # kg
        F_stick_et = np.array([0, -23, -29, -46, 26, 40, 83])                      # N
    
    else:
        # Data from Post_Flight_Datasheet_18_03_V3.xlsx 
        # shift in centre of gravity
        de_cgsh     = np.array([0.2, -0.25])                           # deg  
        dde         = de_cgsh[-1] - de_cgsh[0]                         # deg
        alpha_cgsh  = np.array([5.05, 5.116666667])                    # deg
        hp_cgsh     = np.array([11063.33333, 11043.33333]) * 0.3048    # m
        Vc_cgsh     = np.array([161.1666667, 161.6666667]) * 0.5144444 # m/s
        Tm_cgsh     = np.array([-12.35, -12.2]) + 273.15               # K
        tot_FU_cgsh = np.array([1037, 1058.5]) * 0.453592              # kg
        
        dxcg        = -0.05516095                                      # thomas' calculation
        
        # stationary measurement Elevator Trim Curve
        alpha_et   = np.array([5.266666667, 6.333333333, 7.3, 8.616666667, 4.35, 3.733333333, 3.35])                          # deg
        de_et      = np.array([0.3, -0.2, -0.6, -1.25, 0.6, 0.8, 1.1])                                                        # deg
        hp_et      = np.array([12088.33333, 12306.66667, 12523.33333, 12808.33333, 11840, 11473.33333, 11033.33333]) * 0.3048 # m
        Vc_et      = np.array([160.3333333, 148.8333333, 139, 130.5, 173, 180.5, 189.6666667]) * 0.5144444                    # m/s
        Tm_et      = np.array([-14.2, -15.2, -16.2, -17.2, -13, -12, -10.5]) + 273.15                                         # K
        le_ff_et   = np.array([402.5, 399, 396, 390, 408.5, 414, 420]) * 0.000125998                                          # kg/s
        re_ff_et   = np.array([447.5, 444.5, 441.5, 435, 456, 462, 468]) * 0.000125998                                        # kg/s
        tot_FU_et  = np.array([894.5, 909.5, 926.5, 945, 963.5, 978, 995]) * 0.453592                                         # kg
        F_stick_et = np.array([1, -21.5, -30.5, -40.5, 29.5, 47.5, 65])                                                       # N
    
    
    # cg-shift calculations for Cmde
    (p_cgsh, M_cgsh, T_cgsh, a_cgsh, 
     VTAS_cgsh, rho_cgsh, VEAS_cgsh) = flight_variables(hp_cgsh, Vc_cgsh, Tm_cgsh)
    
    m_cgsh      = ramp_mass - tot_FU_cgsh
    W_cgsh      = m_cgsh*g
    
    CL_cgsh   = 2*W_cgsh/(rho_cgsh*S*VTAS_cgsh*VTAS_cgsh)
    CN        = (max(CL_cgsh) + min(CL_cgsh))/2
    dde_rad   = np.deg2rad(dde)
    Cmde      = -(1/dde_rad)*CN*(dxcg/c)        # elevator effectiveness [ ]
    
    # elevator trim curve calculations for Cma
    (p_et, M_et, T_et, a_et, 
     VTAS_et, rho_et, VEAS_et) = flight_variables(hp_et, Vc_et, Tm_et)
    
    m_et      = ramp_mass - tot_FU_et
    W_et      = m_et*g
    
    # Reduced equivalent airspeed
    Ws = 60500                             # N
    V_et_red = VEAS_et * np.sqrt(Ws/W_et)  # m/s
    
    dde_alpha, intersect = np.polyfit(alpha_et, de_et, 1)
    
    Cma       = -Cmde*dde_alpha
    
    # determine real thrust coefficient
    T_ISA_et    = Temp0 + Tgrad*(hp_et)                                       # K
    delta_T_et  = T_et - T_ISA_et                                               # -
    
    thrust_calc_in_file = open("matlab.dat", "w")
    for i in range(len(delta_T_et)):
        thrust_calc_in_file.write("{} {} {} {} {}\n".format(float(hp_et[i]), float(M_et[i]), float(delta_T_et[i]), float(le_ff_et[i]), float(re_ff_et[i])))
    thrust_calc_in_file.close()
    
    subprocess.call([path_to_thrust_file])
    
    thrust_le = []
    thrust_re = []
    
    with open("thrust.dat") as f:
        all_lines = f.readlines()
    f.close()
    
    for line in all_lines:
        thrust_le.append(float(line.strip().split()[0]))
        thrust_re.append(float(line.strip().split()[1]))
    
    thrust_le = np.array(thrust_le)
    thrust_re = np.array(thrust_re)
    
    thrust_tot = thrust_le + thrust_re
    
    TC = thrust_tot/(0.5*rho_et*VTAS_et**2*S)
    
    # Get standard thrust coefficient by entering the standard fuel flow into 
    # the thrust.exe file
    thrust_calc_in_file = open("matlab.dat", "w")
    thrust_calc_in_file.write("{} {} {} {} {}\n".format(float(np.average(hp_et)), float(np.average(M_et)), float(np.average(delta_T_et)), 0.048, 0.048))
    thrust_calc_in_file.close()
    
    subprocess.call([path_to_thrust_file])
    
    thrust_le = []
    thrust_re = []
    
    with open("thrust.dat") as f:
        all_lines = f.readlines()
    f.close()
    
    thrust_le = float(line.strip().split()[0])
    thrust_re = float(line.strip().split()[1])
    
    thrust_tots = thrust_le + thrust_re
    
    TCs = thrust_tots/(0.5*np.average(rho_et)*np.average(VTAS_et)**2*S)
    
    # reduced elevator deflection
    Cmtc        = -0.0064 # appendix C
    d_e_eq_meas = 0.5
    
    de_et_red   = d_e_eq_meas - 1/Cmde * Cmtc * (TCs-TC)
    
    # reduced stick force curve
    # Force for friction is neglected (see reader and appendix B of assignment)
    F_stick_et_red = F_stick_et * Ws/W
    
    
    if show_plots:
        print("Cmde = {}".format(Cmde))
        print("Cma = {}".format(Cma))
        
        fig, ax = plt.subplots(2, 2)
        
        fig.suptitle("Elevator trim curves")
        
        ax[0,0].plot(alpha_et, de_et, marker="o", linestyle="None", label="Data")
        ax[0,0].plot(alpha_et, dde_alpha*alpha_et + intersect, label="Linear Fit")
        ax[0,0].set(xlabel=r"$\alpha$ [°]", ylabel=r"$\delta_e$ [°]", title="Real deflection vs real angle of attack")
        ax[0,0].legend()
        
        ax[0,1].plot(Vc_et, de_et, marker="o", linestyle="None")
        ax[0,1].set(xlabel=r"$V_c$ [m/s]", ylabel=r"$\delta_e$ [°]", title="Real deflection vs calibrated airspeed")
        
        ax[1,0].plot(V_et_red, de_et, marker="o", linestyle="None")
        ax[1,0].set(xlabel=r"$\tilde{V}_e$ [m/s]", ylabel=r"$\delta_e$ [°]", title="Real deflection vs reduced equivalent airspeed")
        
        ax[1,1].plot(V_et_red, de_et_red, marker="o", linestyle="None")
        ax[1,1].set(xlabel=r"$\tilde{V}_e$ [m/s]", ylabel=r"$\tilde{\delta}_e$ [°]", title="Reduced Elevator-Trim Curve")
        
        
        fig2, ax2 = plt.subplots(2, 2)
        
        fig2.suptitle("Stick force curves")
        
        ax2[0,0].plot(alpha_et, F_stick_et, marker="o", linestyle="None")
        ax2[0,0].set(xlabel=r"$\alpha$ [°]", ylabel="F [N]", title="Real force vs real angle of attack")
        
        ax2[0,1].plot(Vc_et, F_stick_et, marker="o", linestyle="None")
        ax2[0,1].set(xlabel=r"$V_c$ [m/s]", ylabel="F [N]", title="Real force vs calibrated airspeed")
        
        ax2[1,0].plot(V_et_red, F_stick_et, marker="o", linestyle="None")
        ax2[1,0].set(xlabel=r"$\tilde{V}_e$ [m/s]", ylabel="F [N]", title="Real force vs reduced equivalent airspeed")
        
        ax2[1,1].plot(V_et_red, F_stick_et_red, marker="o", linestyle="None")
        ax2[1,1].set(xlabel=r"$\tilde{V}_e$ [m/s]", ylabel="F [N]", title="Reduced Elevator Control Force Curve")
        
    return Cmde, Cma


def stab_coef(tstart, tend, CLarad, CD0, e, data):
    
    print("Determining stability coefficients ...")
    
    # collecting data
    time  = data["time"]
    ind   = np.where(np.logical_and(time>=tstart, time <= tend)) # index
    hp    = data["Dadc1_alt"][ind] * 0.3048                      # m
    le_FU = data["lh_engine_FU"][ind] * 0.453592                 # kg
    re_FU = data["rh_engine_FU"][ind] * 0.453592                 # kg
    VTAS  = data["Dadc1_tas"][ind] * 0.5144444                   # m/s
    
    # total fuel used, mass, weight
    tot_FU    = le_FU + re_FU
    m         = ramp_mass - tot_FU
    
    hp0    = hp[0]          # pressure altitude at the start of the maneuvre [m]
    V0     = VTAS[0]        # true airspeed at the start of the maneuvre [m/sec]
    alpha0 = np.deg2rad(data["vane_AOA"][ind][0]) # angle of attack at the start of the maneuvre [rad]
    th0    = np.deg2rad(data["Ahrs1_Pitch"][ind][0]) # pitch angle at the start of the maneuvre [rad]
    
    # Aircraft mass
    m      = m[0]        # mass [kg]
    
    # aerodynamic properties
    CD0f   = CD0         # Zero lift drag coefficient [ ]
    CLaf   = CLarad      # Slope of CL-alpha curve [ ]
    ## End Initial Conditions ---
    
    
    
    # begin pre-specified values --------------------------------------------------
    # air density [kg/m^3]  
    rho    = rho0 * ((1+(Tgrad * hp0 / Temp0))) ** (-((g / (Tgrad*R)) + 1))   
    W      = m * g            # [N]       (aircraft weight)
    
    # Constant values concerning aircraft inertia
    
    muc    = m / (rho * S * c)
    mub    = m / (rho * S * b)
    KX2    = 0.019
    KZ2    = 0.042
    KXZ    = 0.002
    KY2    = 1.25 * 1.114
    
    # Aerodynamic constants
    
    Cmac   = 0                         # Moment coefficient about the aerodynamic centre [ ]
    CNwa   = CLaf                      # Wing normal force slope [ ]
    CNha   = 2 * np.pi * Ah / (Ah + 2) # Stabiliser normal force slope [ ]
    depsda = 4 / (A + 2)               # Downwash gradient [ ]
    
    # Lift and drag coefficient
    
    CL = 2 * W / (rho * V0 ** 2 * S)                  # Lift coefficient [ ]
    CD = CD0 + (CLaf * alpha0) ** 2 / (np.pi * A * e) # Drag coefficient [ ]
    
    # Stabiblity derivatives
    
    CX0    = W * np.sin(th0) / (0.5 * rho * V0 ** 2 * S)
    CXu    = -0.02792
    CXa    = -0.47966
    CXadot = +0.08330
    CXq    = -0.28170
    CXde   = -0.03728
    
    CZ0    = -W * np.cos(th0) / (0.5 * rho * V0 ** 2 * S)
    CZu    = -0.37616
    CZa    = -5.74340
    CZadot = -0.00350
    CZq    = -5.66290
    CZde   = -0.69612
    
    Cmu    = +0.06990
    Cmadot = +0.17800
    Cmq    = -8.79415
    
    CYb    = -0.7500
    CYbdot =  0     
    CYp    = -0.0304
    CYr    = +0.8495
    CYda   = -0.0400
    CYdr   = +0.2300
    
    Clb    = -0.10260
    Clp    = -0.71085
    Clr    = +0.23760
    Clda   = -0.23088
    Cldr   = +0.03440
    
    Cnb    =  +0.1348
    Cnbdot =   0     
    Cnp    =  -0.0602
    Cnr    =  -0.2061
    Cnda   =  -0.0120
    Cndr   =  -0.0939
    # end pre-specified values ----------------------------------------------------
    
    return (hp0, V0, alpha0, th0, m, e, CD0f, CLaf, W, muc, mub, KX2, KZ2, 
            KXZ, KY2, Cmac, CNwa, CNha, depsda, CL, CD, CX0, CXu, CXa, CXadot, 
            CXq, CXde, CZ0, CZu, CZa, CZadot, CZq, CZde, Cmu, Cmadot, Cmq, CYb, 
            CYbdot, CYp, CYr, CYda, CYdr, Clb, Clp, Clr, Clda, Cldr, Cnb, 
            Cnbdot, Cnp, Cnr, Cnda, Cndr, c, b)