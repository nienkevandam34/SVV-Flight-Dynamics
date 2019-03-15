import numpy as np

import matplotlib.pyplot as plt

import read_mat_data


# reading reference data file
data, unit, description, keys = read_mat_data.read_mat('reference_data.mat')

print("Determining stability coefficients ...")
# timing
tb_statf1  = 1157 # sec
te_statf1  = 1920 # sec
timestamps = np.array([1157, 1297, 1426, 1564, 1787, 1920]) # sec, timestamps in excel file

tb_cgsh    = 3062 # sec
te_cgsh    = 3162 # sec


# Constant values concerning atmosphere and gravity
rho0   = 1.2250          # air density at sea level [kg/m^3] 
Tgrad  = -0.0065         # temperature gradient in ISA [K/m]
Temp0  = 288.15          # temperature at sea level in ISA [K]
R      = 287.05          # specific gas constant [m^2/sec^2K]
g      = 9.81            # [m/sec^2] (gravity constant)
gamma  = 1.4


# Aircraft geometry
S      = 30.00	         # wing area [m^2]
Sh     = 0.2 * S         # stabiliser area [m^2]
Sh_S   = Sh / S	         # [ ]
lh     = 0.71 * 5.968    # tail length [m]
c      = 2.0569	         # mean aerodynamic cord [m]
lh_c   = lh / c	         # [ ]
b      = 15.911	         # wing span [m]
bh     = 5.791	         # stabilser span [m]
A      = b ** 2 / S      # wing aspect ratio [ ]
Ah     = bh ** 2 / Sh    # stabilser aspect ratio [ ]
Vh_V   = 1	             # [ ]
ih     = -2 * np.pi / 180   # stabiliser angle of incidence [rad]



# collecting reference data
time      = data["time"]
stat1_ind = np.where(np.logical_and(time>=tb_statf1, time <= te_statf1)) # index for stationary flight 1
t_stat1   = time[stat1_ind]
aoa       = data["vane_AOA"][stat1_ind]                # deg
hp        = data["Dadc1_alt"][stat1_ind] * 0.3048      # m
Vc        = data["Dadc1_cas"][stat1_ind] * 0.5144444   # m/s
Tm        = data["Dadc1_tat"][stat1_ind] + 273.15      # K
le_FU     = data["lh_engine_FU"][stat1_ind] * 0.453592 # kg
re_FU     = data["rh_engine_FU"][stat1_ind] * 0.453592 # kg
delta_e   = data["elevator_dte"][stat1_ind]            # deg

# 6 points only
index_6 = []
for i in range(6):
    index_6.append(np.where(t_stat1 == timestamps[i])[0][0])
index_6 = np.array(index_6)

aoa       = aoa[index_6]     # deg
hp        = hp[index_6]      # m
Vc        = Vc[index_6]      # m/s
Tm        = Tm[index_6]      # K
le_FU     = le_FU[index_6]   # kg
re_FU     = re_FU[index_6]   # kg
delta_e   = delta_e[index_6] # deg



# pressure
p0 = 101325
p  = p0*(1 + Tgrad*hp/Temp0)**(-g/(Tgrad*R))

# mach number, temperature, speed of sound, VTAS, rho, VEAS
M    = np.sqrt(2/(gamma-1) * ( (1 + p0/p*( (1 + (gamma - 1)/(2*gamma) * (rho0/p0) * Vc*Vc)**(gamma/(gamma - 1)) - 1))**((gamma - 1)/gamma) - 1))
T    = Tm/(1 + M*M*(gamma-1)/2)
a    = np.sqrt(gamma*R*T)
VTAS = M*a
rho  = p/(R*T)
VEAS = VTAS*np.sqrt(rho/rho0)

# ramp mass, total fuel used, mass, weight
ramp_mass = 6689.22            # kg
tot_FU    = le_FU + re_FU
m         = ramp_mass - tot_FU
W         = m*g

# CL
CL = 2*W/(rho*S*VTAS*VTAS)



# prepare for thrust calculations
real_alt = data["Dadc1_bcAlt"][stat1_ind][index_6] * 0.3048        # m
T_ISA    = Temp0 + Tgrad*(real_alt)                                # K
delta_T  = T - T_ISA                                               # -

le_ff    = data["lh_engine_FMF"][stat1_ind][index_6] * 0.000125998 # kg/s
re_ff    = data["rh_engine_FMF"][stat1_ind][index_6] * 0.000125998 # kg/s

thrust_calc_in_file = open("matlab.dat", "w")

for i in range(len(timestamps)):
    thrust_calc_in_file.write("{} {} {} {} {}\n".format(float(hp[i]), float(M[i]), float(delta_T[i]), float(le_ff[i]), float(re_ff[i])))

thrust_calc_in_file.close()


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
e = (CL**2)/(CD-CD0)/(np.pi*A)


# =============================================================================
# fig, ax = plt.subplots()
# ax.plot(aoa, CL, marker="+", markersize=5, linestyle="none", label="reference data")
# ax.set(xlabel=r"$\alpha$ [°]", ylabel=r"$C_L$")
# 
# ax.plot(aoa, CLadeg*aoa+CLa0deg, label="least square fit")
# plt.legend()
# 
# print("\nI think that CLa per deg = {}".format(CLadeg))
# print("I think that CLa per rad = {}".format(CLarad))
# 
# fig, ax = plt.subplots()
# ax.plot(aoa, CD, marker="+", markersize=5, linestyle="none", label="reference data")
# ax.set(xlabel=r"$\alpha$ [°]", ylabel=r"$C_D$")
# 
# aa, bb, cc = np.polyfit(aoa, CD, 2)
# ax.plot(aoa, aa*aoa**2 + bb*aoa + cc, label="degree 2 polynomial fit")
# plt.legend()
# 
# # CL vs CD
# fig, ax = plt.subplots()
# ax.plot(CL, CD, marker="+", markersize=5, linestyle="none", label="reference data")
# ax.set(xlabel=r"$C_L$", ylabel=r"$C_D$")
# plt.legend()
# 
# # CL^2 vs CD (should be straight line)
# fig, ax = plt.subplots()
# ax.plot(CL**2, CD, marker="+", markersize=5, linestyle="none", label="reference data")
# ax.set(xlabel=r"$C_L^2$", ylabel=r"$C_D$")
# 
# 
# ax.plot(CL**2, a2*(CL**2) + CD0, label="degree 1 polynomial fit")
# plt.legend()
# 
# print("I think that CD0 = {}".format(CD0))
# print("I think that e = {}".format(e[-2]))
# 
# plt.show()
# =============================================================================




# Citation 550 - Linear simulation

# xcg = 0.25 * c

# Stationary flight condition
hp0    = hp[0]       # pressure altitude in the stationary flight condition [m]
V0     = VTAS[0]     # true airspeed in the stationary flight condition [m/sec]
alpha0 = np.deg2rad(aoa[0])      # angle of attack in the stationary flight condition [rad]
th0    = np.deg2rad(data["Ahrs1_Pitch"][stat1_ind][0]) # pitch angle in the stationary flight condition [rad]

# Aircraft mass
m      = m[0]        # mass [kg]

# aerodynamic properties
e      = e[-2]       # Oswald factor [ ]
CD0    = CD0         # Zero lift drag coefficient [ ]
CLa    = CLarad      # Slope of CL-alpha curve [ ]

# Longitudinal stability
de_cgsh    = np.array([0, -0.5])                # deg  
dde        = de_cgsh[-1] - de_cgsh[0]           # deg
alpha_cgsh = np.array([5.3, 5.3])               # deg
hp_cgsh    = np.array([5730, 57909]) * 0.3048   # m
Vc_cgsh    = np.array([161, 161]) * 0.5144444   # m/s
Tm_cgsh    = np.array([5.0, 5.0]) + 273.15      # K

p0 = 101325
p_cgsh  = p0*(1 + Tgrad*hp_cgsh/Temp0)**(-g/(Tgrad*R))

M_cgsh    = np.sqrt(2/(gamma-1) * ( (1 + p0/p_cgsh*( (1 + (gamma - 1)/(2*gamma) * (rho0/p0) * Vc_cgsh*Vc_cgsh)**(gamma/(gamma - 1)) - 1))**((gamma - 1)/gamma) - 1))
T_cgsh    = Tm_cgsh/(1 + M_cgsh*M_cgsh*(gamma-1)/2)
a_cgsh    = np.sqrt(gamma*R*T_cgsh)
VTAS_cgsh = M_cgsh*a_cgsh
rho_cgsh  = p_cgsh/(R*T_cgsh)
VEAS_cgsh = VTAS_cgsh*np.sqrt(rho_cgsh/rho0)

ramp_mass_cgsh = 6689.22                         # kg
tot_FU_cgsh    = np.array([881, 910]) * 0.453592 # kg
m_cgsh         = ramp_mass_cgsh - tot_FU_cgsh
W_cgsh         = m_cgsh*g

CL_cgsh = 2*W_cgsh/(rho_cgsh*S*VTAS_cgsh*VTAS_cgsh)

dxcg      = -0.05389589                          # thomas' calculation from excel sheet m
CN        = (max(CL_cgsh) + min(CL_cgsh))/2
dde_rad   = np.deg2rad(dde)
Cmde      = -(1/dde_rad)*CN*(dxcg/c)             # elevator effectiveness [ ], App. C: -1.1642

dde_alpha = (1 - 0)/(3.4 - 5.3)                  # ELEVATOR TRIM CURVE!
Cma       = -Cmde*dde_alpha                      # longitudinal stabilty [ ], App. C:  -0.5626

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
CNwa   = CLa                       # Wing normal force slope [ ]
CNha   = 2 * np.pi * Ah / (Ah + 2) # Stabiliser normal force slope [ ]
depsda = 4 / (A + 2)               # Downwash gradient [ ]

# Lift and drag coefficient

CL = 2 * W / (rho * V0 ** 2 * S)                 # Lift coefficient [ ]
CD = CD0 + (CLa * alpha0) ** 2 / (np.pi * A * e) # Drag coefficient [ ]

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
