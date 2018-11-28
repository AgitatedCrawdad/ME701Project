from iapws import IAPWS97
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt

def massflow(L=1.4,K=5,A=0.0005,q_start=0,q_end=5,Tin=99.7):
#q_dot = 1 # kW
#    q_start = 1
#    q_end = 2.7
    g = 9.81 # m/s^2
#    L = 1.4 # m
#    K = 5
#    A = 0.0005 # cross sectional area m^2
    
    #Subcooled inlet properties of water
#    Tin = 99.4 # deg C
    
    sc_liq = IAPWS97(P=0.101325, T=Tin+273.15)
    
    #Saturated properties of water
    sat_liq = IAPWS97(P=0.101325, x=0)
    
    #Saturated properties of vapor
    sat_vap = IAPWS97(P=0.101325, x=1)
    
    Hvap = sat_vap.h-sat_liq.h #kJ/kg
    cp = sc_liq.cp # kJ/kgK
    
    Tsat = sat_liq.T-273.15
    
    #x = lambda m_dot: (q_dot)/(Hvap*m_dot)
    x = lambda m_dot: (q_dot-m_dot*cp*(Tsat-Tin))/(Hvap*m_dot)
    
    left = lambda m_dot: ((sc_liq.rho-sat_liq.rho)+(sat_liq.rho-sat_vap.rho)/(1+(sat_vap.rho/sat_liq.rho)*((1-x(m_dot))/x(m_dot))))*g*L
    
    right = lambda m_dot: (1+x(m_dot)*((sat_liq.rho/sat_vap.rho)-1))*K*(m_dot**2)/(2*sat_liq.rho*A**2)
    
    mass_fun = lambda m_dot: left(m_dot)-right(m_dot)
    mass_flow = []
    heat = []
    heats = np.linspace(q_start,q_end,1000)
    for i in range(0,1000):
        q_dot = heats[i]
        heat.append(q_dot)
        mass_flow.append(float(optimize.fsolve(mass_fun,0.001)))
        
#    plt.plot(heat,mass_flow)
    return heat,mass_flow