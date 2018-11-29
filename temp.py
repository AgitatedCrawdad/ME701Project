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
        
    plt.plot(heat,mass_flow)
    return heat,mass_flow

def walltemp(z=0.7, L=0.7,K=5,A=0.0005,q_in=1,Tin=69.7,d=0.0254,m_dot = 0.005,Pin=0.101325):
    z = np.linspace(0,0.7,num=100)
    T_wall = []
    T_zplot = []
    for i in range(100):
        sc_liq_in = IAPWS97(P=0.101325, T=Tin+273.15)
        qflux = q_in/(np.pi*d*L)
        #Saturated properties of water
        sat_liq = IAPWS97(P=0.101325, x=0)
        
        #Saturated properties of vapor
        sat_vap = IAPWS97(P=0.101325, x=1)
        
        Hvap = sat_vap.h-sat_liq.h #kJ/kg
        cp = sat_liq.cp # kJ/kgK
#        print(cp)
        Tsat = sat_liq.T-273.15    
        
        M = 18 #Molecular weight
        Pcrit = 220.64 #Critical pressure of water in bars
        
        Pstar = (Pin*10)/Pcrit
        
        T_z = (qflux*np.pi*d*z[i])/(m_dot*cp)+Tin
        
        if T_z > Tsat:
            T_z = Tsat
        T_zplot.append(T_z)
        h_in = sc_liq_in.h
        h_sat = sat_liq.h
        h_gain = (qflux*np.pi*d*z[i])/m_dot;

        
        x_i = (h_in-h_sat)/Hvap
        x_z = (h_gain+h_in-h_sat)/Hvap
        

#        x = x_z - x_i*np.exp((x_z/x_i)-1)
#        print(x)
        x=0
        Re_l = (m_dot*d)/(sat_liq.mu*A)

        Pr_l = sat_liq.Prandt
        
        h_l = 0.023*(sat_liq.k/d)*(Re_l**(0.8))*(Pr_l**0.4)
        
        h_pool = 55 * (Pstar**(0.12))*(qflux**(0.6667))*((-np.log10(Pstar))**(-0.55))*M**(-0.5)
        
        F = (1+x*Pr_l*((sat_liq.rho/sat_vap.rho)-1))**(0.35)
        S = (1+0.055*F**(0.1)*Re_l**(0.16))**(-1)
        
        T = lambda T_w: (qflux-(((F*h_l*(T_w-T_z))**2)+((S*h_pool*(T_w-Tsat))**2))**(0.5))
        if T_z != Tsat:
            T_wall.append(float(optimize.fsolve(T,105)))
        else:
            htp = (((F*h_l)**2)+((S*h_pool)**2))**(0.5)
            T_wall.append((qflux/htp)+Tsat)
            
        
        
        
#        print(T_wall)
    T_wall = np.asarray(T_wall)
    T_zplot = np.asarray(T_zplot)
    plt.plot(T_wall)
    plt.plot(T_zplot)
    plt.show()
#    plt.plot(T_wall-T_zplot)
    
    
if __name__ == "__main__":
#    massflow()
    walltemp()
    
    
    
    
    
    
    
    