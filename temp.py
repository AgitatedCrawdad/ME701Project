from iapws import IAPWS97
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt

def massflow3(height=1.4,L=0.7,K=5,A=0.0005,q_start=0,q_end=5,Tin=99.7):
#    q_dot = 1 # kW
#    q_start = 1
    q_end = 5
    g = 9.81 # m/s^2
#    L = 1.4 # m
#    K = 5
#    A = 0.0005 # cross sectional area m^2
    D = 0.0254
    #Subcooled inlet properties of water
#    Tin = 99.7 # deg C
#    m_dot = 0.0010
    
    sc_liq = IAPWS97(P=0.101325, T=Tin+273.15)
    
    #Saturated properties of water
    sat_liq = IAPWS97(P=0.101325, x=0)
    
    #Saturated properties of vapor
    sat_vap = IAPWS97(P=0.101325, x=1)
    
    Hvap = sat_vap.h-sat_liq.h #kJ/kg


    h_in = sc_liq.h
    h_sat = sat_liq.h
    
    cp = sc_liq.cp # kJ/kgK

    
    Tsat = sat_liq.T-273.15

#    h_gain = ((q_dot)/m_dot)
    
    x_i = (h_in-h_sat)/Hvap
    
#    x_z = (h_gain+h_in-h_sat)/Hvap
    
#    x_range = np.linspace(x_i,x_z,100)
#    x = x_range - x_i*np.exp((x_range/x_i)-1) 
    
#    x = lambda m_dot: (q_dot-m_dot*cp*(Tsat-Tin))/(Hvap*m_dot)
#    x = lambda m_dot: (q_dot)/(Hvap*m_dot)
#    x = lambda m_dot: ((q_dot/m_dot)+h_in-h_sat)/Hvap
    x = lambda m_dot: (((q_dot)/m_dot)+h_in-h_sat)/Hvap - x_i*np.exp((((((q_dot)/m_dot)+h_in-h_sat)/Hvap)/x_i)-1)  
    
    Re_l = lambda m_dot: (((m_dot*D)/(A*sat_liq.mu)))
    Re_v = lambda m_dot: (((m_dot*D)/(A*sat_vap.mu)))
    
    f_l = lambda m_dot: 0.079/(Re_l(m_dot)**0.25)
    f_v = lambda m_dot: 0.079/(Re_v(m_dot)**0.25)
    
    E = lambda m_dot: ((1-x(m_dot))**2)+(x(m_dot)**2)*(sat_liq.rho*f_v(m_dot))/(sat_vap.rho*f_l(m_dot))
    
    F = lambda m_dot: (x(m_dot)**0.78)*((1-x(m_dot))**0.224)
    
    H = ((sat_liq.rho/sat_vap.rho)**(0.91))*((sat_vap.mu/sat_liq.mu)**0.19)*((1-(sat_vap.mu/sat_liq.mu))**0.7)
    
    Fr = lambda m_dot: (m_dot**2)/(g*D*(A**2)*(rhoH(m_dot)**2))
    
    We = lambda m_dot: ((m_dot**2)*D)/((rhoH(m_dot))*sat_liq.sigma)
    
    tpm = lambda m_dot: E(m_dot)+((3.24*F(m_dot)*H/((Fr(m_dot)**0.045)*(We(m_dot)**0.035))))
    
    alpha = lambda m_dot: 1/(1+(((1-x(m_dot))/x(m_dot))*(sat_vap.rho/sat_liq.rho)))
    
    rhoH = lambda m_dot: (sat_liq.rho*(1-alpha(m_dot))+sat_vap.rho*alpha(m_dot))
#    rhoH = lambda m_dot: ((x(m_dot)/sat_vap.rho)+((1-x(m_dot))/sat_liq.rho))**(-1)
    
    dpdz_l = lambda m_dot: K*(2*f_l(m_dot)*(m_dot**2))/(D*(A**2)*sat_liq.rho)
#    left = lambda m_dot: tpm(m_dot)*K*2*(0.079/(Re_l(m_dot)**0.25))*(L/D)*((m_dot**2)/(A**2))*(1/rhoH(m_dot))
#    left = lambda m_dot: K*dpdz_l(m_dot)*(height)*tpm(m_dot)
    
    left = lambda m_dot: dpdz_l(m_dot)*(2*L*tpm(m_dot)+(height-L)*tpm(m_dot)+(height-L))
#    left = lambda m_dot: (1+x(m_dot)*((sat_liq.rho/sat_vap.rho)-1))*K*(m_dot**2)/(2*sat_liq.rho*A**2)
    
#    right = lambda m_dot: ((sc_liq.rho-sat_liq.rho)+(sat_liq.rho-sat_vap.rho)/(1+(sat_vap.rho/sat_liq.rho)*((1-x(m_dot))/x(m_dot))))*g*L
    right = lambda m_dot: g*(sat_liq.rho-sat_vap.rho)*((height)*alpha(m_dot))
    
    
    mass_fun = lambda m_dot: left(m_dot)-right(m_dot)
    mass_flow = []
    heat = []
    heats = np.linspace(q_start,q_end,100)
    for i in range(0,100):
        q_dot = heats[i]
        heat.append(q_dot)
        
        mass_flow.append(float(optimize.fsolve(mass_fun,0.01)))
        
        
#        print(left(mass_flow[i]),right(mass_flow[i]))
        
    plt.plot(heat,mass_flow)
#    print(q_dot,mass_flow[-1],h_in,h_sat,Hvap)
#    print(x(mass_flow))
#    print(alpha(mass_flow))
#    print(tpm(mass_flow[-1]))
#    print(right(mass_flow))
#    print(mass_fun(mass_flow[0]))
    return heat,mass_flow


def massflow2(L=1.4,K=5,A=0.0005,q_start=0,q_end=5,Tin=99.7):
#    q_dot = 1 # kW

    g = 9.81 # m/s^2
#    L = 1.4 # m
#    K = 5
#    A = 0.0005 # cross sectional area m^2
    D = 0.0254
    #Subcooled inlet properties of water
#    Tin = 99.4 # deg C
#    m_dot = 0.0010
    
    sc_liq = IAPWS97(P=0.101325, T=Tin+273.15)
    
    #Saturated properties of water
    sat_liq = IAPWS97(P=0.101325, x=0)
    
    #Saturated properties of vapor
    sat_vap = IAPWS97(P=0.101325, x=1)
    
    Hvap = sat_vap.h-sat_liq.h #kJ/kg
    cp = sc_liq.cp # kJ/kgK

    
    Tsat = sat_liq.T-273.15
    h_in = sc_liq.h
    h_sat = sat_liq.h
#    h_gain = ((q_dot)/m_dot)
    
    x_i = (h_in-h_sat)/Hvap
    
#    x_z = (h_gain+h_in-h_sat)/Hvap
        
#    x = x_z - x_i*np.exp((x_z/x_i)-1) 
    
    #x = lambda m_dot: (q_dot)/(Hvap*m_dot)
    x = lambda m_dot: (((q_dot)/m_dot)+h_in-h_sat)/Hvap - x_i*np.exp((((((q_dot)/m_dot)+h_in-h_sat)/Hvap)/x_i)-1)  
#    x = lambda m_dot: (q_dot-m_dot*cp*(Tsat-Tin))/(Hvap*m_dot)

    
    alpha = lambda m_dot: 1/(1+(((1-x(m_dot))/x(m_dot))*(sat_vap.rho/sat_liq.rho)))
    
    Re = lambda m_dot: (((m_dot*D)/(A*(x(m_dot)*sat_vap.mu+((1-x(m_dot))*sat_liq.mu)))))
    
    rhoH = lambda m_dot: (sat_liq.rho*(1-alpha(m_dot))+sat_vap.rho*alpha(m_dot))
    
    left = lambda m_dot: (sc_liq.rho-rhoH(m_dot))*g*L
    
    right = lambda m_dot: K*2*(0.079/(Re(m_dot)**0.25))*(L/D)*((m_dot**2)/(A**2))*(1/rhoH(m_dot))
    
    mass_fun = lambda m_dot: left(m_dot)-right(m_dot)
    mass_flow = []
    heat = []
    heats = np.linspace(q_start,q_end,100)
    for i in range(0,100):
        q_dot = heats[i]
        heat.append(q_dot)
        
        mass_flow.append(float(optimize.fsolve(mass_fun,0.01)))
#        print(left(mass_flow[i]),right(mass_flow[i]))
#        print(x(mass_flow[i]))
        
    plt.plot(heat,mass_flow)
#    print(q_dot,mass_flow[-1],h_in,h_sat,Hvap)
#    print(x(0.005))
#    print(alpha(x(0.005)))
#    print(left(mass_flow[-1]),right(mass_flow[-1]))
#    print(Re(mass_flow[-1]),rhoH(mass_flow[-1]))
    return heat,mass_flow


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
    h_in = sc_liq.h
    h_sat = sat_liq.h
#    h_gain = ((q_dot)/m_dot)
    
    x_i = (h_in-h_sat)/Hvap
    
    #x = lambda m_dot: (q_dot)/(Hvap*m_dot)
#    x = lambda m_dot: (q_dot-m_dot*cp*(Tsat-Tin))/(Hvap*m_dot)
    x = lambda m_dot: (((q_dot)/m_dot)+h_in-h_sat)/Hvap - x_i*np.exp((((((q_dot)/m_dot)+h_in-h_sat)/Hvap)/x_i)-1)  
    
    left = lambda m_dot: ((sc_liq.rho-sat_liq.rho)+(sat_liq.rho-sat_vap.rho)/(1+(sat_vap.rho/sat_liq.rho)*((1-x(m_dot))/x(m_dot))))*g*L
    
    right = lambda m_dot: (1+x(m_dot)*((sat_liq.rho/sat_vap.rho)-1))*K*(m_dot**2)/(2*sat_liq.rho*A**2)
    
    mass_fun = lambda m_dot: left(m_dot)-right(m_dot)
    mass_flow = []
    heat = []
    heats = np.linspace(q_start,q_end,100)
    for i in range(0,100):
        q_dot = heats[i]
        heat.append(q_dot)
        mass_flow.append(float(optimize.fsolve(mass_fun,0.001)))
        
    plt.plot(heat,mass_flow)
    return heat,mass_flow

def walltemp(z=0.7,K=5,A=0.0005,q_in=2,Tin=90.7,d=0.0254,m_dot = 0.01,Pin=0.101325):
    z = np.linspace(0,z,num=100)
    T_wall = []
    T_zplot = []
    x_thermo = []
    x_levy = []
    alpha_levy = []

    for i in range(100):
        sc_liq_in = IAPWS97(P=0.101325, T=Tin+273.15)
        qflux = (q_in)/(np.pi*d*z[-1])
        
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
        
        T_z = ((qflux)*np.pi*d*z[i])/(m_dot*cp)+Tin
#        print(T_z)
        if T_z > Tsat:
            T_z = Tsat
            
        T_zplot.append(T_z)
        h_in = sc_liq_in.h
        h_sat = sat_liq.h
        h_gain = ((qflux)*np.pi*d*z[i])/m_dot;
#        print(h_gain)

        x_i = (h_in-h_sat)/Hvap

        x_z = (h_gain+h_in-h_sat)/Hvap
        
        x_thermo.append(x_z)
        x = x_z - x_i*np.exp((x_z/x_i)-1) 
#        if T_z == Tsat:
#            x = x_z - x_i*np.exp((x_z/x_i)-1) 

#        else:
#            x = 0

        x_levy.append(x)
        
        alpha = 1/(1+(((1-x)/x)*(sat_vap.rho/sat_liq.rho)))
        alpha_levy.append(alpha)
#            x = 0

        Re_l = (m_dot*d)/(sat_liq.mu*A)

        Pr_l = sat_liq.Prandt
        
        h_l = 0.023*(sat_liq.k/d)*(Re_l**(0.8))*(Pr_l**0.4)
        
        
        h_pool = 55 * (Pstar**(0.12))*((qflux)**(0.6667))*((-np.log10(Pstar))**(-0.55))*M**(-0.5)
#        print(h_l,h_pool)
        F = (1+x*Pr_l*((sat_liq.rho/sat_vap.rho)-1))**(0.35)
        S = (1+0.055*F**(0.1)*Re_l**(0.16))**(-1)
        
        T = lambda T_w: (qflux-(((F*h_l*(T_w-T_z))**2)+((S*h_pool*(T_w-Tsat))**2))**(0.5))
#        Abp = (F*h_l)/(S*h_pool)
#        Aqp = qflux/(S*h_pool*(Tsat-T_z))
        
#        T = lambda T_w: ((Tsat-T_z)/(1+Abp**2))*(1+(1+(1+Abp**2)*(Aqp**2-1))**0.5)+T_z-T_w
        
        if T_z != Tsat:
            T_wall.append(float(optimize.fsolve(T,70)))
        else:
            htp = (((F*h_l)**2)+((S*h_pool)**2))**(0.5)
            T_wall.append((qflux/htp)+Tsat)
#        print(T_z,T_wall[-1])
        
    
#    T_wall = np.asarray(T_wall)
#    T_zplot = np.asarray(T_zplot)
#    plt.plot(x_thermo)
#    plt.plot(x_levy)
#    plt.plot(T_wall)
#    plt.plot(T_zplot)
#    plt.plot(x_plot)
#    plt.show()
    return T_wall, T_zplot,x_thermo,x_levy,alpha_levy,z
#    print(T_wall)

#    plt.plot(T_wall-T_zplot)
    
    
if __name__ == "__main__":

    
#    wut1 = massflow()
#    wut2 = massflow2()
#    massflow3()
    walltemp()
    
    
    
    
    
    
    
    