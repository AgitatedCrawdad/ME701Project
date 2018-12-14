from iapws import IAPWS97
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt

def massflow4(height=1.4,L=0.7,K=5,d=0.0254,q_start=0,q_end=5,Tin=99.7,Pin=0.101325):
    A = np.pi*(d/2)**2
    g = 9.81 # m/s^2
    #Subcooled inlet properties of water   
    sc_liq = IAPWS97(P=Pin, T=Tin+273.15)
    
    #Saturated properties of water
    sat_liq = IAPWS97(P=Pin, x=0)
    
    #Saturated properties of vapor
    sat_vap = IAPWS97(P=Pin, x=1)
    
    Hvap = sat_vap.h-sat_liq.h #kJ/kg
    h_in = sc_liq.h
    h_sat = sat_liq.h
    

    
    x_i = (h_in-h_sat)/Hvap
    
    x = lambda m_dot: (((q_dot)/m_dot)+h_in-h_sat)/Hvap - x_i*np.exp((((((q_dot)/m_dot)+h_in-h_sat)/Hvap)/x_i)-1)  
    
    Re_l = lambda m_dot: (((m_dot*d)/(A*sat_liq.mu)))
    Re_v = lambda m_dot: (((m_dot*d)/(A*sat_vap.mu)))
    
    f_l = lambda m_dot: 0.079/(Re_l(m_dot)**0.25)
    f_v = lambda m_dot: 0.079/(Re_v(m_dot)**0.25)

    xtt = lambda m_dot: (((1-x(m_dot))/(x(m_dot)))**0.9)*((sat_vap.rho/sat_liq.rho)**0.5)*((sat_liq.mu/sat_vap.mu)**0.1)
    
    def tpm(m_dot):
        if Re_l(m_dot)>=2300 and Re_v(m_dot)>=2300:
            C = 20
        elif Re_l(m_dot)<2300 and Re_v(m_dot)>=2300:
            C = 12
        elif Re_l(m_dot)>=2300 and Re_v(m_dot)<2300:
            C = 10
        else:
            C = 5
        if Re_l(m_dot) > 4000:
            tpm_val = 1+(C/xtt(m_dot))+(1/(xtt(m_dot)**2))
        else:
            tpm_val = 1+(C*xtt(m_dot))+(xtt(m_dot)**2)

        
        return tpm_val
        
    alpha = lambda m_dot: 1/(1+(((1-x(m_dot))/x(m_dot))*(sat_vap.rho/sat_liq.rho)))
    
    dpdz_l = lambda m_dot: K*(2*f_l(m_dot)*(m_dot**2))/(d*(A**2)*sat_liq.rho)
    dpdz_v = lambda m_dot: K*(2*f_v(m_dot)*(m_dot**2))/(d*(A**2)*sat_vap.rho)
#    left = lambda m_dot: tpm(m_dot)*K*2*(0.079/(Re_l(m_dot)**0.25))*(L/D)*((m_dot**2)/(A**2))*(1/rhoH(m_dot))
#    left = lambda m_dot: K*dpdz_l(m_dot)*(height)*tpm(m_dot)
    
    def left(m_dot):
        if Re_l(m_dot) > 4000:
            left_val = dpdz_l(m_dot)*(2*L*tpm(m_dot)+(height-L)*tpm(m_dot)+(height-L))
        else:
            left_val = dpdz_v(m_dot)*(2*L*tpm(m_dot)+(height-L)*tpm(m_dot)+(height-L))
            
        return left_val
    
    
    right = lambda m_dot: g*(sat_liq.rho-sat_vap.rho)*((height)*alpha(m_dot))
    
    
    mass_fun = lambda m_dot: left(m_dot)-right(m_dot)
    mass_flow = []
    heat = []
    heats = np.linspace(q_start,q_end,100)
    for i in range(0,100):
        q_dot = heats[i]
        heat.append(q_dot)
        
        if i <5:
            mass_flow.append(float(optimize.fsolve(mass_fun,0.1)))
        else:
            mass_flow.append(float(optimize.fsolve(mass_fun,mass_flow[-1])))      
        
        
    plt.plot(heat,mass_flow)
    return heat,mass_flow


def massflow3(height=1.4,L=0.7,K=5,d=0.0254,q_start=0,q_end=5,Tin=99.7,Pin=0.101325):
    A = np.pi*(d/2)**2
    q_end = 5
    g = 9.81 # m/s^2
    D = 0.0254
    #Subcooled inlet properties of water
    
    sc_liq = IAPWS97(P=Pin, T=Tin+273.15)
    
    #Saturated properties of water
    sat_liq = IAPWS97(P=Pin, x=0)
    
    #Saturated properties of vapor
    sat_vap = IAPWS97(P=Pin, x=1)
    
    Hvap = sat_vap.h-sat_liq.h #kJ/kg


    h_in = sc_liq.h
    h_sat = sat_liq.h
      
    x_i = (h_in-h_sat)/Hvap
    
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
    
    dpdz_l = lambda m_dot: K*(2*f_l(m_dot)*(m_dot**2))/(D*(A**2)*sat_liq.rho)
    
    left = lambda m_dot: dpdz_l(m_dot)*(2*L*tpm(m_dot)+(height-L)*tpm(m_dot)+(height-L))
    
    right = lambda m_dot: g*(sat_liq.rho-sat_vap.rho)*((height)*alpha(m_dot))
    
    
    mass_fun = lambda m_dot: left(m_dot)-right(m_dot)
    mass_flow = []
    heat = []
    heats = np.linspace(q_start,q_end,100)
    for i in range(0,100):
        q_dot = heats[i]
        heat.append(q_dot)
        if i <5:
            mass_flow.append(float(optimize.fsolve(mass_fun,0.01)))
        else:
            mass_flow.append(float(optimize.fsolve(mass_fun,mass_flow[-1])))
        
        
        
        
    plt.plot(heat,mass_flow)
    return heat,mass_flow


def massflow2(L=1.4,K=5,d=0.0254,q_start=0,q_end=5,Tin=99.7,Pin=0.101325):
    A = np.pi*(d/2)**2
    g = 9.81 # m/s^2
    #Subcooled inlet properties of water
#    Tin = 99.4 # deg C
#    m_dot = 0.0010
    
    sc_liq = IAPWS97(P=Pin, T=Tin+273.15)
    
    #Saturated properties of water
    sat_liq = IAPWS97(P=Pin, x=0)
    
    #Saturated properties of vapor
    sat_vap = IAPWS97(P=Pin, x=1)
    
    Hvap = sat_vap.h-sat_liq.h #kJ/kg
    h_in = sc_liq.h
    h_sat = sat_liq.h
    
    x_i = (h_in-h_sat)/Hvap
    
    
    x = lambda m_dot: (((q_dot)/m_dot)+h_in-h_sat)/Hvap - x_i*np.exp((((((q_dot)/m_dot)+h_in-h_sat)/Hvap)/x_i)-1)  

    
    alpha = lambda m_dot: 1/(1+(((1-x(m_dot))/x(m_dot))*(sat_vap.rho/sat_liq.rho)))
    
    Re = lambda m_dot: (((m_dot*d)/(A*(x(m_dot)*sat_vap.mu+((1-x(m_dot))*sat_liq.mu)))))
    
    rhoH = lambda m_dot: (sat_liq.rho*(1-alpha(m_dot))+sat_vap.rho*alpha(m_dot))
    
    left = lambda m_dot: (sc_liq.rho-rhoH(m_dot))*g*L
    
    right = lambda m_dot: K*2*(0.079/(Re(m_dot)**0.25))*(L/d)*((m_dot**2)/(A**2))*(1/rhoH(m_dot))
    
    mass_fun = lambda m_dot: left(m_dot)-right(m_dot)
    mass_flow = []
    heat = []
    heats = np.linspace(q_start,q_end,100)
    for i in range(0,100):
        q_dot = heats[i]
        heat.append(q_dot)

        if i <5:
            mass_flow.append(float(optimize.fsolve(mass_fun,0.01)))
        else:
            mass_flow.append(float(optimize.fsolve(mass_fun,mass_flow[-1])))

    plt.plot(heat,mass_flow)
    return heat,mass_flow


def massflow(L=1.4,K=5,d=0.0005,q_start=0,q_end=5,Tin=99.7,Pin=0.101325):
    g = 9.81 # m/s^2
    A = np.pi*(d/2)**2
    sc_liq = IAPWS97(P=Pin, T=Tin+273.15)
    
    #Saturated properties of water
    sat_liq = IAPWS97(P=Pin, x=0)
    
    #Saturated properties of vapor
    sat_vap = IAPWS97(P=Pin, x=1)
    
    Hvap = sat_vap.h-sat_liq.h #kJ/kg
    
    h_in = sc_liq.h
    h_sat = sat_liq.h
    
    x_i = (h_in-h_sat)/Hvap
    
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

        if i <5:
            mass_flow.append(float(optimize.fsolve(mass_fun,0.01)))
        else:
            mass_flow.append(float(optimize.fsolve(mass_fun,mass_flow[-1])))
            
    plt.plot(heat,mass_flow)
    return heat,mass_flow

def walltemp(z=0.7,K=5,q_in=2,Tin=90.7,d=0.0254,m_dot = 0.01,Pin=0.101325):
    z = np.linspace(0,z,num=100)
    T_wall = []
    T_zplot = []
    x_thermo = []
    x_levy = []
    alpha_levy = []
    A = np.pi*(d/2)**2

    for i in range(100):
        sc_liq_in = IAPWS97(P=Pin, T=Tin+273.15)
        qflux = (q_in)/(np.pi*d*z[-1])
        
        #Saturated properties of water
        sat_liq = IAPWS97(P=Pin, x=0)
        
        #Saturated properties of vapor
        sat_vap = IAPWS97(P=Pin, x=1)
        
        Hvap = sat_vap.h-sat_liq.h #kJ/kg
        cp = sat_liq.cp # kJ/kgK
        Tsat = sat_liq.T-273.15    
        
        M = 18 #Molecular weight
        Pcrit = 220.64 #Critical pressure of water in bars
        
        Pstar = (Pin*10)/Pcrit
        
        T_z = ((qflux)*np.pi*d*z[i])/(m_dot*cp)+Tin

        if T_z > Tsat:
            T_z = Tsat
            
        T_zplot.append(T_z)
        h_in = sc_liq_in.h
        h_sat = sat_liq.h
        h_gain = ((qflux)*np.pi*d*z[i])/m_dot;

        x_i = (h_in-h_sat)/Hvap

        x_z = (h_gain+h_in-h_sat)/Hvap
        
        x_thermo.append(x_z)
        x = x_z - x_i*np.exp((x_z/x_i)-1) 

        x_levy.append(x)
        
        alpha = 1/(1+(((1-x)/x)*(sat_vap.rho/sat_liq.rho)))
        alpha_levy.append(alpha)

        Re_l = (m_dot*d)/(sat_liq.mu*A)

        Pr_l = sat_liq.Prandt
        
        h_l = 0.023*(sat_liq.k/d)*(Re_l**(0.8))*(Pr_l**0.4)
        
        
        h_pool = 55 * (Pstar**(0.12))*((qflux)**(0.6667))*((-np.log10(Pstar))**(-0.55))*M**(-0.5)

        F = (1+x*Pr_l*((sat_liq.rho/sat_vap.rho)-1))**(0.35)
        S = (1+0.055*F**(0.1)*Re_l**(0.16))**(-1)
        
        T = lambda T_w: (qflux-(((F*h_l*(T_w-T_z))**2)+((S*h_pool*(T_w-Tsat))**2))**(0.5))
      
        if T_z != Tsat:
            T_wall.append(float(optimize.fsolve(T,70)))
        else:
            htp = (((F*h_l)**2)+((S*h_pool)**2))**(0.5)
            T_wall.append((qflux/htp)+Tsat)

        
    

    return T_wall, T_zplot,x_thermo,x_levy,alpha_levy,z

    
    
if __name__ == "__main__":

    
    massflow()
    massflow2()
    massflow3()
    massflow4()
    walltemp()
    
    

    
    