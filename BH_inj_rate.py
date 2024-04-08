import numpy as np
import math 

#-------Essential Constants-------------
c_light=3*10**(10)
m_p=1.67*10**(-24)
m_e=9.1*10**(-28)
sigma_T=6.65*10**(-25)

class BH_inj_rate:
    
    def __init__(self, g, x, g_p, R, ind):
        log_gpE=g_p+x
        
        self.g=g
        self.x=x
        self.g_p=g_p
        self.R=R
        self.ind=ind
        self.log_gpE=log_gpE
        self.q_BH=[]
  
        if self.ind==0:# options: "cgs_rate" for 1/s units, "cgs_dens_rate" for 1/(cm^3 s) units, "AM_dens_rate" for ATHEvA 1/(V t) units    
            self.tr=4.*math.pi*c_light*R/(3.*sigma_T)
        elif self.ind==1:
            self.tr=c_light/(sigma_T*R**2)
        elif self.ind==2:
            self.tr=1.
            
        self.slope
        self.A_norm
        self.Q_inj(g, x, g_p, R, ind)

    def slope(self, log_gpE):
        params=[0.6586, 0.65, -0.06073489219556636, 0.9893670856189293, 0.94]
        x0, ss, b, a, ss2= params
        const=2-b*x0**ss2-a/np.exp(1)*x0**(-ss)

        if self.log_gpE<x0:
            self.slope=2.
        else:
            self.slope=a*self.log_gpE**(-ss)*np.exp(-self.log_gpE/x0)+b*self.log_gpE**ss2+const
        
        return self.slope
    
    
    def A_norm(self, x, g_p, R): #Injection Rate Normalization        
        D, d, F, f, E, g, G, J, H, h_1, h_2=[-2.12, 0.8975274693141311, -0.051033221749056536, 0.999057, -111.9, 1.22, 0.1288184966128324, 1.23, 17.5, 1.8, -3.68493]

        self.A_norm=(D*np.exp(-d*(self.log_gpE)**2.)+F*(self.log_gpE)**f+E+J*self.x+G*(self.g_p)**(g)+np.log10((10**15/self.R)**(4)*self.tr))+(H*(self.log_gpE)**h_1)*np.exp(h_2*self.log_gpE)

        return self.A_norm
    
    
    def Q_inj(self, g, x, g_p, R, ind): #Injection Rate For Single proton- Single Photon Interaction
        if type(self.g)==np.float64 or type(self.g)==float:
            self.g=np.array([self.g])
    
        E, gamma, gamma_p=[10**self.x, 10**self.g,  10**self.g_p]
        
        params=[2., (1.23*E)**(-1), 0.47, 0.468, 0.465, 0.95, 1., 0.1, 0*10**(-5), 0.14, 0.25, gamma_p*15., 1.75, 1]    
        int_thres, ge_pk, a1_1, a1_2, a1_3, a2_1, a2_2, a3_1, a3_2, a3_3, a3_4, ge_c, p_lim, AM_flag= params
                    
        A=10**(BH_inj_rate.A_norm(self, self.x, self.g_p, self.R))
        res=np.zeros(len(gamma))
        if gamma_p*E<=2. or gamma_p*E*AM_flag>10**4.:
            pass
        else:
            cont_const=1/np.exp(np.log10(ge_c/ge_pk)**self.slope(np.log10(gamma_p*E))*0.5*(1/a1_2**2-1/a1_3**2))/np.exp(ge_c*a3_2/ge_pk)
            for j in range(0, len(gamma)):
                if gamma[j]>gamma_p*m_p/m_e:
                    pass
                else:
                    power=self.slope
                    if gamma[j]/ge_pk <= 1.:
                        power=2.
                        a1 = a1_1
                        a2 = a2_1
                        a3 = a3_1
                        res[j]=A*np.exp(-(np.log10(gamma[j]/ge_pk))**(power)/( 2.*a1**2))*np.exp(-((ge_pk-gamma[j])*a2/gamma[j])**2.)*np.exp(-gamma[j]*a3/ge_pk)
                    elif gamma[j]<ge_c or power>p_lim:
                        if power>p_lim:
                            a1 = a1_2
                            a2 = a2_2
                            a3 = max(a3_3-0.0665*(gamma_p*E-2.2), 0.007)
                        else:
                            a1 = a1_2
                            a2 = a2_2
                            a3 = a3_2
                        res[j]=A*np.exp(-(np.log10(gamma[j]/ge_pk))**(power)/(2.*a1**2))*np.exp(-((ge_pk-gamma[j])*a2/gamma[j])**2.)*np.exp(-a3*(gamma[j]/ge_pk-1.))
                    else:
                        a1 = a1_3
                        a2 = a2_2
                        a3 = a3_4
                        res[j]=A*np.exp(-(np.log10(gamma[j]/ge_pk))**(power)/(2.*a1**2))*np.exp(-((ge_pk-gamma[j])*a2/gamma[j])**2.)*np.exp(-a3*(gamma[j]/ge_c-1.))*cont_const
        self.q_BH=res
