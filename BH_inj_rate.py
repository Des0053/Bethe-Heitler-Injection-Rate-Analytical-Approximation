import numpy as np
import math 

#-------Essential Constants-------------
c_light=3*10**(10)
m_p=1.67*10**(-24)
m_e=9.1*10**(-28)
sigma_T=6.65*10**(-25)

class BH_inj_rate:
    
    def __init__(self, g, x, g_p, R, ind):
        gpx=g_p+x
        
        self.g=g
        self.x=x
        self.g_p=g_p
        self.R=R
        self.ind=ind
        self.gpx=gpx
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

    def slope(self, gpx):
        params=[0.6586, 0.65, -0.06073489219556636, 0.9893670856189293, 0.94]
        x0, ss, ss2, a_norm, p_s= params
        const=2-ss2*x0**p_s-a_norm/np.exp(1)*x0**(-ss)

        if self.gpx<x0:
            self.slope=2.
        else:
            self.slope=a_norm*self.gpx**(-ss)*np.exp(-self.gpx/x0)+ss2*self.gpx**p_s+const
        
        return self.slope
    
    
    def A_norm(self, x, g_p, R): #Injection Rate Normalization        
        A, B, C, D, E, F, G, p_e, A_1, p_1, c_1=[-2.12, 0.8975274693141311, -0.051033221749056536, 0.999057, -111.9, 1.22, 0.1288184966128324, 1.23, 17.5, 1.8, -3.68493]

        self.A_norm=(A*np.exp(-B*(self.gpx)**2.)+C*(self.gpx)**D+E+p_e*self.x+G*(self.g_p)**(F)+np.log10((10**15/self.R)**(4)*self.tr))+(A_1*(self.gpx)**p_1)*np.exp(c_1*self.gpx)

        return self.A_norm
    
    
    def Q_inj(self, g, x, g_p, R, ind): #Injection Rate For Single proton- Single Photon Interaction
        if type(self.g)==np.float64 or type(self.g)==float:
            self.g=np.array([self.g])
    
        E, gamma, gamma_p=[10**self.x, 10**self.g,  10**self.g_p]
        
        params=[2., (1.23*E)**(-1), 0.47, 0.468, 0.43, 0.465, 0.95, 1., 0.1, 0*10**(-5), 0.14, 0.25, gamma_p*15., 1.75, 1]    
        int_thres, x0, a1, a2, a2_b, a3, cor1a, cor1b, cor2a, cor2b, cor2b_thres, cor2c, x_c, p_lim, AM_flag= params
                    
        A=10**(BH_inj_fun.A_norm(self, self.x, self.g_p, self.R))
        res=np.zeros(len(gamma))
        if gamma_p*E<=2. or gamma_p*E*AM_flag>10**4.:
            pass
        else:
            cont_const=1/np.exp(np.log10(x_c/x0)**self.slope(np.log10(gamma_p*E))*0.5*(1/a2**2-1/a3**2))/np.exp(x_c*cor2b/x0)
            for j in range(0, len(gamma)):
                if gamma[j]>gamma_p*m_p/m_e:
                    pass
                else:
                    power=self.slope
                    if gamma[j]/x0 <= 1.:
                        power=2.
                        cor = cor1a
                        a_slope = a1
                        cor2 = cor2a
                        res[j]=A*np.exp(-(np.log10(gamma[j]/x0))**(power)/( 2.*a_slope**2))*np.exp(-((x0-gamma[j])*cor/gamma[j])**2.)*np.exp(-gamma[j]*cor2/x0)
                    elif gamma[j]<x_c or power>p_lim:
                        if power>p_lim:
                            cor = cor1b
                            cor2 = max(cor2b_thres-0.0665*(gamma_p*E-2.2), 0.007)
                            a_slope = a2
                        else:
                            cor = cor1b
                            cor2 = cor2b
                            a_slope = a2
                        res[j]=A*np.exp(-(np.log10(gamma[j]/x0))**(power)/(2.*a_slope**2))*np.exp(-((x0-gamma[j])*cor/gamma[j])**2.)*np.exp(-cor2*(gamma[j]/x0-1.))
                    else:
                        cor = cor1b
                        a_slope = a3
                        cor2 = cor2c
                        res[j]=A*np.exp(-(np.log10(gamma[j]/x0))**(power)/(2.*a_slope**2))*np.exp(-((x0-gamma[j])*cor/gamma[j])**2.)*np.exp(-cor2*(gamma[j]/x_c-1.))*cont_const
        self.q_BH=res