import numpy as np
import math 

#-------Essential Constants-------------
c_light=3*10**(10)
m_p=1.67*10**(-24)
m_e=9.1*10**(-28)
sigma_T=6.65*10**(-25)

class BH_inj_rate_arr:
    
    def __init__(self, g, x, g_p, R, ind):
        log_gpE=g_p+x
        
        self.g=g  #pair Lorentz factor in logarithm
        self.x=x  #photon energy in m_ec^2 units, in logarithm
        self.g_p=g_p #proton Lorentz factor in logarithm
        self.R=R #source radius in linera scale
        self.ind=ind #index that defines the units of the injection rate
        self.log_gpE=log_gpE #interaction energy γ_p*E in logarithm
        self.q_BH=[]
  
        if self.ind==0:# options: "0" for 1/s units , "1" for 1/(cm^3 s) units, "2" for ATHEvA 1/(V t) units  
            self.tr=4.*math.pi*c_light*R/(3.*sigma_T)
        elif self.ind==1:
            self.tr=c_light/(sigma_T*R**2)
        elif self.ind==2:
            self.tr=1.
            
        self.slope
        self.A_norm
        self.Q_inj(g, x, g_p, R, ind)

    def slope(self, log_gpE): #Injection Rate slope p(γ_p*E)
        params=[0.6586, 0.65, -0.06073489219556636, 0.9893670856189293, 0.94]
        x0, ss, b, a, ss2= params
        const=2-b*x0**ss2-a/np.exp(1)*x0**(-ss)

        if log_gpE<x0:
            s=2.
        else:
            s=a*log_gpE**(-ss)*np.exp(-log_gpE/x0)+b*log_gpE**ss2+const
        
        return s
    
    
    def A_norm(self, x, g_p, R): #Injection Rate Normalization        
        log_gpE=g_p+x
        D, d, F, f, E, g, G, J, H, h_1, h_2=[-2.12, 0.8975274693141311, -0.051033221749056536, 0.999057, -111.9, 1.22, 0.1288184966128324, 1.23, 17.5, 1.8, -3.68493]

        A=(D*np.exp(-d*(log_gpE)**2.)+F*(log_gpE)**f+E+J*x+G*(g_p)**(g)+np.log10((10**15/R)**(4)*self.tr))+(H*(log_gpE)**h_1)*np.exp(h_2*log_gpE)

        return A
    
    
    def Q_inj(self, g, x, g_p, R, ind): #Injection Rate For Single proton- Single Photon Interaction
        if type(self.g)==np.float64 or type(self.g)==float:
            self.g=np.array([self.g])
    
        E, gamma, gamma_p=[10**self.x, 10**self.g,  10**self.g_p]
        
        params=[2., (1.23*E)**(-1), 0.47, 0.468, 0.465, 0.95, 1., 0.1, 0*10**(-5), 0.14, 0.25, gamma_p*15., 1.75, 1]    
        int_thres, ge_pk, a1_1, a1_2, a1_3, a2_1, a2_2, a3_1, a3_2, a3_3, a3_4, gp_c, p_lim, AM_flag= params
        #a1_i are all the branches of the a1 parameter
        #a2_i are all the branches of the a2 parameter
        #a3_i are all the branches of the a3 parameter
        
        res=np.zeros([len(E), len(gamma)])
        for i in range(0,len(E)):
            A=10**(self.A_norm(self.x[i], self.g_p, self.R))
            if gamma_p*E[i]<=2. or gamma_p*E[i]*AM_flag>10**4.:
                pass
            else:
                cont_const=1/np.exp(np.log10(gp_c/ge_pk[i])**self.slope(np.log10(gamma_p*E[i]))*0.5*(1/a1_2**2-1/a1_3**2))/np.exp(gp_c*a3_2/ge_pk[i])
                for j in range(0, len(gamma)):
                    if gamma[j]>gamma_p*m_p/m_e:
                        pass
                    else:
                        power=self.slope(np.log10(gamma_p*E[i]))
                        if gamma[j]/ge_pk[i] <= 1.:
                            power=2.
                            a1 = a1_1
                            a2 = a2_1
                            a3 = a3_1
                            res[i][j]=A*np.exp(-(np.log10(gamma[j]/ge_pk[i]))**(power)/( 2.*a1**2))*np.exp(-((ge_pk[i]-gamma[j])*a2/gamma[j])**2.)*np.exp(-gamma[j]*a3/ge_pk[i])
                        elif gamma[j]<gp_c or power>p_lim:
                            if power>p_lim:
                                a1 = a1_2
                                a2 = a2_2
                                a3 = max(a3_3-0.0665*(gamma_p*E[i]-2.2), 0.007)
                            else:
                                a1 = a1_2
                                a2 = a2_2
                                a3 = a3_2
                            res[i][j]=A*np.exp(-(np.log10(gamma[j]/ge_pk[i]))**(power)/(2.*a1**2))*np.exp(-((ge_pk[i]-gamma[j])*a2/gamma[j])**2.)*np.exp(-a3*(gamma[j]/ge_pk[i]-1.))
                        else:
                            a1 = a1_3
                            a2 = a2_2
                            a3 = a3_4
                            res[i][j]=A*np.exp(-(np.log10(gamma[j]/ge_pk[i]))**(power)/(2.*a1**2))*np.exp(-((ge_pk[i]-gamma[j])*a2/gamma[j])**2.)*np.exp(-a3*(gamma[j]/gp_c-1.))*cont_const
        if len(E)==1:
            self.q_BH=res[0]
        else:
            self.q_BH=res
