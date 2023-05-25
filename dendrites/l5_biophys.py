# -*- coding: utf-8 -*-
"""
Created on Fri May 12 12:59:10 2023

@author: xiao208
"""


import numpy as np



class NaTa_t:
    def __init__(self, v, E = 50.0):
        self.m = self.minf(v)
        self.h = self.hinf(v)
        self.qt = 2.3**((34-21)/10)
        self.E = E

    def update(self, v, dt):
        # self.m = self.m + (self.minf(v)-self.m)/self.mtau(v)*dt
        # self.h = self.h + (self.hinf(v)-self.h)/self.htau(v)*dt
        self.m = self.m + (1-np.exp(-dt/self.mtau(v)))*(self.minf(v)-self.m)
        self.h = self.h + (1-np.exp(-dt/self.htau(v)))*(self.hinf(v)-self.h)


    def malpha(self, v):
        return (0.182 * (v- (-38.)))/(1-(np.exp(-(v + 38.)/6.)))

    def d_malpha(self, v):
        return 0.182*((1. - np.exp(-(v - (-38.))/6.)) - 1/6.*(v - (-38.))*np.exp(-(v - (-38.))/6.))/(1. - np.exp(-(v - (-38.))/6.))**2

    def mbeta(self, v):
        return (-0.124 * (v + 38.))/(1-(np.exp((v + 38.)/6.)))

    def d_mbeta(self, v):
        return -0.124*((1. - np.exp((v + 38.)/6.)) + 1/6.*(v + 38.)*np.exp((v + 38.)/6.))/(1. - np.exp((v + 38.)/6.))**2

    def mtau(self,v):
        return (1/(self.malpha(v) + self.mbeta(v)))/self.qt
    
    def minf(self, v):
        return self.malpha(v)/(self.malpha(v) + self.mbeta(v))
    
    def halpha(self, v):
        return (-0.015 * (v + 66.))/(1-(np.exp((v + 66.)/6.)))

    def d_halpha(self, v):
        return -0.015*((1. - np.exp((v + 66.)/6.)) + 1/6.*(v + 66.)*np.exp((v + 66.)/6.))/(1. - np.exp((v + 66.)/6.))**2

    def hbeta(self, v):
        return (0.015 * (v + 66.))/(1-(np.exp(-(v + 66.)/6.)))

    def d_hbeta(self, v):
        return 0.015*((1. - np.exp(-(v - (-66.))/6.)) - 1/6.*(v - (-66.))*np.exp(-(v - (-66.))/6.))/(1. - np.exp(-(v - (-66.))/6.))**2

    def htau(self, v):
        return (1/(self.halpha(v) + self.hbeta(v)))/self.qt
    
    def hinf(self, v):
        return self.halpha(v)/(self.halpha(v) + self.hbeta(v))

    def g_s(self, m, h):
        """
        scaling term for conductance i = gbar*g_s(t, v)

        """
        return m**3*h

    def d_minf(self, v):
        ma = self.malpha(v)
        mb = self.mbeta(v)
        d_ma = self.d_malpha(v)
        d_mb = self.d_mbeta(v)
        return (mb*d_ma-ma*d_mb)/((ma+mb)**2)

    def d_mtau(self, v):
        ma = self.malpha(v)
        mb = self.mbeta(v)
        d_ma = self.d_malpha(v)
        d_mb = self.d_mbeta(v)
        return -1/self.qt*(d_ma + d_mb)/((ma + mb)**2)

    def d_hinf(self, v):
        ha = self.halpha(v)
        hb = self.hbeta(v)
        d_ha = self.d_halpha(v)
        d_hb = self.d_hbeta(v)
        return (hb*d_ha-ha*d_hb)/((ha+hb)**2)

    def d_htau(self, v):
        ha = self.halpha(v)
        hb = self.hbeta(v)
        d_ha = self.d_halpha(v)
        d_hb = self.d_hbeta(v)
        return -1/self.qt*(d_ha + d_hb)/((ha + hb)**2)

    def m_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        d_minf = self.d_minf(v)
        mtau = self.mtau(v)
        d_mtau = self.d_mtau(v)
        minf = self.minf(v)
        return (d_minf/mtau)-(minf-self.m)/(mtau**2)*d_mtau

    def m_b(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        return 1/self.mtau(v)

    def h_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        d_hinf = self.d_hinf(v)
        htau = self.htau(v)
        d_htau = self.d_htau(v)
        hinf = self.hinf(v)
        return (d_hinf/htau)-(hinf-self.h)/(htau**2)*d_htau

    def h_b(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        return 1/self.htau(v)

    def m_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_m

        """
        return 3*(m**2*h)*(self.E - v)

    def h_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_h

        """
        return (m**3)*(self.E - v)


class NaTs2_t:
    def __init__(self, v, E = 50.0):
        self.m = self.minf(v)
        self.h = self.hinf(v)
        self.qt = 2.3**((34-21)/10)
        self.E = E

    def update(self, v, dt):
        # self.m = self.m + (self.minf(v)-self.m)/self.mtau(v)*dt
        # self.h = self.h + (self.hinf(v)-self.h)/self.htau(v)*dt
        self.m = self.m + (1-np.exp(-dt/self.mtau(v)))*(self.minf(v)-self.m)
        self.h = self.h + (1-np.exp(-dt/self.htau(v)))*(self.hinf(v)-self.h)


    def malpha(self, v):
        # i = np.where(v==-32)[0]
        # v[i] = v[i] + 0.0001
        return (0.182 * (v- (-32.)))/(1-(np.exp(-(v + 32.)/6.)))

    def d_malpha(self, v):
        # i = np.where(v==-32)[0]
        # v[i] = v[i] + 0.0001
        return 0.182*((1. - np.exp(-(v - (-32.))/6.)) - 1/6.*(v - (-32.))*np.exp(-(v - (-32.))/6.))/(1. - np.exp(-(v - (-32.))/6.))**2

    def mbeta(self, v):
        # i = np.where(v==-32)[0]
        # v[i] = v[i] + 0.0001
        return (-0.124 * (v + 32.))/(1-(np.exp((v + 32.)/6.)))

    def d_mbeta(self, v):
        # i = np.where(v==-32)[0]
        # v[i] = v[i] + 0.0001
        return -0.124*((1. - np.exp((v + 32.)/6.)) + 1/6.*(v + 32.)*np.exp((v + 32.)/6.))/(1. - np.exp((v + 32.)/6.))**2

    def mtau(self,v):
        return (1/(self.malpha(v) + self.mbeta(v)))/self.qt
    
    def minf(self, v):
        return self.malpha(v)/(self.malpha(v) + self.mbeta(v))
    
    def halpha(self, v):
        # i = np.where(v==-60)[0]
        # v[i] = v[i] + 0.0001
        return (-0.015 * (v + 60.))/(1-(np.exp((v + 60.)/6.)))

    def d_halpha(self, v):
        # i = np.where(v==-60)[0]
        # v[i] = v[i] + 0.0001
        return -0.015*((1. - np.exp((v + 60.)/6.)) + 1/6.*(v + 60.)*np.exp((v + 60.)/6.))/(1. - np.exp((v + 60.)/6.))**2

    def hbeta(self, v):
        # i = np.where(v==-60)[0]
        # v[i] = v[i] + 0.0001
        return (0.015 * (v + 60.))/(1-(np.exp(-(v + 60.)/6.)))

    def d_hbeta(self, v):
        # i = np.where(v==-60)[0]
        # v[i] = v[i] + 0.0001
        return 0.015*((1. - np.exp(-(v - (-60.))/6.)) - 1/6.*(v - (-60.))*np.exp(-(v - (-60.))/6.))/(1. - np.exp(-(v - (-60.))/6.))**2

    def htau(self, v):
        return (1/(self.halpha(v) + self.hbeta(v)))/self.qt
    
    def hinf(self, v):
        return self.halpha(v)/(self.halpha(v) + self.hbeta(v))

    def g_s(self, m, h):
        """
        scaling term for conductance i = gbar*g_s(t, v)

        """
        return m**3*h

    def d_minf(self, v):
        ma = self.malpha(v)
        mb = self.mbeta(v)
        d_ma = self.d_malpha(v)
        d_mb = self.d_mbeta(v)
        return (mb*d_ma-ma*d_mb)/((ma+mb)**2)

    def d_mtau(self, v):
        ma = self.malpha(v)
        mb = self.mbeta(v)
        d_ma = self.d_malpha(v)
        d_mb = self.d_mbeta(v)
        return -1/self.qt*(d_ma + d_mb)/((ma + mb)**2)

    def d_hinf(self, v):
        ha = self.halpha(v)
        hb = self.hbeta(v)
        d_ha = self.d_halpha(v)
        d_hb = self.d_hbeta(v)
        return (hb*d_ha-ha*d_hb)/((ha+hb)**2)

    def d_htau(self, v):
        ha = self.halpha(v)
        hb = self.hbeta(v)
        d_ha = self.d_halpha(v)
        d_hb = self.d_hbeta(v)
        return -1/self.qt*(d_ha + d_hb)/((ha + hb)**2)

    def m_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        d_minf = self.d_minf(v)
        mtau = self.mtau(v)
        d_mtau = self.d_mtau(v)
        minf = self.minf(v)
        return (d_minf/mtau)-(minf-m)/(mtau**2)*d_mtau

    def m_b(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        return 1/self.mtau(v)

    def h_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        d_hinf = self.d_hinf(v)
        htau = self.htau(v)
        d_htau = self.d_htau(v)
        hinf = self.hinf(v)
        return (d_hinf/htau)-(hinf-h)/(htau**2)*d_htau

    def h_b(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        return 1/self.htau(v)

    def m_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_m

        """
        return 3*(m**2*h)*(self.E - v)

    def h_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_h

        """
        return (m**3)*(self.E - v)


class Nap_Et2:
    def __init__(self, v, E = 50.0):
        self.m = self.minf(v)
        self.h = self.hinf(v)
        self.qt = 2.3**((34-21)/10)
        self.E = E

    def update(self, v, dt):
        # self.m = self.m + (self.minf(v)-self.m)/self.mtau(v)*dt
        # self.h = self.h + (self.hinf(v)-self.h)/self.htau(v)*dt
        self.m = self.m + (1-np.exp(-dt/self.mtau(v)))*(self.minf(v)-self.m)
        self.h = self.h + (1-np.exp(-dt/self.htau(v)))*(self.hinf(v)-self.h)

    def malpha(self, v):
        # i = np.where(v==-38)[0]
        # v[i] = v[i] + 0.0001
        return (0.182 * (v- (-38.)))/(1-(np.exp(-(v + 38.)/6.)))

    def d_malpha(self, v):
        # i = np.where(v==-38)[0]
        # v[i] = v[i] + 0.0001
        return 0.182*((1. - np.exp(-(v - (-38.))/6.)) - 1/6.*(v - (-38.))*np.exp(-(v - (-38.))/6.))/(1. - np.exp(-(v - (-38.))/6.))**2

    def mbeta(self, v):
        # i = np.where(v==-38)[0]
        # v[i] = v[i] + 0.0001
        return (-0.124 * (v + 38.))/(1-(np.exp((v + 38.)/6.)))

    def d_mbeta(self, v):
        # i = np.where(v==-38)[0]
        # v[i] = v[i] + 0.0001
        return -0.124*((1. - np.exp((v + 38.)/6.)) + 1/6.*(v + 38.)*np.exp((v + 38.)/6.))/(1. - np.exp((v + 38.)/6.))**2

    def mtau(self,v):
        return 6*(1/(self.malpha(v) + self.mbeta(v)))/self.qt
    
    def minf(self, v):
        return 1/(1+np.exp((v + 52.6)/(-4.6))) # self.m_alpha/(self.m_alpha + self.m_beta)
    
    def halpha(self, v):
        # i = np.where(v==-17)[0]
        # v[i] = v[i] + 0.0001
        return -2.88e-6 * (v + 17)/(1-np.exp((v+17)/4.63))

    def d_halpha(self, v):
        # i = np.where(v==-17)[0]
        # v[i] = v[i] + 0.0001
        return -2.88e-6*((1-np.exp((v+17)/4.63)) + 1/4.63*(v + 17)*np.exp((v + 17)/4.63))/(1 - np.exp((v+17)/4.63))**2

    def hbeta(self, v):
        # i = np.where(v==-17)[0]
        # v[i] = v[i] + 0.0001
        return 6.94e-6 * (v+ 64.4)/(1 - np.exp(-(v + 64.4)/2.63))

    def d_hbeta(self, v):
        # i = np.where(v==-17)[0]
        # v[i] = v[i] + 0.0001
        return 6.94e-6*((1-np.exp(-(v + 64.4)/2.63)) - 1/2.63*(v+ 64.4)*np.exp(-(v + 64.4)/2.63))/(1 - np.exp(-(v + 64.4)/2.63))**2

    def htau(self, v):
        return (1/(self.halpha(v) + self.hbeta(v)))/self.qt
    
    def hinf(self, v):
        return 1/(1+np.exp((v + 48.8)/10)) # self.h_alpha/(self.h_alpha + self.h_beta)

    def g_s(self, m, h):
        """
        scaling term for conductance i = gbar*g_s(t, v)

        """
        return m**3*h

    def d_minf(self, v):
        return 1/4.6*np.exp(-(v + 52.6)/4.6)/(1 + np.exp(-(v + 52.6)/4.6))**2

    def d_mtau(self, v):
        ma = self.malpha(v)
        mb = self.mbeta(v)
        d_ma = self.d_malpha(v)
        d_mb = self.d_mbeta(v)
        return -1/self.qt*(d_ma + d_mb)/((ma + mb)**2)

    def d_hinf(self, v):
        return 1/10*np.exp(-(v + 48.8)/10)/(1 + np.exp(-(v + 48.8)/10))**2

    def d_htau(self, v):
        ha = self.halpha(v)
        hb = self.hbeta(v)
        d_ha = self.d_halpha(v)
        d_hb = self.d_hbeta(v)
        return -1/self.qt*(d_ha + d_hb)/((ha + hb)**2)

    def m_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        d_minf = self.d_minf(v)
        mtau = self.mtau(v)
        d_mtau = self.d_mtau(v)
        minf = self.minf(v)
        return (d_minf/mtau)-(minf-m)/(mtau**2)*d_mtau

    def m_b(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        return 1/self.mtau(v)

    def h_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        d_hinf = self.d_hinf(v)
        htau = self.htau(v)
        d_htau = self.d_htau(v)
        hinf = self.hinf(v)
        return (d_hinf/htau)-(hinf-h)/(htau**2)*d_htau

    def h_b(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        return 1/self.htau(v)

    def m_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_m

        """
        return 3*(m**2*h)*(self.E - v)

    def h_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_h

        """
        return (m**3)*(self.E - v)

class SKv3_1:
    def __init__(self, v, E = -85.0):
        self.m = self.minf(v)
        self.qt = 2.3**((34-21)/10)
        self.E = E

    def update(self, v, dt):
        self.m = self.m + (1-np.exp(-dt/self.mtau(v)))*(self.minf(v)-self.m)

    def mtau(self,v):
        return 4/(1 + np.exp((v + 46.56)/(-44.14)))
    
    def minf(self, v):
        return 1/(1+np.exp((v -18.7)/(-9.7))) # self.m_alpha/(self.m_alpha + self.m_beta)

    def d_minf(self, v):
        return 1/9.7*np.exp(-(v - 18.7)/9.7)/(1 + np.exp(-(v - 18.7)/9.7))**2

    def d_mtau(self, v):
        return 4/44.14*np.exp(-(v + 46.56)/44.14)/(1 + np.exp(-(v + 46.56)/44.14))**2

    def g_s(self, m):
        """
        scaling term for conductance i = gbar*g_s(t, v)

        """
        return m

    def m_a(self, v, m):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        d_minf = self.d_minf(v)
        mtau = self.mtau(v)
        d_mtau = self.d_mtau(v)
        minf = self.minf(v)
        return (d_minf/mtau)-(minf-m)/(mtau**2)*d_mtau

    def m_b(self, v, m):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        return 1/self.mtau(v)

    def m_c(self, v, m):
        """
        scaling term for partial gradiant computation s_dv = c*s_m

        """
        return self.E - v

class Im:
    def __init__(self, v, E = -85.0):
        self.m = self.minf(v)
        self.qt = 2.3**((34-21)/10)
        self.E = E

    def update(self, v, dt):
        self.m = self.m + (1-np.exp(-dt/self.mtau(v)))*(self.minf(v)-self.m)

    def minf(self, v):
        return self.malpha(v)/(self.malpha(v) + self.mbeta(v))


    def mtau(self, v):
    	return 1/(self.malpha(v) + self.mbeta(v))/self.qt
    
    
    def malpha(self, v):
    	return 3.3e-3 *np.exp(0.1*(v + 35))

    
    def mbeta(self, v):
    	return 3.3e-3 *np.exp(-0.1*(v + 35))
    
    
    def d_minf(self, v):
        ma = self.malpha(v)
        mb = self.mbeta(v)
        d_ma = self.d_malpha(v)
        d_mb = self.d_mbeta(v)
        return (mb*d_ma-ma*d_mb)/((ma+mb)**2)
    
    
    def d_mtau(self, v):
        ma = self.malpha(v)
        mb = self.mbeta(v)
        d_ma = self.d_malpha(v)
        d_mb = self.d_mbeta(v)
        return -1/self.qt*(d_ma + d_mb)/((ma + mb)**2)
    
    
    def d_malpha(self, v):
    	return 3.3e-4*np.exp(0.1*(v + 35))
    
    
    def d_mbeta(self, v):
    	return -3.3e-4*np.exp(-0.1*(v + 35))

    def g_s(self, m):
        """
        scaling term for conductance i = gbar*g_s(t, v)

        """
        return m

    def m_a(self, v, m):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        d_minf = self.d_minf(v)
        mtau = self.mtau(v)
        d_mtau = self.d_mtau(v)
        minf = self.minf(v)
        return (d_minf/mtau)-(minf-m)/(mtau**2)*d_mtau

    def m_b(self, v, m):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        return 1/self.mtau(v)

    def m_c(self, v, m):
        """
        scaling term for partial gradiant computation s_dv = c*s_m

        """
        return self.E - v

class Ih:
    def __init__(self, v, E = -45.0):
        self.m = self.minf(v)
        self.qt = 2.3**((34-21)/10)
        self.E = E

    def update(self, v, dt):
        self.m = self.m + (1-np.exp(-dt/self.mtau(v)))*(self.minf(v)-self.m)

    def malpha(self, v):
        return 0.00643*(v + 154.9)/(np.exp((v + 154.9)/11.9) - 1)


    def d_malpha(self, v):
        return 0.00643*(np.exp((v + 154.9)/11.9) - 1 - 1/11.9*(v + 154.9)*
                        np.exp((v + 154.9)/11.9))/(np.exp((v + 154.9)/11.9) - 1)**2
    
    
    def mbeta(self, v):
        return 0.193*np.exp(v/33.1)
    

    def d_mbeta(self, v):
        return 0.193/33.1*np.exp(v/33.1)
    
    
    def minf(self, v):
        return self.malpha(v)/(self.malpha(v) + self.mbeta(v))
    
    
    def mtau(self, v):
        return 1/(self.malpha(v) + self.mbeta(v))

    def g_s(self, m):
        """
        scaling term for conductance i = gbar*g_s(t, v)

        """
        return m

    def m_a(self, v, m):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        d_ma = self.d_malpha(v)
        d_mb = self.d_mbeta(v)
        return d_ma*(1-m)-d_mb*m


    def m_b(self, v, m):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        ma = self.malpha(v)
        mb = self.mbeta(v)
        return ma + mb

    def m_c(self, v, m):
        """
        scaling term for partial gradiant computation s_dv = c*s_m

        """
        return self.E - v

class Ca_HVA:
    def __init__(self, v, E = 127.0):
        self.m = self.minf(v)
        self.h = self.hinf(v)
        self.qt = 2.3**((34-21)/10)
        self.E = E

    def update(self, v, dt):
        # self.m = self.m + (self.minf(v)-self.m)/self.mtau(v)*dt
        # self.h = self.h + (self.hinf(v)-self.h)/self.htau(v)*dt
        self.m = self.m + (1-np.exp(-dt/self.mtau(v)))*(self.minf(v)-self.m)
        self.h = self.h + (1-np.exp(-dt/self.htau(v)))*(self.hinf(v)-self.h)


    def malpha(self, v):
        # i = np.where(v==-27)[0]
        # v[i] = v[i] + 0.0001
        return (0.055 * (v + 27.))/(1-(np.exp(-(v + 27.)/3.8)))

    def d_malpha(self, v):
        # i = np.where(v==-27)[0]
        # v[i] = v[i] + 0.0001
        return 0.055*((1 - np.exp(-(v + 27)/3.8)) - 1/3.8*(v + 27)*np.exp(-(v + 27)/3.8))/(1 - np.exp(-(v + 27)/3.8))**2

    def mbeta(self, v):
        # i = np.where(v==-27)[0]
        # v[i] = v[i] + 0.0001
        return 0.94 * np.exp(-(v + 75)/17)

    def d_mbeta(self, v):
        # i = np.where(v==-27)[0]
        # v[i] = v[i] + 0.0001
        return -0.94/17 * np.exp(-(v + 75)/17)

    def mtau(self,v):
        return 1/(self.malpha(v) + self.mbeta(v))

    def minf(self, v):
        return self.malpha(v)/(self.malpha(v) + self.mbeta(v))
    
    def halpha(self, v):
        return 0.000457*np.exp(-(v + 13)/50)

    def d_halpha(self, v):
        return -0.000457/50*np.exp(-(v + 13)/50)

    def hbeta(self, v):
        return 0.0065/(1 + np.exp(-(v + 15)/28))

    def d_hbeta(self, v):
        return - 0.0065/28*np.exp(-(v + 15)/28)/(1 + np.exp(-(v + 15)/28))**2

    def htau(self, v):
        return 1/(self.halpha(v) + self.hbeta(v))

    def hinf(self, v):
        return self.halpha(v)/(self.halpha(v) + self.hbeta(v))

    def g_s(self, m, h):
        """
        scaling term for conductance i = gbar*g_s(t, v)

        """
        return m**2*h

    def d_minf(self, v):
        ma = self.malpha(v)
        mb = self.mbeta(v)
        d_ma = self.d_malpha(v)
        d_mb = self.d_mbeta(v)
        return (mb*d_ma-ma*d_mb)/((ma+mb)**2)

    def d_mtau(self, v):
        ma = self.malpha(v)
        mb = self.mbeta(v)
        d_ma = self.d_malpha(v)
        d_mb = self.d_mbeta(v)
        return -1*(d_ma + d_mb)/((ma + mb)**2)

    def d_hinf(self, v):
        ha = self.halpha(v)
        hb = self.hbeta(v)
        d_ha = self.d_halpha(v)
        d_hb = self.d_hbeta(v)
        return (hb*d_ha-ha*d_hb)/((ha+hb)**2)

    def d_htau(self, v):
        ha = self.halpha(v)
        hb = self.hbeta(v)
        d_ha = self.d_halpha(v)
        d_hb = self.d_hbeta(v)
        return -1*(d_ha + d_hb)/((ha + hb)**2)

    def m_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        d_minf = self.d_minf(v)
        mtau = self.mtau(v)
        d_mtau = self.d_mtau(v)
        minf = self.minf(v)
        return (d_minf/mtau)-(minf-m)/(mtau**2)*d_mtau

    def m_b(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        return 1/self.mtau(v)

    def h_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        d_hinf = self.d_hinf(v)
        htau = self.htau(v)
        d_htau = self.d_htau(v)
        hinf = self.hinf(v)
        return (d_hinf/htau)-(hinf-h)/(htau**2)*d_htau

    def h_b(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        return 1/self.htau(v)

    def m_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_m

        """
        return 2*(m*h)*(self.E - v)

    def h_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_h

        """
        return (m**2)*(self.E - v)

class Ca_LVAst:
    def __init__(self, v, E = 127.0):
        self.m = self.minf(v)
        self.h = self.hinf(v)
        self.qt = 2.3**((34-21)/10)
        self.E = E

    def update(self, v, dt):
        self.m = self.m + (1-np.exp(-dt/self.mtau(v)))*(self.minf(v)-self.m)
        self.h = self.h + (1-np.exp(-dt/self.htau(v)))*(self.hinf(v)-self.h)

    def minf(self, v):
        return 1/(1 + np.exp((v + 30)/(-6))) # self.m_alpha/(self.m_alpha + self.m_beta)

    def mtau(self,v):
        return (5.0000 + 20.0000/(1+np.exp((v + 25.000)/5)))/self.qt

    def d_minf(self, v):
        return 1/6*np.exp(-(v + 30)/6)/(1 + np.exp(-(v + 30)/6))**2

    def d_mtau(self, v):
        return -4/self.qt*np.exp((v + 25)/5)/(1 + np.exp((v + 25)/5))**2

    def hinf(self, v):
        return 1/(1 + np.exp((v + 80)/6.4)) # self.m_alpha/(self.m_alpha + self.m_beta)

    def htau(self,v):
        return (20.0 + 50.0/(1+np.exp((v + 40.000)/7)))/self.qt

    def d_hinf(self, v):
        return -1/6.4*np.exp((v + 80)/6.4)/(1 + np.exp((v + 80)/6.4))**2

    def d_htau(self, v):
        return -50/7/self.qt*np.exp((v + 40)/7)/(1 + np.exp((v + 40)/7))**2

    def g_s(self, m, h):
        """
        scaling term for conductance i = gbar*g_s(t, v)

        """
        return m**2*h

    def m_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        d_minf = self.d_minf(v)
        mtau = self.mtau(v)
        d_mtau = self.d_mtau(v)
        minf = self.minf(v)
        return (d_minf/mtau)-(minf-m)/(mtau**2)*d_mtau

    def m_b(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        return 1/self.mtau(v)

    def h_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        d_hinf = self.d_hinf(v)
        htau = self.htau(v)
        d_htau = self.d_htau(v)
        hinf = self.hinf(v)
        return (d_hinf/htau)-(hinf-h)/(htau**2)*d_htau

    def h_b(self, v,m,h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        return 1/self.htau(v)

    def m_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_m

        """
        return 2*(m*h)*(self.E - v)

    def h_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_h

        """
        return (m**2)*(self.E - v)

class K_Pst:
    def __init__(self, v, E = -85.0):
        self.m = self.minf(v)
        self.h = self.hinf(v)
        self.qt = 2.3**((34-21)/10)
        self.E = E

    def update(self, v, dt):
        # self.m = self.m + (self.minf(v)-self.m)/self.mtau(v)*dt
        # self.h = self.h + (self.hinf(v)-self.h)/self.htau(v)*dt
        self.m = self.m + (1-np.exp(-dt/self.mtau(v)))*(self.minf(v)-self.m)
        self.h = self.h + (1-np.exp(-dt/self.htau(v)))*(self.hinf(v)-self.h)

    def mtau(self,v):
        v = v + 10
        if hasattr(v, "__len__"):
            mT = np.zeros(v.shape)
            for k, vi in enumerate(v):
                if hasattr(vi, "__len__"):
                    for kk, vii in enumerate(vi):
                        if vii < -50:
                            mT[k][kk] = (1.25 + 175.03*np.exp(0.026*vii))/self.qt
                        else:
                            mT[k][kk] = (1.25 + 13 * np.exp(-0.026*vii))/self.qt
                else:
                    if vi < -50:
                        mT[k] = (1.25 + 175.03*np.exp(0.026*vi))/self.qt
                    else:
                        mT[k] = (1.25 + 13 * np.exp(-0.026*vi))/self.qt
        else:
            if v < -50:
                mT = (1.25 + 175.03*np.exp(0.026*v))/self.qt
            else:
                mT = (1.25 + 13 * np.exp(-0.026*v))/self.qt

        return mT
    
    def minf(self, v):
        v = v + 10
        return 1/(1 + np.exp(-(v + 1)/12))

    def d_minf(self, v):
        v = v + 10
        return 1/12*np.exp(-(v + 1)/12)/(1 + np.exp(-(v + 1)/12))**2

    def d_mtau(self, v):
        v = v + 10
        if hasattr(v, "__len__"):
            d_mT = np.zeros(v.shape)
            for k, vi in enumerate(v):
                if hasattr(vi, "__len__"):
                    for kk, vii in enumerate(vi):
                        if vii < -50:
                            d_mT[k][kk] = 1/self.qt*175.03*0.026*np.exp(0.026*vii)
                        else:
                            d_mT[k][kk] = 1/self.qt*(-0.026)*13*np.exp(-0.026*vii)
                else:
                    if vi < -50:
                        d_mT[k] = 1/self.qt*175.03*0.026*np.exp(0.026*vi)
                    else:
                        d_mT[k] = 1/self.qt*(-0.026)*13*np.exp(-0.026*vi)
        else:
            if v < -50:
                d_mT = 1/self.qt*175.03*0.026*np.exp(0.026*v)
            else:
                d_mT = 1/self.qt*(-0.026)*13*np.exp(-0.026*v)
        return d_mT

    def htau(self, v):
        v = v + 10
        return (360+(1010+24*(v+55))*np.exp(-((v+75)/48)**2))/self.qt

    def hinf(self, v):
        v = v + 10
        return 1/(1 + np.exp((v+54)/11))

    def d_hinf(self, v):
        v = v + 10
        return -1/11*np.exp((v + 54)/11)/(1 + np.exp((v + 54)/11))**2

    def d_htau(self, v):
        v = v + 10
        return (24*np.exp(-((v + 75)/48)**2) + (24*(v + 55) + 1370)*(np.exp(-((v + 75)/48)**2)*(-2/48**2*(v + 75))))/self.qt

    def g_s(self, m, h):
        """
        scaling term for conductance i = gbar*g_s(t, v)

        """
        return m**2*h

    def m_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        d_minf = self.d_minf(v)
        mtau = self.mtau(v)
        d_mtau = self.d_mtau(v)
        minf = self.minf(v)
        return (d_minf/mtau)-(minf-m)/(mtau**2)*d_mtau

    def m_b(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        return 1/self.mtau(v)

    def h_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        d_hinf = self.d_hinf(v)
        htau = self.htau(v)
        d_htau = self.d_htau(v)
        hinf = self.hinf(v)
        return (d_hinf/htau)-(hinf-h)/(htau**2)*d_htau

    def h_b(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        return 1/self.htau(v)

    def m_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_m

        """
        return 2*(m*h)*(self.E - v)

    def h_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_h

        """
        return (m**2)*(self.E - v)

class K_Tst:
    def __init__(self, v, E = -85.0):
        self.m = self.minf(v)
        self.h = self.hinf(v)
        self.qt = 2.3**((34-21)/10)
        self.E = E

    def update(self, v, dt):
        self.m = self.m + (1-np.exp(-dt/self.mtau(v)))*(self.minf(v)-self.m)
        self.h = self.h + (1-np.exp(-dt/self.htau(v)))*(self.hinf(v)-self.h)

    def mtau(self,v):
        v = v + 10
        return (0.34+0.92*np.exp(-((v+71)/59)**2))/self.qt

    def minf(self, v):
        v = v + 10
        return 1/(1 + np.exp(-(v+0)/19))

    def d_minf(self, v):
        v = v + 10
        return 1/19*np.exp(-(v)/19)/(1 + np.exp(-(v)/19))**2

    def d_mtau(self, v):
        v = v + 10
        return (-2*0.92/(59**2)*np.exp(-((v + 71)/59)**2)*(v + 71))/self.qt

    def htau(self, v):
        v = v + 10
        return (8+49*np.exp(-((v+73)/23)**2))/self.qt

    def hinf(self, v):
        v = v + 10
        return 1/(1 + np.exp((v+66)/10))

    def d_hinf(self, v):
        v = v + 10
        return -1/10*np.exp((v + 66)/10)/(1 + np.exp((v + 66)/10))**2

    def d_htau(self, v):
        v = v + 10
        return (-2*49/(23**2)*np.exp(-((v + 73)/23)**2)*(v + 73))/self.qt

    def g_s(self, m, h):
        """
        scaling term for conductance i = gbar*g_s(t, v)

        """
        return m**4*h

    def m_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        d_minf = self.d_minf(v)
        mtau = self.mtau(v)
        d_mtau = self.d_mtau(v)
        minf = self.minf(v)
        return (d_minf/mtau)-(minf-m)/(mtau**2)*d_mtau

    def m_b(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        return 1/self.mtau(v)

    def h_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        d_hinf = self.d_hinf(v)
        htau = self.htau(v)
        d_htau = self.d_htau(v)
        hinf = self.hinf(v)
        return (d_hinf/htau)-(hinf-h)/(htau**2)*d_htau

    def h_b(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        return 1/self.htau(v)

    def m_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_m

        """
        return 4*(m**3*h)*(self.E - v)

    def h_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_h

        """
        return (m**4)*(self.E - v)

class na:
    def __init__(self, v, E = 50.0):
        self.m = self.minf(v)
        self.h = self.hinf(v)
        self.qt = 2.3**((34-21)/10)
        self.E = E

    def update(self, v, dt):
        # self.m = self.m + (self.minf(v)-self.m)/self.mtau(v)*dt
        # self.h = self.h + (self.hinf(v)-self.h)/self.htau(v)*dt
        self.m = self.m + (1-np.exp(-dt/self.mtau(v)))*(self.minf(v)-self.m)
        self.h = self.h + (1-np.exp(-dt/self.htau(v)))*(self.hinf(v)-self.h)


    def malpha(self, v):
        return (0.182 * (v + 38.))/(1-(np.exp(-(v + 38.)/9.)))

    def d_malpha(self, v):
        return 0.182*((1. - np.exp(-(v + 38.)/9.)) - 1/9.*(v + 38.)*np.exp(-(v + 38.)/9.))/(1. - np.exp(-(v + 38.)/9.))**2

    def mbeta(self, v):
        return (-0.124 * (v + 38.))/(1-(np.exp((v + 38.)/9.)))

    def d_mbeta(self, v):
        return -0.124*((1. - np.exp((v + 38.)/9.)) + 1/9.*(v + 38.)*np.exp((v + 38.)/9.))/(1. - np.exp((v + 38.)/9.))**2

    def mtau(self,v):
        return (1/(self.malpha(v) + self.mbeta(v)))/self.qt
    
    def minf(self, v):
        return self.malpha(v)/(self.malpha(v) + self.mbeta(v))
    
    def halpha(self, v):
        return (-0.024 * (v + 65.))/(1-(np.exp((v + 65.)/6.)))

    def d_halpha(self, v):
        return -0.024*((1. - np.exp((v + 65.)/6.)) + 1/6.*(v + 65.)*np.exp((v + 65.)/6.))/(1. - np.exp((v + 65.)/6.))**2

    def hbeta(self, v):
        return (0.02 * (v + 65.))/(1-(np.exp(-(v + 65.)/6.)))

    def d_hbeta(self, v):

        return 0.02*((1. - np.exp(-(v - (-65.))/6.)) - 1/6.*(v - (-65.))*np.exp(-(v - (-65.))/6.))/(1. - np.exp(-(v - (-65.))/6.))**2

    def htau(self, v):
        return (1/(self.halpha(v) + self.hbeta(v)))/self.qt
    
    def hinf(self, v):
        return self.halpha(v)/(self.halpha(v) + self.hbeta(v))

    def g_s(self, m, h):
        """
        scaling term for conductance i = gbar*g_s(t, v)

        """
        return m**3*h

    def d_minf(self, v):
        ma = self.malpha(v)
        mb = self.mbeta(v)
        d_ma = self.d_malpha(v)
        d_mb = self.d_mbeta(v)
        return (mb*d_ma-ma*d_mb)/((ma+mb)**2)

    def d_mtau(self, v):
        ma = self.malpha(v)
        mb = self.mbeta(v)
        d_ma = self.d_malpha(v)
        d_mb = self.d_mbeta(v)
        return -1/self.qt*(d_ma + d_mb)/((ma + mb)**2)

    def d_hinf(self, v):
        ha = self.halpha(v)
        hb = self.hbeta(v)
        d_ha = self.d_halpha(v)
        d_hb = self.d_hbeta(v)
        return (hb*d_ha-ha*d_hb)/((ha+hb)**2)

    def d_htau(self, v):
        ha = self.halpha(v)
        hb = self.hbeta(v)
        d_ha = self.d_halpha(v)
        d_hb = self.d_hbeta(v)
        return -1/self.qt*(d_ha + d_hb)/((ha + hb)**2)

    def m_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        d_minf = self.d_minf(v)
        mtau = self.mtau(v)
        d_mtau = self.d_mtau(v)
        minf = self.minf(v)
        return (d_minf/mtau)-(minf-m)/(mtau**2)*d_mtau

    def m_b(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        return 1/self.mtau(v)

    def h_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        d_hinf = self.d_hinf(v)
        htau = self.htau(v)
        d_htau = self.d_htau(v)
        hinf = self.hinf(v)
        return (d_hinf/htau)-(hinf-h)/(htau**2)*d_htau

    def h_b(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        return 1/self.htau(v)

    def m_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_m

        """
        return 3*(m**2*h)*(self.E - v)

    def h_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_h

        """
        return (m**3)*(self.E - v)

class kv:
    def __init__(self, v, E = -85.0):
        self.m = self.minf(v)
        self.qt = 2.3**((34-21)/10)
        self.E = E

    def update(self, v, dt):
        # self.m = self.m + (self.minf(v)-self.m)/self.mtau(v)*dt
        # self.h = self.h + (self.hinf(v)-self.h)/self.htau(v)*dt
        self.m = self.m + (1-np.exp(-dt/self.mtau(v)))*(self.minf(v)-self.m)

    def malpha(self, v):
        return 0.02*(v - 25.)/(1. - np.exp(-(v - 25.)/9.))

    def d_malpha(self, v):
    	return 0.02*(1. - np.exp(-(v - 25.)/9.) - 1/9.*(v - 25.)*
    			np.exp(-(v - 25.)/9.))/(1. - np.exp(-(v - 25.)/9.))**2
    
    
    def mbeta(self, v):
    	return -0.002*(v - 25.)/(1. - np.exp((v - 25.)/9.))
    
    
    def d_mbeta(self, v):
    	return -0.002* (1. - np.exp((v - 25.)/9.) + 1/9.*(v - 25.)*
    			np.exp((v - 25.)/9.))/(1. - np.exp((v - 25.)/9.))**2

    def mtau(self,v):
        return (1/(self.malpha(v) + self.mbeta(v)))/self.qt
    
    def minf(self, v):
        return self.malpha(v)/(self.malpha(v) + self.mbeta(v))

    def d_minf(self, v):
        ma = self.malpha(v)
        mb = self.mbeta(v)
        d_ma = self.d_malpha(v)
        d_mb = self.d_mbeta(v)
        return (mb*d_ma-ma*d_mb)/((ma+mb)**2)

    def d_mtau(self, v):
        ma = self.malpha(v)
        mb = self.mbeta(v)
        d_ma = self.d_malpha(v)
        d_mb = self.d_mbeta(v)
        return -1/self.qt*(d_ma + d_mb)/((ma + mb)**2)

    def g_s(self, m):
        """
        scaling term for conductance i = gbar*g_s(t, v)

        """
        return m

    def m_a(self, v, m):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        d_minf = self.d_minf(v)
        mtau = self.mtau(v)
        d_mtau = self.d_mtau(v)
        minf = self.minf(v)
        return (d_minf/mtau)-(minf-m)/(mtau**2)*d_mtau

    def m_b(self, v, m):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        return 1/self.mtau(v)

    def m_c(self, v, m):
        """
        scaling term for partial gradiant computation s_dv = c*s_m

        """
        return (self.E - v)

class kad:
    def __init__(self, v, E = -85.0):
        self.m = self.minf(v)
        self.h = self.hinf(v)
        self.qt = 2.3**((34-21)/10)
        self.E = E

    def update(self, v, dt):
        # self.m = self.m + (self.minf(v)-self.m)/self.mtau(v)*dt
        # self.h = self.h + (self.hinf(v)-self.h)/self.htau(v)*dt
        self.m = self.m + (1-np.exp(-dt/self.mtau(v)))*(self.minf(v)-self.m)
        self.h = self.h + (1-np.exp(-dt/self.htau(v)))*(self.hinf(v)-self.h)

    def malpha(self, v):
        zeta = -1.8 - 1/(1 + np.exp((v + 40)/5))
        return np.exp(1e-3*37.78*zeta*(v + 1))

    def d_malpha(self, v):
        zeta = -1.8 - 1/(1 + np.exp((v + 40)/5))
        return 1e-3*37.78*self.malpha(v)*(zeta + (v + 1)/5*np.exp((v + 40)/5)/(1 + np.exp((v + 40)/5))**2)

    def mbeta(self, v):
        zeta = -1.8 - 1/(1 + np.exp((v + 40)/5))
        return np.exp(1e-3*37.78*0.29*zeta*(v + 1))

    def d_mbeta(self, v):
        zeta = -1.8 - 1/(1 + np.exp((v + 40)/5))
        return 1e-3*37.78*0.29*self.malpha(v)*(zeta + (v + 1)/5*np.exp((v + 40)/5)/(1 + np.exp((v + 40)/5))**2)

    def mtau(self,v):
        return self.mbeta(v)/(self.qt*0.1*(1 + self.malpha(v)))
    
    def minf(self, v):
        return 1/(1 + self.malpha(v))
    
    def halpha(self, v):
        return np.exp(1e-3*37.78*3*(v + 56))

    def d_halpha(self, v):
        return 1e-3*37.78*3*self.halpha(v)


    def htau(self, v):
        return 0.26*(v+50)

    def hinf(self, v):
        return 1/(1 + self.halpha(v)) # self.h_alpha/(self.h_alpha + self.h_beta)

    def g_s(self, m, h):
        """
        scaling term for conductance i = gbar*g_s(t, v)

        """
        return m*h

    def d_minf(self, v):
        return -1/((self.malpha(v)+1)**2)*self.d_malpha(v)

    def d_mtau(self, v):
        ma = self.malpha(v)
        mb = self.mbeta(v)
        d_ma = self.d_malpha(v)
        d_mb = self.d_mbeta(v)
        return 1/(0.1*self.qt)*(d_mb/(1 + ma) - mb*d_ma/((1 + ma)**2))

    def d_hinf(self, v):
        return -1/((self.halpha(v)+1)**2)*self.d_halpha(v)

    def d_htau(self, v):
        return 0.26

    def m_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        d_minf = self.d_minf(v)
        mtau = self.mtau(v)
        d_mtau = self.d_mtau(v)
        minf = self.minf(v)
        return (d_minf/mtau)-(minf-m)/(mtau**2)*d_mtau

    def m_b(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        return 1/self.mtau(v)

    def h_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        d_hinf = self.d_hinf(v)
        htau = self.htau(v)
        d_htau = self.d_htau(v)
        hinf = self.hinf(v)
        return (d_hinf/htau)-(hinf-h)/(htau**2)*d_htau

    def h_b(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        return 1/self.htau(v)

    def m_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_m

        """
        return h*(self.E - v)

    def h_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_h

        """
        return m*(self.E - v)

class kap:
    def __init__(self, v, E = -85.0):
        self.m = self.minf(v)
        self.h = self.hinf(v)
        self.qt = 2.3**((34-21)/10)
        self.E = E

    def update(self, v, dt):
        # self.m = self.m + (self.minf(v)-self.m)/self.mtau(v)*dt
        # self.h = self.h + (self.hinf(v)-self.h)/self.htau(v)*dt
        self.m = self.m + (1-np.exp(-dt/self.mtau(v)))*(self.minf(v)-self.m)
        self.h = self.h + (1-np.exp(-dt/self.htau(v)))*(self.hinf(v)-self.h)

    def malpha(self, v):
        zeta = -1.5 - 1/(1 + np.exp((v + 40)/5))
        return np.exp(1e-3*37.78*zeta*(v - 11))

    def d_malpha(self, v):
        zeta = -1.5 - 1/(1 + np.exp((v + 40)/5))
        return 1e-3*37.78*self.malpha(v)*(zeta + (v - 11)/5*np.exp((v + 40)/5)/(1 + np.exp((v + 40)/5))**2)

    def mbeta(self, v):
        zeta = -1.8 - 1/(1 + np.exp((v + 40)/5))
        return np.exp(1e-3*37.78*0.55*zeta*(v - 11))

    def d_mbeta(self, v):
        zeta = -1.8 - 1/(1 + np.exp((v + 40)/5))
        return 1e-3*37.78*0.55*self.malpha(v)*(zeta + (v - 11)/5*np.exp((v + 40)/5)/(1 + np.exp((v + 40)/5))**2)

    def mtau(self,v):
        return self.mbeta(v)/(self.qt*0.05*(1 + self.malpha(v)))
    
    def minf(self, v):
        return 1/(1 + self.malpha(v))
    
    def halpha(self, v):
        return np.exp(1e-3*37.78*3*(v + 56))

    def d_halpha(self, v):
        return 1e-3*37.78*3*self.halpha(v)


    def htau(self, v):
        return 0.26*(v+50)

    def hinf(self, v):
        return 1/(1 + self.halpha(v)) # self.h_alpha/(self.h_alpha + self.h_beta)

    def g_s(self, m, h):
        """
        scaling term for conductance i = gbar*g_s(t, v)

        """
        return m*h

    def d_minf(self, v):
        return -1/((self.malpha(v)+1)**2)*self.d_malpha(v)

    def d_mtau(self, v):
        ma = self.malpha(v)
        mb = self.mbeta(v)
        d_ma = self.d_malpha(v)
        d_mb = self.d_mbeta(v)
        return 1/(0.05*self.qt)*(d_mb/(1 + ma) - mb*d_ma/((1 + ma)**2))

    def d_hinf(self, v):
        return -1/((self.halpha(v)+1)**2)*self.d_halpha(v)

    def d_htau(self, v):
        return 0.26

    def m_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        d_minf = self.d_minf(v)
        mtau = self.mtau(v)
        d_mtau = self.d_mtau(v)
        minf = self.minf(v)
        return (d_minf/mtau)-(minf-m)/(mtau**2)*d_mtau

    def m_b(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        return 1/self.mtau(v)

    def h_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        d_hinf = self.d_hinf(v)
        htau = self.htau(v)
        d_htau = self.d_htau(v)
        hinf = self.hinf(v)
        return (d_hinf/htau)-(hinf-h)/(htau**2)*d_htau

    def h_b(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        return 1/self.htau(v)

    def m_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_m

        """
        return h*(self.E - v)

    def h_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_h

        """
        return m*(self.E - v)

class ca:
    def __init__(self, v, E = 120.0):
        self.m = self.minf(v)
        self.h = self.hinf(v)
        self.qt = 2.3**((34-21)/10)
        self.E = E

    def update(self, v, dt):
        # self.m = self.m + (self.minf(v)-self.m)/self.mtau(v)*dt
        # self.h = self.h + (self.hinf(v)-self.h)/self.htau(v)*dt
        self.m = self.m + (1-np.exp(-dt/self.mtau(v)))*(self.minf(v)-self.m)
        self.h = self.h + (1-np.exp(-dt/self.htau(v)))*(self.hinf(v)-self.h)


    def malpha(self, v):
        return (0.055 * (v + 27.))/(1-(np.exp(-(v + 27.)/3.8)))

    def d_malpha(self, v):
        return 0.055*((1. - np.exp(-(v + 27.)/3.8)) - 1/3.8*(v + 27.)*np.exp(-(v + 27.)/3.8))/(1. - np.exp(-(v + 27.)/3.8))**2

    def mbeta(self, v):
        return -0.94*np.exp(-(v + 75)/17)

    def d_mbeta(self, v):
        return -0.94/27*np.exp(-(v + 75)/17)

    def mtau(self,v):
        return (1/(self.malpha(v) + self.mbeta(v)))/self.qt
    
    def minf(self, v):
        return self.malpha(v)/(self.malpha(v) + self.mbeta(v))
    
    def halpha(self, v):
        return 0.000457*np.exp(-(v + 13)/50)

    def d_halpha(self, v):
        return -0.000457/50*np.exp(-(v + 13)/50)

    def hbeta(self, v):
        return 0.0065/(1 + np.exp(-(v + 15)/28))

    def d_hbeta(self, v):
        return 0.0065/28*np.exp(-(v + 15)/28)/(1 + np.exp(-(v + 15)/28))**2

    def htau(self, v):
        return (1/(self.halpha(v) + self.hbeta(v)))/self.qt
    
    def hinf(self, v):
        return self.halpha(v)/(self.halpha(v) + self.hbeta(v))

    def g_s(self, m, h):
        """
        scaling term for conductance i = gbar*g_s(t, v)

        """
        return m**2*h

    def d_minf(self, v):
        ma = self.malpha(v)
        mb = self.mbeta(v)
        d_ma = self.d_malpha(v)
        d_mb = self.d_mbeta(v)
        return (mb*d_ma-ma*d_mb)/((ma+mb)**2)

    def d_mtau(self, v):
        ma = self.malpha(v)
        mb = self.mbeta(v)
        d_ma = self.d_malpha(v)
        d_mb = self.d_mbeta(v)
        return -1/self.qt*(d_ma + d_mb)/((ma + mb)**2)

    def d_hinf(self, v):
        ha = self.halpha(v)
        hb = self.hbeta(v)
        d_ha = self.d_halpha(v)
        d_hb = self.d_hbeta(v)
        return (hb*d_ha-ha*d_hb)/((ha+hb)**2)

    def d_htau(self, v):
        ha = self.halpha(v)
        hb = self.hbeta(v)
        d_ha = self.d_halpha(v)
        d_hb = self.d_hbeta(v)
        return -1/self.qt*(d_ha + d_hb)/((ha + hb)**2)

    def m_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        d_minf = self.d_minf(v)
        mtau = self.mtau(v)
        d_mtau = self.d_mtau(v)
        minf = self.minf(v)
        return (d_minf/mtau)-(minf-m)/(mtau**2)*d_mtau

    def m_b(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        return 1/self.mtau(v)

    def h_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        d_hinf = self.d_hinf(v)
        htau = self.htau(v)
        d_htau = self.d_htau(v)
        hinf = self.hinf(v)
        return (d_hinf/htau)-(hinf-h)/(htau**2)*d_htau

    def h_b(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        return 1/self.htau(v)

    def m_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_m

        """
        return 2*(m*h)*(self.E - v)

    def h_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_h

        """
        return (m**2)*(self.E - v)

class it:
    def __init__(self, v, E = 120.0):
        self.m = self.minf(v)
        self.h = self.hinf(v)
        self.qt = 2.3**((34-21)/10)
        self.E = E

    def update(self, v, dt):
        self.m = self.m + (1-np.exp(-dt/self.mtau(v)))*(self.minf(v)-self.m)
        self.h = self.h + (1-np.exp(-dt/self.htau(v)))*(self.hinf(v)-self.h)

    def minf(self, v):
        return 1/(1 + np.exp((v + 50)/(-7.4))) # self.m_alpha/(self.m_alpha + self.m_beta)

    def mtau(self,v):
        return 3 + 1/(np.exp((v + 25)/20) + np.exp(-(v + 100)/15))

    def d_minf(self, v):
        return 1/7.4*np.exp(-(v + 50)/7.4)/((1 + np.exp(-(v + 50)/7.4))**2)

    def d_mtau(self, v):
        return -(np.exp((v + 25)/20) - np.exp(-(v + 100)/15))/((np.exp((v + 25)/20) + np.exp(-(v + 100)/15))**2)

    def hinf(self, v):
        return 1/(1 + np.exp((v + 78)/5.0)) # self.m_alpha/(self.m_alpha + self.m_beta)

    def htau(self,v):
        return 85 + 1/(np.exp((v + 46)/4) + np.exp(-(v + 405)/50))

    def d_hinf(self, v):
        return -1/5*np.exp((v + 78)/5)/((1 + np.exp((v + 78)/5))**2)

    def d_htau(self, v):
        return -(np.exp((v + 46)/4) - np.exp(-(v + 405)/50))/((np.exp((v + 46)/4) + np.exp(-(v + 405)/50))**2)

    def g_s(self, m, h):
        """
        scaling term for conductance i = gbar*g_s(t, v)

        """
        return m**2*h

    def m_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        d_minf = self.d_minf(v)
        mtau = self.mtau(v)
        d_mtau = self.d_mtau(v)
        minf = self.minf(v)
        return (d_minf/mtau)-(minf-m)/(mtau**2)*d_mtau

    def m_b(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        return 1/self.mtau(v)

    def h_a(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        d_hinf = self.d_hinf(v)
        htau = self.htau(v)
        d_htau = self.d_htau(v)
        hinf = self.hinf(v)
        return (d_hinf/htau)-(hinf-h)/(htau**2)*d_htau

    def h_b(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dh = a*s_v - b*s_h

        """
        return 1/self.htau(v)

    def m_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_m

        """
        return 2*(m*h)*(self.E - v)

    def h_c(self, v, m, h):
        """
        scaling term for partial gradiant computation s_dv = c*s_h

        """
        return (m**2)*(self.E - v)

class SK_E2:
    def __init__(self, ca =  1e-4, E = -85.0):
        self.m = self.minf(ca)
        self.qt = 2.3**((34-21)/10)
        self.E = E

    def update(self, ca, dt):
        # self.m = self.m + (self.minf(v)-self.m)/self.mtau(v)*dt
        # self.h = self.h + (self.hinf(v)-self.h)/self.htau(v)*dt
        self.m = self.m + (1-np.exp(-dt/self.mtau(ca)))*(self.minf(ca)-self.m)

    def mtau(self,ca):
        return 1
    
    def minf(self, ca):
        if hasattr(ca, '__len__'):
            for k, ca_i in enumerate(ca):
                if ca_i <= 1e-7:
                    ca[k] = 1e-7
        else:
            if ca <= 1e-7:
                ca = 1e-7
        return 1/(1 + (0.00043/ca)**4.8)

    def g_s(self, m):
        """
        scaling term for conductance i = gbar*g_s(t, v)

        """
        return m

    def d_minf(self, ca, dca):
        return 0.00043*4.8*(0.00043/ca)**3.8*dca/(ca**2)*(self.minf(ca))**2

    def d_mtau(self, ca, dca):
        return 0

    def m_a(self, ca, dca):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        d_minf = self.d_minf(ca, dca)
        mtau = self.mtau(ca)
        return d_minf/mtau

    def m_b(self, ca = 0, dca = 0):
        """
        scaling term for partial gradiant computation s_dm = a*s_v - b*s_m

        """
        return 1

    def m_c(self, v, m):
        """
        scaling term for partial gradiant computation s_dv = c*s_m

        """
        return self.E - v


class kBK:
    def __init__(self, v, ca =  1e-4, E = -85.0):
        self.m = self.minf(v, ca)
        self.qt = 2.3**((34-21)/10)
        self.E = E

    def update(self, ca, v, dt):
        # self.m = self.m + (self.minf(v)-self.m)/self.mtau(v)*dt
        # self.h = self.h + (self.hinf(v)-self.h)/self.htau(v)*dt
        self.m = self.m + (1-np.exp(-dt/self.mtau(ca)))*(self.minf(ca)-self.m)

    def mtau(self, v, ca):
        return 1
    
    def minf(self, v, ca):
        ca = ca * 100
        if hasattr(ca, '__len__'):
            for k, ca_i in enumerate(ca):
                if ca_i <= 1e-7:
                    ca[k] = 1e-7
        else:
            if ca <= 1e-7:
                ca = 1e-7
        P0ca = 1/ ( 1 + (2e-3/ca) )
        Vhca = -46.08 + ( (155.67 + 46.08 ) / ( 1 + ((2e-3/ca)**(-0.94208))) )
        return P0ca/ ( 1 + np.exp( (Vhca-v)/17) )

    def g_s(self, m):
        """
        scaling term for conductance i = gbar*g_s(t, v)

        """
        return m

class CaDynamics_E2:
    def __init__(self, ica = 0):
        self.minCai = 1e-4
        self.depth = 0.1
        self.ca = self.minCai

    def update(self, ica, ca, gamma, decay, dt):
        ca = ca + dt*(-10000*ica*gamma/(2*96500*self.depth) - (ca - self.minCai)/decay)
        return ca


