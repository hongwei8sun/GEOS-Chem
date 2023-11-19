#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import xarray as xr

from tqdm import tqdm
import os

import scipy.constants
AVO = scipy.constants.Avogadro

# In[2]:


GF_THRESHOLD = 0.0e+0
GF_RANGE     = 1.0e-8
GF_DELTAHBYR = 10156.0e+0
GF_T0        = 360.0e+0
GF_TC        = 905.0e+0

#ATM : Standard atmosphere [Pa]  (Source: NIST, 2014)
ATM     = 1.01325e+5
#AIRMW : Average molecular weight of dry air [g/mol]
AIRMW    = 28.97e+0



folder = '/n/home12/hongwei/HONGWEI/GC_aerosol/rundirs/Volcano/Add_H2SO4/geosfp_4x5_gc_timing_Seb_Sigma0.5/'
file = 'GEOSChem.Restart.NO_H2SO4G.20150101_0000z.nc4'
df = xr.open_dataset(folder+file)



# In[5]:


lon = df['lon']
lat = df['lat']
lev = df['lev']

Nx = len(lon)
Ny = len(lat)
Nz = len(lev)


# In[6]:


# PCENTER_P: dry pressure

hyam = df['hyam']
hybm = df['hybm']
P0 = df['P0']
P = (hyam + hybm * P0) # [hPa]


# In[7]:


# PCENTER = PCENTER_P + e
# moist air pressure = dry air pressure + vapor pressure
# e: partial water vapor pressure

q = df["Met_SPHU1"] # g/kg
q = q/1000.0 # kg/kg

e = q*0.0

for iz in range(Nz):
    # P: pressure [hPa]
    e[0,iz,:,:] = q[0,iz,:,:]*P[iz] / (0.622 + 0.378*q[0,iz,:,:])
    

# In[8]:


# water vapor mixing ratio [mol/mol]

SpeciesRst_H2O = (q/18.0) / ((1-q)/28.96) # [kg/kg] => [mol/mol]

SpeciesRst_H2O.shape


# In[9]:


# AD: dry air mass [kg] = volume * density = area * dz * density = area * abs(-dp)/g

Area = df['AREA'] # [m2]
dP = df['Met_DELPDRY']*100.0 # [Pa]
g = 9.8 # [m/s2]

AD = dP*0.0
for iz in range(Nz):
    AD[0,iz,:,:] = Area[:,:] * dP[0,iz,:,:] / g
    

# In[10]:


# TCENTER: temperature
T = df['Met_TMPU1'] # [K]


# In[11]:


# Spc_SO4: 

Spc_SO4 = df['SpeciesRst_SO4'] # mol mol-1

Spc_SO4.shape


# ## InTroposphere is set to be 0 now, which need to be update in the future!!!

# In[12]:


SO4_MW_G = 96.0 # g/mol
H2O_MW_G = 18.0

# Calculate H2SO4 gas phase prefactors
GF_INVT0 = 1.0/GF_T0
GF_LOGP0 = -1.0*GF_DELTAHBYR*GF_INVT0 + 16.259
GF_BFACTOR = 0.38/(GF_TC - GF_T0)
GF_ATMCONV = np.log(ATM)


InTroposphere = False

TCENTER = T
PCENTER_P = e*0.0
# P
PCENTER = e*0.0
# PCENTER_P + e
        
for iz in range(Nz):
    PCENTER_P[0,iz,:,:] = P[iz]
    PCENTER[0,iz,:,:] = PCENTER_P[0,iz,:,:] + e[0,iz,:,:]
    
PCENTER.shape, PCENTER_P.shape, e.shape, P.shape


# In[13]:


# calcualte the gas-liquid partition for H2SO4:
    # TCENTER = State_Met%T(I,J,L).
    # PCENtER = State_Met%PMID(I,J,L): moist air pressure.
    # PCENTER_P = State_Met%PMID_DRY(I,J,L): Dry air partial pressure for part. P calc.
    # AD = State_Met%AD(I,J,L), ! Dry air mass [kg] in grid box
    # InTroposphere = State_Met%InTroposphere(I,J,L)
    # Spc_SO4 = Spc(I,J,L,id_SO4), [mol/mol]
    
    # SO4_MW_G = State_Chm%SpcData(id_SO4)%Info%emMW_g ! g/mol
    

# INVAIR = AIRMW / AD
# print(INVAIR.shape)

# H2SO4SUM = Spc_SO4 * INVAIR / SO4_MW_G
H2SO4SUM = Spc_SO4
print(H2SO4SUM.shape)

GF_PP = H2SO4SUM * PCENTER_P
print(GF_PP.shape)

GF_INVT = 1.0/TCENTER
print(GF_INVT.shape)

GF_CFACTOR = 1.0 + np.log(GF_T0*GF_INVT) - GF_T0*GF_INVT
print(GF_CFACTOR.shape)

GF_LOGPSULFATE = GF_LOGP0 + GF_DELTAHBYR * (GF_INVT0 - GF_INVT + GF_BFACTOR*GF_CFACTOR)
print(GF_LOGPSULFATE.shape)

GF_LOGPSULFATE = GF_LOGPSULFATE + GF_ATMCONV
print(GF_LOGPSULFATE.shape)

GF_PVAP = 1.0e-2 * np.exp(GF_LOGPSULFATE)
print(GF_PVAP.shape)

GF_DIFF = (GF_PVAP+GF_THRESHOLD) - GF_PP
print(GF_DIFF.shape)


AERFRAC = GF_DIFF*0.0 + 1.0

for ix in tqdm(range(Nx)):
    for iy in range(Ny):
        for iz in range(Nz):
            
            if PCENTER[0,iz,iy,ix]>=100.0:
                AERFRAC[0,iz,iy,ix] = 1.0
                
            elif InTroposphere:
                AERFRAC[0,iz,iy,ix] = 1.0
                
            else:
                if GF_DIFF[0,iz,iy,ix] < 0.0:
                    AERFRAC[0,iz,iy,ix] = 1.0
                elif GF_DIFF[0,iz,iy,ix] < GF_RANGE:
                    AERFRAC[0,iz,iy,ix] = 1.0 - (GF_DIFF[0,iz,iy,ix]/GF_RANGE)
                else:
                    AERFRAC[0,iz,iy,ix] = 0.0


# In[14]:


def func(TCENTER, PCENTER, PCENTER_P, AD, InTroposphere, Spc_SO4):
    # TCENTER = State_Met%T(I,J,L).
    # PCENtER = State_Met%PMID(I,J,L): moist air pressure.
    # PCENTER_P = State_Met%PMID_DRY(I,J,L): Dry air partial pressure for part. P calc.
    # AD = State_Met%AD(I,J,L), ! Dry air mass [kg] in grid box
    # InTroposphere = State_Met%InTroposphere(I,J,L)
    # Spc_SO4 = Spc(I,J,L,id_SO4), [mol/mol]
    
    # SO4_MW_G = State_Chm%SpcData(id_SO4)%Info%emMW_g ! g/mol
    
    # INVAIR = AIRMW / AD
            
    if PCENTER>=100.0:
        AERFRAC = 1.0
    elif InTroposphere:
        AERFRAC = 1.0
    else:
        # H2SO4SUM = Spc_SO4 * INVAIR / SO4_MW_G
        H2SO4SUM = Spc_SO4
        
        GF_PP = H2SO4SUM * PCENTER_P
        GF_INVT = 1.0/TCENTER
        GF_CFACTOR = 1.0 + np.log(GF_T0*GF_INVT) - GF_T0*GF_INVT
        
        GF_LOGPSULFATE = GF_LOGP0 + GF_DELTAHBYR * (GF_INVT0 - GF_INVT + GF_BFACTOR*GF_CFACTOR)
        
        GF_LOGPSULFATE = GF_LOGPSULFATE + GF_ATMCONV
        
        GF_PVAP = 1.0e-2 * np.exp(GF_LOGPSULFATE)
        GF_DIFF = (GF_PVAP+GF_THRESHOLD) - GF_PP
        
        if GF_DIFF < 0.0:
            AERFRAC = 1.0
        elif GF_DIFF < GF_RANGE:
            AERFRAC = 1.0 - (GF_DIFF/GF_RANGE)
        else:
            AERFRAC = 0.0


# ## Partition H2SO4 into Gas and Liquid phases


Spc_H2SO4G = Spc_SO4*0.0
Spc_H2SO4G = Spc_SO4 * (1-AERFRAC)
Spc_SO4_new = Spc_SO4 * AERFRAC

AERFRAC.shape, Spc_SO4.shape


# # Asign liquid SO4 to 40-bin

def CARSLAW_DENSITY(M_H2SO4,M_HNO3,TCENTER):
    
    CS = M_H2SO4
    CN = M_HNO3
    T = TCENTER
    
    DENSS=1000.0 +123.64 *CS-5.6e-4 *CS *T**2 -29.54 *CS**1.5 +1.814e-4 *CS**1.5 *T**2 +2.343 *CS**2 -1.487e-3 *CS**2*T -1.324e-5 *CS**2 *T**2

    DENSN=1000.0e+0 +85.107 *CN-5.043e-4 *CN*T**2 -18.96 *CN**1.5 +1.427e-4 * CN**1.5*T**2+1.458 *CN**2-1.198e-3 *CN**2 *T -9.703e-6 *CN**2 *T**2

    # Error trap for zeros (ckeller, 12/29/17)
    if( CS == 0.0 and CN == 0.0 ):
        SLA_rho = 0.0
    else:
        SLA_rho =1.0e+0 / ((1.0e+0/DENSS*CS/(CS+CN) +1.0e+0 / DENSN*CN / (CS+CN) ))
    
    return SLA_rho


# In[22]:


# ATM : Standard atmosphere [Pa]  (Source: NIST, 2014)
# RSTARG : Molar gas constant [J/K/mol]

KS = [0.0, -21.661e+0, 2724.2e+0, 51.81e+0, -15732.0e+0, 47.004e+0, -6969.0e+0 ,-4.6183e+0]
ATM = 1.01325e+5

M_HNO3 = 0.0e+0
TNY = 1.0e-28

RSTARG = 8.3144598e+0
M_HNO3 = 0.0e+0
    
def TERNARY(PCENTER_in, TCENTER_in, H2OSUM_in, H2SO4SUM_in):
    
    PCENTER = max(PCENTER_in,5.0e+0)
    TCENTER = TCENTER_in
    H2OSUM_IN = max(H2OSUM_in,5.0e-7)
    H2SO4SUM_IN = H2SO4SUM_in
    
    # Calculate partial pressure of H2O & HNO3
    # PCENTER is in hPa, so need to convert ATM from Pa to hPa
    PATMH2O  = H2OSUM_IN  * PCENTER / (ATM*1e-2)

    # Carslaw only valid for 2e-5 < PPH2O < 2e-3 (hPa)
    PATMH2O = max(PATMH2O,1.9738465e-8)
    PATMH2O = min(PATMH2O,1.9738465e-6)
    
    # Determine H2SO4/H2O pure solution concentration
    # Mole fraction of H2SO4 in binary solution
    TMP1 = (KS[1]+KS[2]/TCENTER)**2.0e+0-4.0e+0*(KS[3]+KS[4]/TCENTER)*(KS[5]+KS[6]/TCENTER+KS[7]*np.log(TCENTER)-np.log(PATMH2O))
    if( TMP1 > 0.0 ):
        X_H2SO4_BIN = 1.0e+0/(2.0e+0*(KS[3]+KS[4]/TCENTER))*(-KS[1]-KS[2]/TCENTER-(TMP1)**0.5e+0)
    else:
        X_H2SO4_BIN = 0.0
    
    # Molality (mol H2SO4/kg H2O) in binary solution
    M_H2SO4_BIN = 55.51e+0*X_H2SO4_BIN/(1.0e+0-X_H2SO4_BIN)
    M_H2SO4 = M_H2SO4_BIN
    W_H2SO4 = M_H2SO4_BIN*0.098076e+0/(1.0e+0+M_H2SO4_BIN*0.098076e+0)

    SLA_RHO = CARSLAW_DENSITY(M_H2SO4,M_HNO3,TCENTER)
    
    # Moles of H2SO4 per m3 air
    MOLDENS_H2SO4 = 100.0e+0*PCENTER*H2SO4SUM_IN/(RSTARG*TCENTER)
    
    # Aerosol mass density in kg/m3 aerosol
    SLA_RHO = CARSLAW_DENSITY(M_H2SO4,M_HNO3,TCENTER)
    
    # Aerosol volume in m3/m3 air
    if( W_H2SO4 < TNY or SLA_RHO < TNY ):
        SLA_VOL = 0.0
    else:
        SLA_VOL = (MOLDENS_H2SO4*98.076/W_H2SO4/SLA_RHO)*1.0e-3
        
    return SLA_VOL


# In[17]:


# SLA_VR     : SLA volume-effective radius conversion
# SLA_RR     : SLA effective-liquid radius conversion
# SLA_VOL    : Aerosol volume (m3/m3)

SLA_VR = (0.357e-6)*(10.0e+0**(12.0e+0*0.249))
SLA_RR = np.exp(-0.173)


# In[18]:


n_aer_bin = 40
den_h2so4=1.8E-12 # pure h2so4 density in g/um^3, used for calculating h2so4 mass/particle !eth_af_dryS
aer_R0=3.9376E-4  # smallest bin's dry sulfate radius in um
aer_pi = np.pi
aer_Vrat = 2.0e+0

aer_mass = np.zeros(n_aer_bin)
aer_dry_rad = np.zeros(n_aer_bin)

for k in range(n_aer_bin):
    if(k==0):
        aer_mass[k]=den_h2so4*4./3.*aer_pi*aer_R0**3 # mass H2SO4/particle in g
        aer_dry_rad[k] = aer_R0
    else:
        aer_mass[k]=aer_mass[k-1]*aer_Vrat # mass H2SO4/particle in g
        aer_dry_rad[k] = (3.0*aer_mass[k]/(4.0*aer_pi*den_h2so4))**(1.0/3.0)


# In[19]:


Spc_Bin_0 = Spc_SO4*0.0
Spc_Bin_1 = Spc_SO4*0.0
Spc_Bin_2 = Spc_SO4*0.0
Spc_Bin_3 = Spc_SO4*0.0
Spc_Bin_4 = Spc_SO4*0.0
Spc_Bin_5 = Spc_SO4*0.0
Spc_Bin_6 = Spc_SO4*0.0
Spc_Bin_7 = Spc_SO4*0.0
Spc_Bin_8 = Spc_SO4*0.0
Spc_Bin_9 = Spc_SO4*0.0

Spc_Bin_10 = Spc_SO4*0.0
Spc_Bin_11 = Spc_SO4*0.0
Spc_Bin_12 = Spc_SO4*0.0
Spc_Bin_13 = Spc_SO4*0.0
Spc_Bin_14 = Spc_SO4*0.0
Spc_Bin_15 = Spc_SO4*0.0
Spc_Bin_16 = Spc_SO4*0.0
Spc_Bin_17 = Spc_SO4*0.0
Spc_Bin_18 = Spc_SO4*0.0
Spc_Bin_19 = Spc_SO4*0.0

Spc_Bin_20 = Spc_SO4*0.0
Spc_Bin_21 = Spc_SO4*0.0
Spc_Bin_22 = Spc_SO4*0.0
Spc_Bin_23 = Spc_SO4*0.0
Spc_Bin_24 = Spc_SO4*0.0
Spc_Bin_25 = Spc_SO4*0.0
Spc_Bin_26 = Spc_SO4*0.0
Spc_Bin_27 = Spc_SO4*0.0
Spc_Bin_28 = Spc_SO4*0.0
Spc_Bin_29 = Spc_SO4*0.0

Spc_Bin_30 = Spc_SO4*0.0
Spc_Bin_31 = Spc_SO4*0.0
Spc_Bin_32 = Spc_SO4*0.0
Spc_Bin_33 = Spc_SO4*0.0
Spc_Bin_34 = Spc_SO4*0.0
Spc_Bin_35 = Spc_SO4*0.0
Spc_Bin_36 = Spc_SO4*0.0
Spc_Bin_37 = Spc_SO4*0.0
Spc_Bin_38 = Spc_SO4*0.0
Spc_Bin_39 = Spc_SO4*0.0


# In[ ]:


# H2OSUM = SpeciesRst_H2O * INVAIR / H2O_MW_G
# H2SO4SUM = Spc_SO4_new * INVAIR / SO4_MW_G
H2OSUM = SpeciesRst_H2O
H2SO4SUM = Spc_SO4_new

Sigma = 0.55 # 1.82

Wts = np.zeros(n_aer_bin)

for ix in tqdm(range(Nx)):
    for iy in range(Ny):
        for iz in range(Nz):
            
            # (1) calculate the effective radius
            PCENTER_in = PCENTER[0,iz,iy,ix]
            TCENTER_in = TCENTER[0,iz,iy,ix]
            H2OSUM_in = H2OSUM[0,iz,iy,ix]
            H2SO4SUM_in = H2SO4SUM[0,iz,iy,ix]
            
            if H2SO4SUM_in>0.0 and PCENTER_in<100.0:
                
                VOL_SLA = TERNARY(PCENTER_in, TCENTER_in, H2OSUM_in, H2SO4SUM_in)
            
                RAD_AER_BOX = SLA_VR*SLA_RR*(VOL_SLA**0.249e+0) # m
                
                # (2) set initial Rm and sigma
                Rm = RAD_AER_BOX*1e6/np.exp(5.0/2.0*Sigma**2) # [um]
                
                # (3) set initial unimodal lognormal distribution
                for I_Bin in range(40):
                    r  = aer_dry_rad[I_Bin] # [um]
                    aer_vol_dry = 4/3 *np.pi * r**3

                    if(I_Bin==0):
                        Dr = aer_dry_rad[I_Bin+1] - aer_dry_rad[I_Bin]
                    elif(I_Bin==39):
                        Dr = aer_dry_rad[I_Bin] - aer_dry_rad[I_Bin-1]
                    else:
                        Dr = (aer_dry_rad[I_Bin+1] - aer_dry_rad[I_Bin-1])/2.0

                    # total particle mass = mass per particle * particle number density
                    # Assume density is same in different bins, so mass distribution is
                    # same as volume distribution

                    # Wts[I_Bin] = aer_mass[I_Bin] * Dr * 1/(Sigma*np.sqrt(2.0*np.pi)) * 1/r * np.exp( -1* (np.log(r)-np.log(Rm))**2 / (2.0*Sigma**2) )
                    Wts[I_Bin] = (den_h2so4*aer_vol_dry/98.076*AVO) *Dr *1/(Sigma*np.sqrt(2.0*np.pi)) *1/r *np.exp( -1*(np.log(r)-np.log(Rm))**2 / (2.0*Sigma**2) )
            
                Spc_Bin_0[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[0]/np.sum(Wts)
                Spc_Bin_1[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[1]/np.sum(Wts)
                Spc_Bin_2[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[2]/np.sum(Wts)
                Spc_Bin_3[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[3]/np.sum(Wts)
                Spc_Bin_4[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[4]/np.sum(Wts)
                Spc_Bin_5[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[5]/np.sum(Wts)
                Spc_Bin_6[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[6]/np.sum(Wts)
                Spc_Bin_7[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[7]/np.sum(Wts)
                Spc_Bin_8[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[8]/np.sum(Wts)
                Spc_Bin_9[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[9]/np.sum(Wts)
            
                Spc_Bin_10[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[10]/np.sum(Wts)
                Spc_Bin_11[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[11]/np.sum(Wts)
                Spc_Bin_12[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[12]/np.sum(Wts)
                Spc_Bin_13[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[13]/np.sum(Wts)
                Spc_Bin_14[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[14]/np.sum(Wts)
                Spc_Bin_15[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[15]/np.sum(Wts)
                Spc_Bin_16[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[16]/np.sum(Wts)
                Spc_Bin_17[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[17]/np.sum(Wts)
                Spc_Bin_18[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[18]/np.sum(Wts)
                Spc_Bin_19[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[19]/np.sum(Wts)  
            
                Spc_Bin_20[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[20]/np.sum(Wts)
                Spc_Bin_21[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[21]/np.sum(Wts)
                Spc_Bin_22[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[22]/np.sum(Wts)
                Spc_Bin_23[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[23]/np.sum(Wts)
                Spc_Bin_24[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[24]/np.sum(Wts)
                Spc_Bin_25[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[25]/np.sum(Wts)
                Spc_Bin_26[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[26]/np.sum(Wts)
                Spc_Bin_27[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[27]/np.sum(Wts)
                Spc_Bin_28[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[28]/np.sum(Wts)
                Spc_Bin_29[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[29]/np.sum(Wts)
            
                Spc_Bin_30[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[30]/np.sum(Wts)
                Spc_Bin_31[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[31]/np.sum(Wts)
                Spc_Bin_32[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[32]/np.sum(Wts)
                Spc_Bin_33[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[33]/np.sum(Wts)
                Spc_Bin_34[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[34]/np.sum(Wts)
                Spc_Bin_35[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[35]/np.sum(Wts)
                Spc_Bin_36[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[36]/np.sum(Wts)
                Spc_Bin_37[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[37]/np.sum(Wts)
                Spc_Bin_38[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[38]/np.sum(Wts)
                Spc_Bin_39[0,iz,iy,ix] = Spc_SO4_new[0,iz,iy,ix] * Wts[39]/np.sum(Wts)
            
                Spc_SO4_new[0,iz,iy,ix] = 0.0
            

# # Add new variables to NetCDF file


Spc_H2SO4G.attrs['long_name'] = 'Dry mixing ratio of species H2SO4G'
Spc_H2SO4G.attrs['units'] = 'mol mol-1 dry'
Spc_H2SO4G.attrs['averaging_method'] = 'instantaneous'

Spc_SO4_new.attrs['long_name'] = 'Dry mixing ratio of species SO4'
Spc_SO4_new.attrs['units'] = 'mol mol-1 dry'
Spc_SO4_new.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_0.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 1'
Spc_Bin_0.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_0.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_1.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 2'
Spc_Bin_1.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_1.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_2.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 3'
Spc_Bin_2.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_2.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_3.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 4'
Spc_Bin_3.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_3.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_4.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 5'
Spc_Bin_4.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_4.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_5.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 6'
Spc_Bin_5.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_5.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_6.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 7'
Spc_Bin_6.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_6.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_7.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 8'
Spc_Bin_7.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_7.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_8.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 9'
Spc_Bin_8.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_8.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_9.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 10'
Spc_Bin_9.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_9.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_10.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 11'
Spc_Bin_10.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_10.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_11.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 12'
Spc_Bin_11.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_11.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_12.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 13'
Spc_Bin_12.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_12.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_13.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 14'
Spc_Bin_13.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_13.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_14.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 15'
Spc_Bin_14.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_14.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_15.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 16'
Spc_Bin_15.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_15.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_16.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 17'
Spc_Bin_16.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_16.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_17.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 18'
Spc_Bin_17.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_17.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_18.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 19'
Spc_Bin_18.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_18.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_19.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 20'
Spc_Bin_19.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_19.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_20.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 21'
Spc_Bin_20.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_20.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_21.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 22'
Spc_Bin_21.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_21.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_22.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 23'
Spc_Bin_22.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_22.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_23.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 24'
Spc_Bin_23.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_23.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_24.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 25'
Spc_Bin_24.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_24.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_25.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 26'
Spc_Bin_25.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_25.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_26.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 27'
Spc_Bin_26.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_26.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_27.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 28'
Spc_Bin_27.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_27.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_28.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 29'
Spc_Bin_28.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_28.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_29.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 30'
Spc_Bin_29.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_29.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_30.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 31'
Spc_Bin_30.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_30.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_31.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 32'
Spc_Bin_31.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_31.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_32.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 33'
Spc_Bin_32.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_32.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_33.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 34'
Spc_Bin_33.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_33.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_34.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 35'
Spc_Bin_34.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_34.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_35.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 36'
Spc_Bin_35.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_35.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_36.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 37'
Spc_Bin_36.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_36.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_37.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 38'
Spc_Bin_37.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_37.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_38.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 39'
Spc_Bin_38.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_38.attrs['averaging_method'] = 'instantaneous'

Spc_Bin_39.attrs['long_name'] = 'Dry mixing ratio of species H2SO4 in Bin 40'
Spc_Bin_39.attrs['units'] = 'mol mol-1 dry'
Spc_Bin_39.attrs['averaging_method'] = 'instantaneous'

Spc_H2SO4G.shape


# In[ ]:


df2 = df
df2['SpeciesRst_H2SO4G'] = Spc_H2SO4G
df2['SpeciesRst_SO4'] = Spc_SO4_new

df2['SpeciesRst_AERSctSul001'] = Spc_Bin_0
df2['SpeciesRst_AERSctSul002'] = Spc_Bin_1
df2['SpeciesRst_AERSctSul003'] = Spc_Bin_2
df2['SpeciesRst_AERSctSul004'] = Spc_Bin_3
df2['SpeciesRst_AERSctSul005'] = Spc_Bin_4
df2['SpeciesRst_AERSctSul006'] = Spc_Bin_5
df2['SpeciesRst_AERSctSul007'] = Spc_Bin_6
df2['SpeciesRst_AERSctSul008'] = Spc_Bin_7
df2['SpeciesRst_AERSctSul009'] = Spc_Bin_8
df2['SpeciesRst_AERSctSul010'] = Spc_Bin_9

df2['SpeciesRst_AERSctSul011'] = Spc_Bin_10
df2['SpeciesRst_AERSctSul012'] = Spc_Bin_11
df2['SpeciesRst_AERSctSul013'] = Spc_Bin_12
df2['SpeciesRst_AERSctSul014'] = Spc_Bin_13
df2['SpeciesRst_AERSctSul015'] = Spc_Bin_14
df2['SpeciesRst_AERSctSul016'] = Spc_Bin_15
df2['SpeciesRst_AERSctSul017'] = Spc_Bin_16
df2['SpeciesRst_AERSctSul018'] = Spc_Bin_17
df2['SpeciesRst_AERSctSul019'] = Spc_Bin_18
df2['SpeciesRst_AERSctSul020'] = Spc_Bin_19

df2['SpeciesRst_AERSctSul021'] = Spc_Bin_20
df2['SpeciesRst_AERSctSul022'] = Spc_Bin_21
df2['SpeciesRst_AERSctSul023'] = Spc_Bin_22
df2['SpeciesRst_AERSctSul024'] = Spc_Bin_23
df2['SpeciesRst_AERSctSul025'] = Spc_Bin_24
df2['SpeciesRst_AERSctSul026'] = Spc_Bin_25
df2['SpeciesRst_AERSctSul027'] = Spc_Bin_26
df2['SpeciesRst_AERSctSul028'] = Spc_Bin_27
df2['SpeciesRst_AERSctSul029'] = Spc_Bin_28
df2['SpeciesRst_AERSctSul030'] = Spc_Bin_29

df2['SpeciesRst_AERSctSul031'] = Spc_Bin_30
df2['SpeciesRst_AERSctSul032'] = Spc_Bin_31
df2['SpeciesRst_AERSctSul033'] = Spc_Bin_32
df2['SpeciesRst_AERSctSul034'] = Spc_Bin_33
df2['SpeciesRst_AERSctSul035'] = Spc_Bin_34
df2['SpeciesRst_AERSctSul036'] = Spc_Bin_35
df2['SpeciesRst_AERSctSul037'] = Spc_Bin_36
df2['SpeciesRst_AERSctSul038'] = Spc_Bin_37
df2['SpeciesRst_AERSctSul039'] = Spc_Bin_38
df2['SpeciesRst_AERSctSul040'] = Spc_Bin_39

df2.to_netcdf("GEOSChem.Restart.H2SO4G.40Bins.20150101_0000z.nc4")



