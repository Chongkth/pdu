import pandas as pd
import sys
import numpy as np
import scipy.constants as constants

def RFc(df):
    MeV = constants.eV*1e6
    hbar = constants.physical_constants["reduced Planck constant times c in MeV fm"][0]
    ialpha = constants.physical_constants["inverse fine-structure constant"][0]
    RF = [] #|RF(R)|^2 = 1/fm
    RF2 = [] #10**(RF(R)/R) 
    plist = []
    Xlist = []
    PDU = []
    for n in df.index: 
        Qc = df.loc[n].Qa_calc*1000 #MeV to keV, may differ in database between calculated Qa and experimental Qa
        l = df.loc[n].l
        At_1_2 = df.loc[n].At_1_2
        #print(df.loc[[n], ['logRF(R)']].values[0][0])
#        if df.loc[[n], ['logRF(R)']].values != '':
#            RF.append(df.loc[[n], ['logRF(R)']].values[0][0])
#            plist.append(df.loc[n].rho)
#            Xlist.append(df.loc[n].chi)
#            continue
        if Qc == '' or l == '' or At_1_2 == '': # or '#'in df.loc[n].A_br or df.loc[n].A_br.strip() == '': #or '#' in Qc 
            RF.append('')
            RF2.append('')
            plist.append('')
            Xlist.append('')
            PDU.append('')
            continue
        elif float(Qc) <= 0:
            RF.append('')
            RF2.append('')
            plist.append('')
            Xlist.append('')
            PDU.append('')
            continue
        try:
            logAt_1_2 = np.log10(At_1_2)
        except TypeError as e:
            print(e)
            print(At_1_2)
            sys.exit()
        Qc = np.longdouble(Qc)/1000
        #d = get_daughters(df, df.loc[n].Z, df.loc[n].N)
        #if len(d) == 0:
        #    RF.append('')
        #    continue
        #d = d[0]
        Zd = df.loc[n].Z.astype(int) - 2
        Zc = 2
        Ad = df.loc[n].A.astype(int) - 4
        Ac = 4
        R = 1.2*(Ad**(1/3) + Ac**(1/3))    
        uu = Ad*Ac*938.90595/df.loc[n].A
        v = np.sqrt(2*Qc/(uu)) #*constants.c
        
        p = np.sqrt(2*Qc*uu)*R/hbar 
        X = 2*Zd*Zc*hbar/ialpha*np.sqrt(uu/(2*Qc))/hbar

        Beta = np.arccos(np.sqrt(np.array(Qc)*R*1e-15*MeV*4*np.pi*constants.epsilon_0/(constants.e**2*Zd*Zc)))
        BB = np.arccos(np.sqrt(p/X))
        #l = 0
        H = np.sqrt(1/np.tan(Beta))*np.exp(X*(Beta - np.sin(Beta)*np.cos(Beta)) - l*(l+1)/X*np.tan(Beta))
    #    H1 = np.real(hankel1(X, p))
        logT1_2exp = np.log10(df.loc[n].At_1_2)
        logRF = 1/2*(-logT1_2exp +np.log10(np.log(2)/v*np.abs(H)**2) + np.log10(6.582119*10.0**(-22)) - np.log10(197.327)) 
    #    logRFH = 1/2*(-logT1_2exp +np.log10(np.log(2)/v*np.abs(H1)**2))
        #if l !=0:
        #    print(l)
        #    print(l*(l+1)/X*np.tan(Beta))
        #    print(logRF, '\n')
        RF.append(logRF)
        RF2.append(10**logRF/R)
        PDU.append(8**(0.5)*np.pi**(-7/4)*(0.574)**(-3/4)/R**3)
        plist.append(p)
        Xlist.append(X)
    ratio = []
    for i in range(len(RF2)):
        if RF2[i] == '':
            ratio.append('')
        else:
            ratio.append(RF2[i]/PDU[i])
    
    df['logRF(R)'] = RF
    df['F(R)'] = RF2
    df['rho'] = plist
    df['chi'] = Xlist
    df['PDU'] = PDU
    df['F(R)/PDU'] = ratio
    return df

df = pd.read_excel('NuclearDB.xlsx').fillna('')
df = df.set_index('EL_nu')

df = RFc(df) 
df = df[(df.Z == 84) | (df.Z == 86) | (df.Z == 88)]
df = df[(df.Z%2 == 0) & (df.N%2 == 0)]
df.to_excel('RFRdatabase.xlsx')