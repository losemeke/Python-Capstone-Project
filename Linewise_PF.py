import time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from data_intake import NB, NS, NG, NLB, NT, NSH, TOLER,ITMAX, BIND, PG, QG, \
    PD, QD, DEL, VSH, FB, TB, ZL, BL, TAP, SNO, SUS

    
NSI = int(BIND[NS-1,0])
LNSI = np.concatenate((np.arange(1,NSI) , np.arange(NSI+1,NB+1)))
LNVI = np.array(np.arange(NG+1,NB+1))
LRI = np.concatenate(((np.arange(1,4*NT+1)), (4*NT+LNSI), (4*NT+NB+LNVI)))
# LRI = LRI[np.newaxis,:]
LCI = np.concatenate((LNVI, (LNSI + NB), (np.arange(1,4*NT+1)) + 2*NB))
# LCI = LCI[np.newaxis,:]


FBI  = BIND[FB-1,0].astype(int)
TBI  = BIND[TB-1,0].astype(int)
PAI  = 2*NB+0*NT+np.arange(1,NT);
QAI  = 2*NB+1*NT+np.arange(1,NT);
PBI  = 2*NB+2*NT+np.arange(1,NT);
QBI  = 2*NB+3*NT+np.arange(1,NT);

MF   = sparse.lil_matrix(((np.ones((NT,1)) * np.arange(1,NB+1)) 
                          == (FBI @ np.ones((1,NB))))*1).T
MT   = sparse.lil_matrix(((np.ones((NT,1)) * np.arange(1,NB+1)) 
                          == (TBI @ np.ones((1,NB))))*1).T
MS   = sparse.lil_matrix(((np.ones((NSH,1)) * np.arange(1,NB+1)) 
                          == (BIND[SNO-1,0] @ np.ones((1,NB))))*1).T

YS = MF@BL + MT@BL + MS@SUS*1j + MF@((1/ZL)*(1-TAP)/TAP**2) + \
    MT @((1/ZL)*(TAP-1)/(TAP))
ZL = ZL * TAP
RL = ZL.real
XL = ZL.imag 
ZLM2 = abs(ZL)**2

U   = np.concatenate((VSH[0:NG], np.ones((NLB,1))))**2  
# DEL = np.zeros((NB,1))
Pa  = np.zeros((NT,1))  
Qa  = np.zeros((NT,1)) 
Pb  = np.zeros((NT,1))  
Qb  = np.zeros((NT,1)) 

LJAC = sparse.lil_matrix((2*NB+3*NT+NT,2*NB+3*NT+NT))
LMIS = sparse.lil_matrix((4*NT+NB-1+NLB,1))




def LINEJAC():
    ## FUA
    LJAC[0*NT+0*NB+0:0*NT+0*NB+NT,0*NB+0*NT+0:0*NB+0*NT+NB] \
    = sparse.lil_matrix.multiply(
        MF.T, ((2*U[FBI-1,0] + 2*(Pa*RL + Qa*XL - U[TBI-1,0]/2)) 
               @ np.ones((1,NB))))
    
    LJAC[0*NT+0*NB+0:0*NT+0*NB+NT,0*NB+0*NT+0:0*NB+0*NT+NB] \
    = LJAC[0*NT+0*NB+0:0*NT+0*NB+NT,0*NB+0*NT+0:0*NB+0*NT+NB] + \
    sparse.lil_matrix.multiply(MT.T,(-U[FBI-1,0] @ np.ones((1,NB))))
    
    LJAC[0*NT+0*NB+0:0*NT+0*NB+NT,2*NB+0*NT+0:2*NB+0*NT+NT] \
    = sparse.spdiags((2*U[FBI-1,0]*RL+2*Pa*ZLM2).T,0,NT,NT)
    
    LJAC[0*NT+0*NB+0:0*NT+0*NB+NT,2*NB+1*NT+0:2*NB+1*NT+NT] \
    = sparse.spdiags((2*U[FBI-1,0]*XL+2*Qa*ZLM2).T,0,NT,NT)
    
    
    ## FUB
    LJAC[1*NT+0*NB+0:1*NT+0*NB+NT,0*NB+0*NT+0:0*NB+0*NT+NB] \
    = sparse.lil_matrix.multiply(
        MT.T, ((2*U[TBI-1,0] + 2*(Pb*RL + Qb*XL - U[FBI-1,0]/2))
               @ np.ones((1,NB))))
    
    LJAC[1*NT+0*NB+0:1*NT+0*NB+NT,0*NB+0*NT+0:0*NB+0*NT+NB] \
    = LJAC[1*NT+0*NB+0:1*NT+0*NB+NT,0*NB+0*NT+0:0*NB+0*NT+NB] + \
    sparse.lil_matrix.multiply(MF.T,(-U[TBI-1,0] @ np.ones((1,NB))))
    
    LJAC[1*NT+0*NB+0:1*NT+0*NB+NT,2*NB+2*NT+0:2*NB+2*NT+NT] \
    = sparse.spdiags((2*U[TBI-1,0]*RL+2*Pb*ZLM2).T,0,NT,NT)
    
    LJAC[1*NT+0*NB+0:1*NT+0*NB+NT,2*NB+3*NT+0:2*NB+3*NT+NT] \
    = sparse.spdiags((2*U[TBI-1,0]*XL+2*Qb*ZLM2).T,0,NT,NT)
    
    
    ## FDA
    LJAC[2*NT+0*NB+0:2*NT+0*NB+NT,0*NB+0*NT+0:0*NB+0*NT+NB] \
    = sparse.lil_matrix.multiply(
        MF.T, np.tan(DEL[TBI-1,0]-DEL[FBI-1,0])*np.ones((1,NB)))
    
    LJAC[2*NT+0*NB+0:2*NT+0*NB+NT,1*NB+0*NT+0:1*NB+0*NT+NB] \
    = sparse.lil_matrix.multiply(
        -MF.T, (((U[FBI-1,0]+Pa*RL+Qa*XL)/np.cos(DEL[TBI-1,0]-DEL[FBI-1,0])**2)
                *np.ones((1,NB))))
    
    LJAC[2*NT+0*NB+0:2*NT+0*NB+NT,1*NB+0*NT+0:1*NB+0*NT+NB] \
    = LJAC[2*NT+0*NB+0:2*NT+0*NB+NT,1*NB+0*NT+0:1*NB+0*NT+NB] + \
        sparse.lil_matrix.multiply(
            MT.T, (((U[FBI-1,0]+Pa*RL+Qa*XL)\
                    /np.cos(DEL[TBI-1,0]-DEL[FBI-1,0])**2)*np.ones((1,NB))))
    
    LJAC[2*NT+0*NB+0:2*NT+0*NB+NT,2*NB+0*NT+0:2*NB+0*NT+NT] \
    = sparse.spdiags((RL*np.tan(DEL[TBI-1,0]-DEL[FBI-1,0])-XL).T,0,NT,NT)
    
    LJAC[2*NT+0*NB+0:2*NT+0*NB+NT,2*NB+1*NT+0:2*NB+1*NT+NT] \
    = sparse.spdiags((XL*np.tan(DEL[TBI-1,0]-DEL[FBI-1,0])+RL).T,0,NT,NT)
    
    
    ## FDB
    LJAC[3*NT+0*NB+0:3*NT+0*NB+NT,0*NB+0*NT+0:0*NB+0*NT+NB] \
    = sparse.lil_matrix.multiply(
        MT.T, np.tan(DEL[FBI-1,0]-DEL[TBI-1,0])*np.ones((1,NB)))
    
    LJAC[3*NT+0*NB+0:3*NT+0*NB+NT,1*NB+0*NT+0:1*NB+0*NT+NB] \
    = sparse.lil_matrix.multiply(
        MF.T, (((U[TBI-1,0]+Pb*RL+Qb*XL)/np.cos(DEL[FBI-1,0]-DEL[TBI-1,0])**2)
                *np.ones((1,NB))))
    
    LJAC[3*NT+0*NB+0:3*NT+0*NB+NT,1*NB+0*NT+0:1*NB+0*NT+NB] \
    = LJAC[3*NT+0*NB+0:3*NT+0*NB+NT,1*NB+0*NT+0:1*NB+0*NT+NB] - \
        sparse.lil_matrix.multiply(
            MT.T, (((U[TBI-1,0]+Pb*RL+Qb*XL) \
                    /np.cos(DEL[FBI-1,0]-DEL[TBI-1,0])**2)*np.ones((1,NB))))
    
    LJAC[3*NT+0*NB+0:3*NT+0*NB+NT,2*NB+2*NT+0:2*NB+2*NT+NT] \
    = sparse.spdiags((RL*np.tan(DEL[FBI-1,0]-DEL[TBI-1,0])-XL).T,0,NT,NT)
    
    LJAC[3*NT+0*NB+0:3*NT+0*NB+NT,2*NB+3*NT+0:2*NB+3*NT+NT] \
    = sparse.spdiags((XL*np.tan(DEL[FBI-1,0]-DEL[TBI-1,0])+RL).T,0,NT,NT)
    
    
    ## FP
    LJAC[4*NT+0*NB+0:4*NT+0*NB+NB,0*NB+0*NT+0:0*NB+0*NT+NB] \
    = sparse.spdiags((-YS.real).T,0,NB,NB)
    
    LJAC[4*NT+0*NB+0:4*NT+0*NB+NB,2*NB+0*NT+0:2*NB+0*NT+NT] \
    = MF
    
    LJAC[4*NT+0*NB+0:4*NT+0*NB+NB,2*NB+2*NT+0:2*NB+2*NT+NT] \
    = MT
    
    
    ## FQ
    LJAC[4*NT+1*NB+0:4*NT+1*NB+NB,0*NB+0*NT+0:0*NB+0*NT+NB] \
    = sparse.spdiags((YS.imag).T,0,NB,NB)
    
    LJAC[4*NT+1*NB+0:4*NT+1*NB+NB,2*NB+1*NT+0:2*NB+1*NT+NT] \
    = MF
    
    LJAC[4*NT+1*NB+0:4*NT+1*NB+NB,2*NB+3*NT+0:2*NB+3*NT+NT] \
    = MT
    
LINEJAC()

def LINEMIS():

    LMIS[np.arange(0,(1*NT))] = U[FBI-1,0]**2 + \
    2*U[FBI-1,0]*(Pa*RL+Qa*XL-U[TBI-1,0]/2) + (Pa**2+Qa**2)*ZLM2
    
    LMIS[np.arange(NT,(2*NT))] = U[TBI-1,0]**2 + \
    2*U[TBI-1,0]*(Pb*RL+Qb*XL-U[FBI-1,0]/2) + (Pb**2+Qb**2)*ZLM2

    LMIS[np.arange(2*NT,(3*NT))] =  \
    (U[FBI-1,0]+Pa*RL+Qa*XL)*np.tan(DEL[TBI-1,0]-DEL[FBI-1,0]) - (Pa*XL-Qa*RL)
    
    LMIS[np.arange(3*NT,(4*NT))] =  \
    (U[TBI-1,0]+Pb*RL+Qb*XL)*np.tan(DEL[FBI-1,0]-DEL[TBI-1,0]) - (Pb*XL-Qb*RL)
 
    LMIS[np.arange(4*NT,4*NT+1*(NB-1))] = \
    PG[LNSI-1] - PD[LNSI-1] + MF[LNSI-1,:]*Pa + MT[LNSI-1,:]*Pb - \
    (U[LNSI-1]*YS[LNSI-1]).real
    
    LMIS[np.arange(4*NT+(NB-1),4*NT+(NB-1)+(NB-NG))] = \
    QG[LNVI-1] - QD[LNVI-1] + MF[LNVI-1,:]*Qa + MT[LNVI-1,:]*Qb + \
    (U[LNVI-1]*YS[LNVI-1]).imag

LINEMIS()

IT= 0
ERROR = np.zeros((1,1))
DXX = (spsolve(-LJAC[np.ix_(LRI-1,LCI-1)],LMIS) \
          [np.newaxis,:]).T
ERROR[IT] = max(abs(LMIS)).toarray()[0,0]
# ERROR = np.append(ERROR, max(abs(LMIS)).toarray()[0,0])
t1 = time.time()


while (ERROR[IT]>TOLER and IT <= ITMAX):
    LINEJAC()
    
    DX = (spsolve(-LJAC[np.ix_(LRI-1,LCI-1)],LMIS) \
          [np.newaxis,:]).T
    #DX = (sparse.inv(-LJAC[np.ix_(LRI-1,LCI-1)])) @ LMIS
    
    if (IT == 0):
        DXX = DX
    else:
        DXX = np.append(DXX, DX, axis=1)
        
    U[LNVI-1]   = U[LNVI-1]   + DX[0:NB-NG]
    DEL[LNSI-1] = DEL[LNSI-1] + DX[(NB-NG)+0:(NB-NG)+NB-1]
    Pa = Pa + DX[(NB-NG+NB-1)+0*NT:(NB-NG+NB-1)+1*NT]
    Qa = Qa + DX[(NB-NG+NB-1)+1*NT:(NB-NG+NB-1)+2*NT]
    Pb = Pb + DX[(NB-NG+NB-1)+2*NT:(NB-NG+NB-1)+3*NT]
    Qb = Qb + DX[(NB-NG+NB-1)+3*NT:(NB-NG+NB-1)+4*NT]
    
    LINEMIS()
    IT = IT + 1;
    ERROR = np.append(ERROR, max(abs(LMIS)).toarray()[0,0])
    
t2 = time.time()
print("\nNumber of Iterations: ",IT+1)
print("Time taken: ",(t2-t1)*1000," ms")

V   = np.sqrt(U);

    
    
