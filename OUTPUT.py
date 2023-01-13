import numpy as np
from math import cos,sin, pi
from data_intake import NB, NS, NG, NT, BIND, PG, QG, PD, QD, FB, TB, \
    BL, TAP, YL, BSN, BNAM, PBASE, QGMAX, QGMIN, NCKT, RAT
from Linewise_PF import DEL, V

    
# exec(open(r'data_intake.py').read()) # Change the directory to where your data file is


PI = 3.1415927
PG_TOT = 0
QG_TOT = 0
PD_TOT = 0
QD_TOT = 0
BSNMAX = 0

SLOSS = complex(0, 0)

SLFLOW = (0 + 0j) * np.zeros((NT, 2), dtype=int)
SBUS = (0 + 0j) * np.zeros((NB, 1), dtype=int)

#BIND = np.zeros((int(BSNMAX), 1))

for k in range(0, NT):
    i1 = int(BIND[FB[k]-1])-1
    j1 = int(BIND[TB[k]-1])-1
    vi = V[i1] * complex(cos(DEL[i1]), sin(DEL[i1]))
    vj = V[j1] * complex(cos(DEL[j1]), sin(DEL[j1]))
    SLFLOW[k, 0] = vi * np.conj((vi - vj) * YL[k] / TAP[k])
    SLFLOW[k, 0] = SLFLOW[k, 0] + vi * np.conj(vi * BL[k])
    SLFLOW[k, 0] = SLFLOW[k, 0] + vi * np.conj(vi * YL[k] * (1 / (TAP[k] * TAP[k]) - 1 / TAP[k]))

    SLFLOW[k, 1] = vj * np.conj((vj - vi) * YL[k] / TAP[k])
    SLFLOW[k, 1] = SLFLOW[k, 1] + vj * np.conj(vj * BL[k])
    SLFLOW[k, 1] = SLFLOW[k, 1] + vj * np.conj(vj * YL[k] * (1 - 1 / TAP[k]))

    SBUS[i1] = SBUS[i1] + SLFLOW[k, 0]
    SBUS[j1] = SBUS[j1] + SLFLOW[k, 1]

    SLOSS = SLOSS + SLFLOW[k, 0] + SLFLOW[k, 1]

i1 =  int(BIND[NS-1])-1
PG[i1] = np.real(SBUS[i1])

QG[0:NG] = np.imag(SBUS[0:NG])

print("OUTPUT OF THE LOAD FLOW FDLF METHOD PROGRAM ")
print("                         -------------------------------------------")
print("                                       GENERATOR STATUS")
print("                                       ----------------")
print("------------------------------------------------------------------------------------------")
print("SNO B.NO BUS NAME ---PD--- ---QD--- ---PG--- ---QG--- --QGMX-- --QGMN-- --V SH-- ---DEL---")
print("SNO B.NO BUS NAME   (MW)    (MVAR)    (MW)    (MVAR)   (MVAR)   (MVAR)    (PU)   (DEGREES)")
print("------------------------------------------------------------------------------------------")

for k in range(0, NG):
    print(" ", k, " ", BSN[k], " ", BNAM[k], "%8.2f" % (PD[k] * PBASE), \
          "%8.2f" % (QD[k] * PBASE), "%8.2f" % (PG[k] * PBASE), "%8.2f" % \
              (QG[k] * PBASE), "%8.2f" % (QGMAX[k] * PBASE), "%8.2f" % \
                  (QGMIN[k] * PBASE), "%8.4f" % (V[k]),"%8.4f" % \
                      (DEL[k] * (180.0 / pi)))
    PG_TOT = PG_TOT + PG[k]
    QG_TOT = QG_TOT + QG[k]
    PD_TOT = PD_TOT + PD[k]
    QD_TOT = QD_TOT + QD[k]

print("------------------------------------------------------------------------------------------")

print("                LOAD BUS STATUS")
print("                ---------------")
print("-----------------------------------------------------")
print("SNO B.NO BUS NAME ---PD--- ---QD--- --V SH-- --DEL---")
print("                    (MW)    (MVAR)    (PU)   (DEGREE)")

for k in range(NG , NB ):
    print(" ", k, " ", BSN[k], " ",
          BNAM[k], "%8.2f" % (PD[k] * PBASE), "%8.2f" % (QD[k] * PBASE), \
               "%8.4f" % (V[k]), "%8.4f" % (DEL[k] * (180.0 / pi)))
    PD_TOT = PD_TOT + PD[k]
    QD_TOT = QD_TOT + QD[k]

print("-----------------------------------------------------")

print("TRANSMISSION LINE FLOW  ---  FLOWS ARE PER CIRCUIT ")
print("-------------------------------------------------- ")
print("----------------------------------------------------------------------------------")
print("SNo FBNO --NAME-- TBNO --NAME-- NCKT      FROM FLOW         TO  FLOW          MVA ")
print("                                     -------- -------- -------- --------     FLOW ")
print("                                        REAL  REACTIVE   REAL  REACTIVE       PER ")
print("                                        (MW)   (MVAR)    (MW)   (MVAR)        CKT ")
print("----------------------------------------------------------------------------------")

for k in range(0, NT):
    k1 = int(BIND[FB[k]-1])-1
    k2 = int(BIND[TB[k]-1])-1
    print(" ", k, " ", BSN[k1], " ", BNAM[k1], " ", BSN[k2], " ", BNAM[k2], " ", " ", NCKT[k],
          "%8.2f" % (np.real(SLFLOW[k, 0]) * PBASE), \
              "%8.2f" % (np.imag(SLFLOW[k, 0]) * PBASE), \
                  "%8.2f" % (np.real(SLFLOW[k, 1]) * PBASE), \
                      "%8.2f" % (np.imag(SLFLOW[k, 1]) * PBASE),\
                          "%8.2f" % (abs(SLFLOW[k, 0]) * PBASE), end='  ')
    if abs(SLFLOW[k, 0]) < (0.1 * RAT[k]):
        print("UNDER - %4.2f "% (100 * abs(SLFLOW[k, 0]) / (RAT[k])))

    if abs(SLFLOW[k, 1]) > 1.1 * RAT[k]:
        print("OVER  - %4.2f "% (100 * abs(SLFLOW[k, 0]) / (RAT[k])))

    if abs(SLFLOW[k, 0]) > 0.1 * RAT[k] and abs(SLFLOW[k, 0]) < 1.1 * RAT[k]:
        print(" ")

print("----------------------------------------------------------------------------------")

print("SUMMARY OF THE RESULTS")
print("----------------------")

print("----------------------------------------------")
print("TOTAL REAL GENERATION     = %15.5f"% (PG_TOT * PBASE))
print("TOTAL REAL LOAD           = %15.5f"% (PD_TOT * PBASE))
print("TOTAL REAL LOSS           = %15.5f"% (np.real(SLOSS) * PBASE))
print("----------------------------------------------")
print("TOTAL REACTIVE GENERATION = %15.5f"% (QG_TOT * PBASE))
print("TOTAL REACTIVE LOAD       = %15.5f"% (QD_TOT * PBASE))
print("TOTAL REACTIVE LOSS       = %15.5f"% (np.imag(SLOSS) * PBASE))
# print("TOTAL SHUNT GENERATION    = %15.5f", sh_gen*PBASE);
print("----------------------------------------------")
print("TRANSMISSION LOSS         = %15.5f" %(np.real(SLOSS) * PBASE))
print("----------------------------------------------")
