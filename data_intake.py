#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 22:10:11 2021

@author: vineet
"""
import numpy;
# import math;

global NB, NBB, NS, NG, NLB, NTR, NTRL, NT, NSHC, NSVS, NSHR, NSH, NREG
global VSLACK, TOLER, PBASE, VLMAX, VLMIN, ITMAX
global BIND, BSN, BNAM, PG, QG, PD, QD, V, DEL, QGMAX, QGMIN, VSH                
global FB, TB, NCKT, YL, ZL, BL, TAP, RAT, KV, LEN, TAPMAX, TAPMIN, TAPSTP        
global SNO, SUS, SUSMAX, SUSMIN, SUSSTP
global YB
global LAM
global PGMAX, PGMIN, VGMAX, VGMIN, ap, bp, cp 
global ofp, ifp                          

# temp1 = input("Input  File Name: ");
# temp2 = input("Output File Name: ");

temp1 = '6bus.i'
temp2 = 'test1'

ifp = open(temp1,'r');
ofp = open(temp2,'w');
    
with ifp as fi:
    lines = fi.readlines();
    line0 = lines[0];
    line1 = lines[1];
    line2 = lines[2];
    line3 = lines[3];
    line4 = lines[4];
    line5 = lines[5];
    
ttt1 = format(line4.strip());    
ttt2 = ttt1.split();    
ttt = list(map(int, ttt2));
    
NB = ttt[0];
NBB = ttt[1];
NS = ttt[2];
NG = ttt[3];
temp = ttt[4];
temp = ttt[5];
NLB = ttt[6];
NTR = ttt[7];
NTRL = ttt[8];
NSHC = ttt[9];
NSVS = ttt[10];
NSHR = ttt[11];
NREG = ttt[12];
    
NT = NTR + NTRL;
NSH = NSHC + NSVS + NSHR;

ttt1 = format(line5.strip());
ttt2 = ttt1.split();
ttt = list(map(float, ttt2));

VSLACK = ttt[0];
TOLER = ttt[1];
PBASE = ttt[2];
VLMAX = ttt[3];
VLMIN = ttt[4];
ITMAX = int(ttt[5]);

BSN   = numpy.zeros((NB,1),dtype=int);  
BNAM  = list(map(str,numpy.arange(NB)));
PD    = numpy.zeros((NB,1),dtype=float);  
QD    = numpy.zeros((NB,1),dtype=float);
PG    = numpy.zeros((NB,1),dtype=float);  
QG    = numpy.zeros((NB,1),dtype=float);
QGMAX = numpy.zeros((NB,1),dtype=float);  
QGMIN = numpy.zeros((NB,1),dtype=float);
VSH   = numpy.zeros((NB,1),dtype=float);

V     = numpy.ones((NB,1));   
DEL   = numpy.zeros((NB,1));

BSNMAX = 0;
   
with ofp as fo:    
    print("INPUT  FILE NAME: ", temp1, file=fo);
    print("OUTPUT FILE NAME: ", temp2, file=fo);
    print("SYSTEM: {}".format(line0.strip()), file=fo);
    print("YEAR: {}".format(line1.strip()), file=fo);    
    print("CASE: {}".format(line2.strip()), file=fo); 
    print("NUMBER: {}".format(line3.strip()), file=fo); 
    print("NUMBER OF BUSES : ", NB, file=fo);
    print("SLACK BUS NUMBER : ", NS, file=fo);
    print("NUMBER OF GENERATORS : ", NG, file=fo);
    print("NUMBER OF LOAD BUSES : ", NLB, file=fo);
    print("NUMBER OF TRANSFORMERS : ", NTR, file=fo);
    print("NUMBER OF TRANSMISSION LINES : ", NTRL, file=fo);
    print("NUMBER OF SHUNT CAPACITORS : ", NSHC, file=fo);
    print("NUMBER OF SWITCHABLE CAPACITORS : ", NSVS, file=fo);
    print("NUMBER OF SHUNT REACTORS : ", NSHR, file=fo);
    print("SLACK BUS VOLATGE : {:8.4f}".format(VSLACK), file=fo);
    print("TOLERANCE (MW) : {:8.4f}".format(TOLER*PBASE), file=fo);
    print("BASE (MVA)    60 : {:8.4f}".format(PBASE), file=fo);
    print("MINIMUM LOAD BUS VOLTAGE : {:8.4f}".format(VLMIN), file=fo);
    print("MAXIMUM LOAD BUS VOLTAGE : {:8.4f}".format(VLMAX), file=fo);
    print("MAXIMUM NUMBER OF ITERATIONS : ", ITMAX, file=fo);
    print("DETAILED OUTPUT IN PER UNIT  ", file=fo);
    print("GENERATOR BUS DATA", file=fo);
    print("SNO R.NO   B.NO BUS NAME ---PG--- ---PD--- ---QD--- --QGMX--- ---QGMN--- --V SH--", file=fo);    
   
    for k in range(0, NG):
        line = list(map(str, (format(lines[7+k])).split()));
        BSN[k] = int(line[2]);
        BNAM[k] = line[3];
        PD[k] = float(line[4]);
        QD[k] = float(line[5]);
        PG[k] = float(line[6]);      
        QG[k] = 0;  
        QGMAX[k] = float(line[7]);      
        QGMIN[k] = float(line[8]);
        VSH[k] = float(line[9]);      
        V[k] = VSH[k];
    
        if BSNMAX < BSN[k]:
            BSNMAX = BSN[k];
       
        print(" ", k+1,"  ", k+1,"", BSN[k],"  ", BNAM[k],"%8.4f" %PG[k],"%8.4f" %PD[k],"%8.4f" %QD[k],"%8.4f" %QGMAX[k],"%8.4f" %QGMIN[k],"%8.4f" %VSH[k], file=fo);
    
        PD[k] = PD[k]/PBASE;     
        QD[k] = QD[k]/PBASE;
        PG[k] = PG[k]/PBASE;     
        QG[k] = QG[k]/PBASE;
        QGMAX[k] = QGMAX[k]/PBASE;     
        QGMIN[k] = QGMIN[k]/PBASE;

    
    print("LOAD BUS DATA", file=fo);
    print("SNO R.NO   B.NO BUS NAME ---PD--- ---QD---", file=fo);
    
    for k in range(NG, NB):
        line = list(map(str, (format(lines[7+k])).split()));
        BSN[k] = int(line[2]);
        BNAM[k] = line[3];
        PD[k] = float(line[4]);
        QD[k] = float(line[5]);

        if BSNMAX < BSN[k]:
            BSNMAX = BSN[k];

        print(" ", k+1, "  ", k+1, "", BSN[k], "  ", BNAM[k], "%8.4f" %PD[k], "%8.4f" %QD[k], file=fo);
        
        PD[k] = PD[k]/PBASE;
        QD[k] = QD[k]/PBASE;
    
    BIND = numpy.zeros((int(BSNMAX),1));
    for k in range(0,NB):
        BIND[BSN[k]-1] = k+1;
        
    ZL = (0+0j)  * numpy.zeros((NT,1));
    YL = (0+0j)  * numpy.zeros((NT,1));
    BL = (0+0j)  * numpy.zeros((NT,1));
    
    FB = numpy.zeros((NT,1), dtype=int);   
    TB = numpy.zeros((NT,1), dtype=int);
    NCKT = numpy.zeros((NT,1), dtype=int);  
    TAP = numpy.ones((NT,1));
    RAT = numpy.zeros((NT,1));   
    KV = numpy.zeros((NT,1));
    LEN = numpy.zeros((NT,1)); 
    
    TAPMAX = numpy.ones((NT,1));
    TAPMIN = numpy.ones((NT,1));
    TAPSTP = numpy.ones((NT,1));
    
    print("TRANSFORMER DATA FOR TOTAL NOS OF CIRCUITS", file=fo);
    print("SNO F.NO   T.NO NCKT --R PU-- --X PU-- ---AN--- --RAT---", file=fo);
    
    m = 0;
    l = 1;
    
    for k in range(0, NTR):
        line = list(map(float, (format(lines[7+NB+m])).split()));
        m = 2+m;
        
        FB[k] = int(line[1]);
        TB[k] = int(line[2]);
        NCKT[k] = int(line[3]);
        ZL[k] = complex(line[4],line[5]);
        TAP[k] = float(line[6]);
        BL[k] = complex(0,0);
        RAT[k] = float(line[7])/PBASE;
        LEN[k] = 0;
        KV[k] = 0;
        
        BL[k] = BL[k]*NCKT[k];
        ZL[k] = ZL[k]/NCKT[k];
        
        YL[k] = 1/ZL[k];
        
        lineTAP = list(map(float, (format(lines[7+NB+l])).split())); 
        l = 2+l;
        
        TAPMAX[k] = float(lineTAP[0]);
        TAPMIN[k] = float(lineTAP[1]);
        TAPSTP[k] = float(lineTAP[2]);
        
        print("", k+1, FB[k], TB[k], NCKT[k], "%8.4f" %ZL[k].real, "%8.4f" %ZL[k].imag, "%8.4f" %TAP[k], "%8.4f" %RAT[k], file=fo);
        
    print("LINE DATA FOR TOTAL NOS OF CIRCUITS", file=fo);
    print("SNO F.NO   T.NO NCKT --R PU-- --X PU-- --HLC--- --RAT--- ---AN---", file=fo);
    
    m = 0;
    for k in range(NTR, NT):
        line = list(map(float, (format(lines[7+NB+NTR+NTR+m])).split()));
        m = m+1;
        
        FB[k] = line[1];            
        TB[k]  = line[2];
        NCKT[k] = line[3];
        ZL[k] = complex(line[4],line[5]);
        BL[k] = complex(0,line[6]);
        RAT[k] = line[7]/PBASE;      
        LEN[k] = line[8];
        KV[k] = line[9];
      
        BL[k] = BL[k] * NCKT[k];        
        ZL[k] = ZL[k]/NCKT[k];
        
        if LEN[k]>0 and KV[k]>0:
            ZBASE = KV[k] * KV[k]  / PBASE;       
            ZL[k] = ZL[k] * LEN[k] / ZBASE;
            BL[k] = BL[k] * LEN[k] * ZBASE;
            
        YL[k] = 1/ZL[k];
        
        print("", k+1, FB[k], TB[k], NCKT[k], "%8.4f" %ZL[k].real, "%8.4f" %ZL[k].imag, "%8.4f" %BL[k].real, "%8.4f" %RAT[k], "%8.4f" %TAP[k], file=fo);
    
    SNO = numpy.zeros((NSH,1), dtype=int);    
    SUS = numpy.zeros((NSH,1));
    SUSMAX = numpy.zeros((NSH,1));    
    SUSMIN = numpy.zeros((NSH,1));
    SUSSTP = numpy.zeros((NSH,1));  
    
    print("SHUNT CAPACITOR DATA", file=fo);
    print("SNO B.NO  -MVAR-pu-", file=fo);
    
    for k in range(0, NSHC):
        line = list(map(float, (format(lines[7+NB+NTR+NTR+NTRL+k])).split()));
        SNO[k] = line[1];            
        SUS[k] = line[2]/PBASE;  
        
        print("", k+1, SNO[k], "%8.4f" %SUS[k], file=fo);
    
    print("SWITCHABLE CAPACITOR DATA", file=fo);
    print("SNO B.NO --MAX--- --MIN--- --STEP-- -ACTUAL-", file=fo);
    print("         --MVAR-- --MVAR-- --MVAR-- --MVAR--", file=fo);

    for k in range(NSHC, (NSHC+NSVS)):
        line = list(map(float, (format(lines[7+NB+NTR+NTR+NTRL+NSHC+1+k])).split()));
        SNO[k] = line[1];            
        SUSMAX[k] = line[2] / PBASE;
        SUSMIN[k] = line[3] / PBASE;
        SUSSTP[k] = line[4] / PBASE;
        SUS[k] = line[5] / PBASE;   
        
        print("", k+1, SNO[k], "%8.4f" %SUSMAX[k], "%8.4f" %SUSMIN[k], "%8.4f" %SUSSTP[k], "%8.4f" %SUS[k], file=fo);
   
    
    print("SHUNT REACTOR DATA", file=fo);
    print("SNO B.NO  --MVAR--", file=fo);
    
    m = 0;
    for k in range((NSHC+NSVS), NSH):
        line = list(map(float, (format(lines[7+NB+NTR+NTR+NTRL+NSHC+1+NSVS+m])).split())); 
        m = m+1;
        
        SNO[k] = line[1];            
        SUS[k] = -line[2]/PBASE;  
        
        print("", k+1, SNO[k], "%8.4f" %SUS[k], file=fo);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    