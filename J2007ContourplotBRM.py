# -*- coding: utf-8 -*-
#!/usr/bin/python
#
#   J2007ContourplotBRM.py
#
#   Created by Tjark Miener and Bruce Allen on 06.02.18.
#   Copyright (c) 2018 Tjark Miener and Bruce Allen. All rights reserved.
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

#This function is doing the same as the function range(), but with floats.
def frange(start, stop, step):
    x = start
    while x <= stop:
        yield x
        x += step

def main():
    #Change MHZ for working with the 2000MHz data!
    MHZ = '1500'
    file = open("J2007_BRM%sChi.txt" % MHZ, "r")
    data = np.loadtxt(file, dtype=str, delimiter='\n')
    
    #Initialization of the arrays.
    e2_len=29
    psi2_len=31
    chi_squared=[[-1. for _ in range(psi2_len)] for _ in range(e2_len)]
    e2_Array=[[-1. for _ in range(psi2_len)] for _ in range(e2_len)]
    psi2_Array=[[-1. for _ in range(psi2_len)] for _ in range(e2_len)]
    #Reading the data and writing it into the arrays.
    counterA=0
    counterB=0
    for line in xrange(len(data)):
        d=data[line].split( )
        e2_Array[counterA][counterB]=float(d[0])
        psi2_Array[counterA][counterB]=float(d[1])
        chi_squared[counterA][counterB]=float(d[2])
        #Updating the counters.
        counterB+=1
        if (counterB==psi2_len):
            counterA+=1
            counterB=0
    file.close()

    #Creating e2 and psi2 arrays for the plt.contourf().
    e2_start=0.1
    e2_end=1.51
    e2_it=0.05
    if MHZ=='1500':
        psi2_start=40.0
        psi2_end=70.0
    else:
        psi2_start=45.0
        psi2_end=75.0
    psi2_it=1.0

    e2 = []
    for i in frange(e2_start,e2_end,e2_it):
        e2.append(i)
    e2 = np.asarray(e2)

    psi2 = []
    for x in frange(psi2_start,psi2_end,psi2_it):
        psi2.append(x)
    psi2 = np.asarray(psi2)

    #To get the x-axis and y-axis of the heatmap right, the 2-dim. array has to be transpose.
    chi_squared=np.asarray(chi_squared)
    chi_squared=chi_squared.T

    #Searching for the minimum and the corresponding values (e2,psi2).
    chi_squared_min=chi_squared.min()
    chi_squaredindex=np.where(chi_squared==chi_squared_min)
    e2_minPos=int(chi_squaredindex[1])
    psi2_minPos=int(chi_squaredindex[0])
    e2_min=float(e2_Array[e2_minPos][psi2_minPos])
    psi2_min=float(psi2_Array[e2_minPos][psi2_minPos])

    #Printing out the values of the minimum.
    print "Value_min: " + str(chi_squared_min)
    print "E2_min: " + str(e2_min)
    print "Psi2_min: " + str(psi2_min)
    #Marking the minimum with a cross in the contourplot and writing the values besides.
    stringChi=r'+ $\, \chi^{2}_{min}=%.3f$' % chi_squared_min
    stringE2=r'$\quad \ |\vec{E}_{BG}|^{2}=%.2f$' % e2_min
    stringPsi2=r'$\quad \ \psi_{BG}=%d^{\circ}$' % psi2_min
    plt.text(e2_min,psi2_min,stringChi,fontsize=15)
    plt.text(e2_min,psi2_min-2,stringE2,fontsize=15)
    plt.text(e2_min,psi2_min-4,stringPsi2,fontsize=15)

    #Creating the contourplot.
    e2, psi2 = np.meshgrid(e2, psi2)
    clevs = np.arange(0.5,10.55,0.1)
    plt.contourf(e2,psi2,chi_squared,clevs,cmap='RdYlGn')
    cbar=plt.colorbar(ticks=[0.5,2.5,4.5,6.5,8.5,10.5])
    cbar.ax.set_ylabel('$ \quad\chi^2 $', rotation=0, fontsize=24)
    plt.title("%sMHz" % MHZ, fontsize=24)
    plt.xlim(0.1,1.5)
    if MHZ=='1500':
        plt.ylim(40,70)
    else:
        plt.ylim(45,75)
    plt.xlabel(r'$ |\vec{E}_{BG}|^{2} $ [mJy]', fontsize=16)
    plt.ylabel(r'$ \psi_{BG} $ [deg]', fontsize=16)

    #Saving the heatmap.
    plt.savefig("J2007_%sMHz_ChiBRMPlot.png" % MHZ, dpi = 600)

#Execute the main function.
if __name__ == "__main__":
    main()
