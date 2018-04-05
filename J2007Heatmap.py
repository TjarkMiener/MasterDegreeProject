# -*- coding: utf-8 -*-
#!/usr/bin/python
#
#   J2007Heatmap.py
#
#   Created by Tjark Miener and Bruce Allen on 27.08.17.
#   Copyright (c) 2017 Tjark Miener and Bruce Allen. All rights reserved.
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
    file = open("J2007_%sMHz_Chi.txt" % MHZ, "r")   #"J2007_%sMHz_Chi_finer.txt" for the finer grid
    data = np.loadtxt(file, dtype=str, delimiter='\n')
    #Initialization of the arrays.
    arr_len=181     #201 for the finer grid
    chi_squared=[[-1. for _ in range(arr_len)] for _ in range(arr_len)]
    phi0_minArray=[[-1. for _ in range(arr_len)] for _ in range(arr_len)]
    psi0_minArray=[[-1. for _ in range(arr_len)] for _ in range(arr_len)]
    alpha_array=[[-1. for _ in range(arr_len)] for _ in range(arr_len)]
    zeta_array=[[-1. for _ in range(arr_len)] for _ in range(arr_len)]
    
    #Reading the data and writing it into the arrays.
    counterA=0
    counterB=0
    for line in xrange(len(data)):
        d=data[line].split( )
        alpha_array[counterA][counterB]=float(d[0])
        zeta_array[counterA][counterB]=float(d[1])
        phi0_minArray[counterA][counterB]=int(d[2])
        psi0_minArray[counterA][counterB]=float(d[3])
        chi_squared[counterA][counterB]=float(d[4])
        
        #Updating the counters.
        counterB+=1
        if (counterB==arr_len):
            counterA+=1
            counterB=0

    file.close()
    #Creating alpha and xi arrays for the plt.pcolormesh().
    if MHZ=='1500':
        alpha_start=0.0     #100.0 for the finer grid
        alpha_end=180.0     #120.0 for the finer grid
        zeta_start=0.0      #95.0 for the finer grid
        zeta_end=180.0      #115.0 for the finer grid
    else:
        alpha_start=0.0     #105.0 for the finer grid
        alpha_end=180.0     #125.0 for the finer grid
        zeta_start=0.0      #100.0 for the finer grid
        zeta_end=180.0      #120.0 for the finer grid
    #Length of iteration
    it_len=1.0        #0.1 for the finer grid
    alpha = []
    for i in frange(alpha_start,alpha_end,it_len):
        alpha.append(i)
    alpha = np.asarray(alpha)
    zeta = []
    for x in frange(zeta_start,zeta_end,it_len):
        zeta.append(x)
    zeta = np.asarray(zeta)
    #To get the x-axis and y-axis of the heatmap right, the 2-dim. array has to be transpose.
    chi_squared=np.asarray(chi_squared)
    chi_squared=chi_squared.T

    #Searching for the minimum and the corresponding values (alpha,zeta,phi_0,psi_0).
    chi_squared_min=chi_squared.min()
    chi_squaredindex=np.where(chi_squared==chi_squared_min)
    alpha_minPos=int(chi_squaredindex[1])
    zeta_minPos=int(chi_squaredindex[0])
    alpha_min=float(alpha_array[alpha_minPos][zeta_minPos])
    zeta_min=float(zeta_array[alpha_minPos][zeta_minPos])
    phi0_min=int(phi0_minArray[alpha_minPos][zeta_minPos])
    psi0_min=float(psi0_minArray[alpha_minPos][zeta_minPos])

    #Printing out the values of the minimum.
    print "Value_min: " + str(chi_squared_min)
    print "Alpha_min: " + str(alpha_min)
    print "Zeta_min: " + str(zeta_min)
    print "Phi0_min: " + str(phi0_min)
    print "Psi0_min: " + str(psi0_min)

    #Marking the minimum with a cross in the heatmap and writing the values besides.
    stringChi=r'+ $\, \chi^{2}_{min}=%.3f$' % chi_squared_min
    stringAlpha=r'$\quad \ \alpha_{min}=%d^{\circ}$' % alpha_min     #%.1f^{\circ} for the finer grid
    stringZeta=r'$\quad \ \zeta_{min}=%d^{\circ}$' % zeta_min        #%.1f^{\circ} for the finer grid
    stringPhi=r'$\quad \ \phi_{min}=%d^{\circ}$' % phi0_min
    stringPsi=r'$\quad \ \psi_{min}=%d^{\circ}$' % psi0_min          #%.1f^{\circ} for the finer grid
    plt.text(alpha_min,zeta_min,stringChi,fontsize=15)
    plt.text(alpha_min,zeta_min-10,stringAlpha,fontsize=15)    #zeta_min-1 for the finer grid
    plt.text(alpha_min,zeta_min-20,stringZeta,fontsize=15)     #zeta_min-2 for the finer grid
    plt.text(alpha_min,zeta_min-30,stringPhi,fontsize=15)      #zeta_min-3 for the finer grid
    plt.text(alpha_min,zeta_min-40,stringPsi,fontsize=15)      #zeta_min-4 for the finer grid

    #Creating the heatmap with a logarithmic colorbar.
    alpha, zeta = np.meshgrid(alpha, zeta)
    plt.pcolormesh(alpha, zeta, chi_squared,cmap='RdYlGn',norm=mpl.colors.LogNorm(),vmin=10e-1,vmax=10e2)     #vmax=10e1 for the finer grid
    cbar=plt.colorbar()
    cbar.ax.set_ylabel('$ \quad\chi^2 $', rotation=0, fontsize=24)
    plt.title("%sMHz" % MHZ, fontsize=24)
    plt.xticks([0,30,60,90,120,150,180])    #comment this line out for the finer grid
    plt.yticks([0,30,60,90,120,150,180])    #comment this line out for the finer grid
    plt.xlabel(r'$ \alpha $ [deg]', fontsize=18)
    plt.ylabel(r'$ \zeta $ [deg]', fontsize=18)

    #Saving the heatmap.
    plt.savefig("J2007_%sMHz_ChiPlot.png" % MHZ, dpi = 600)

#Execute the main function.
if __name__ == "__main__":
    main()
