# -*- coding: utf-8 -*-
#!/usr/bin/python
#
#   J2007DisplayData.py
#
#   Created by Tjark Miener and Bruce Allen on 10.08.17.
#   Copyright (c) 2017 Tjark Miener and Bruce Allen. All rights reserved.
#

import numpy as np
import matplotlib.pyplot as plt

#Degrees to Radians conversion.
def deg2rad(val):
    return val*np.pi/180
#Radians to Degrees conversion.
def rad2deg(val):
    return val*180/np.pi
#Equation of the theoretical polarization position angle psi. All arguments in the head of the function are in radians, except the variable phi.
def position_angle(alpha,zeta,phi,phi0,psi0):
    x=-(np.sin(zeta)*np.cos(alpha)-np.cos(zeta)*np.sin(alpha)*np.cos(phi-phi0))
    y=np.sin(phi-phi0)*np.sin(alpha)
    return np.arctan2(y,x) + psi0

def main():
    #Change MHZ for working with the 2000MHz profile!
    MHZ = '1500'
    file = open("J2007_%sMHz_profile.txt" % MHZ, "r")
    data = np.loadtxt(file, dtype=str, delimiter='\n')
    #Reading the data and writing it into arrays.
    for line in xrange(0, len(data)):
        if line==0:
            d=data[line].split( )
            intensitaetI=[float(d[3])]
            linearQ=[float(d[4])]
            linearU=[float(d[5])]
            linearL=[np.sqrt(np.power(float(d[4]),2)+np.power(float(d[5]),2))]
            circularV=[float(d[6])]
            psi=[float(d[7])]
            err=[float(d[8])]
        else:
            d=data[line].split( )
            intensitaetI.append(float(d[3]))
            linearQ.append(float(d[4]))
            linearU.append(float(d[5]))
            linearL.append(np.sqrt(np.power(float(d[4]),2)+np.power(float(d[5]),2)))
            circularV.append(float(d[6]))
            psi.append(float(d[7]))
            err.append(float(d[8]))

    file.close()
    
    #In the following only non-zero value of psi and error of psi are used.
    psi=np.asarray(psi)
    err=np.asarray(err)
    psi_index=np.where(psi!=0.0)
    err_index=np.where(err!=0.0)

    #The pulse phase phi is discretized with the help of the index x (number of measurements).
    x=np.arange(0.0,256.0,1.0)
    phi=2*np.pi*x/256
    #Calculating the theoretical psi with the resulting values of the RVM analysis.
    if MHZ=='1500':
        alpha=deg2rad(111.6)
        zeta=deg2rad(104.1)
        phi0=deg2rad(193)
        psi0=deg2rad(-12.9)
        psi_theo=position_angle(alpha,zeta,phi,phi0,psi0)
        print "J2007_%sMHz_profile.txt" % MHZ
    if MHZ=='2000':
        alpha=deg2rad(115.1)
        zeta=deg2rad(109.7)
        phi0=deg2rad(202)
        psi0=deg2rad(-5.5)
        psi_theo=position_angle(alpha,zeta,phi,phi0,psi0)
        print "J2007_%sMHz_profile.txt" % MHZ
    
    #Creating the plots in figure 1 to see the results.
    plt.figure(figsize=(9,7.2))
    #The first subplot (top) is the measured and the theoretical polarization angle psi as a function of pulse phase.
    plt.subplot(211)
    plt.plot(phi/(2*np.pi),rad2deg(psi_theo),'g--',label='theo PA')
    plt.errorbar((x/256)[psi_index], psi[psi_index], err[err_index], fmt='g.', ecolor='g',label='measured PA')
    plt.xlim(0,1)
    plt.ylim(-90+psi0,90+psi0)
    plt.yticks([-80,-40,0,40,80])
    plt.xlabel('Pulse phase [periods]')
    plt.ylabel(r'$ \psi $ [deg]',fontsize=18)
    plt.xlim(0,1)
    #The second subplot (bottom) is the different radio flux-density S as a function of pulse phase.
    plt.subplot(212)
    plt.plot(x/256,intensitaetI,'g-',label='I')
    plt.plot(x/256,linearL,'y--',label='L')
    plt.plot(x/256,circularV,'k-.',label='V')
    plt.text(0.03,6.2,"%sMHz" % MHZ,fontsize=14)
    plt.xlim(0,1)
    plt.ylim(-1,7)
    plt.yticks([0,3,6])
    lg=plt.legend(loc='upper right')
    lg.draw_frame(False)
    plt.xlabel('Pulse phase [periods]')
    plt.ylabel(r'$ S $ [mJy]',fontsize=16)

#   Saving the plots.
    plt.savefig("J2007_%sMHz_profilePlot.png" % MHZ, dpi = 600)

#Execute the main function.
if __name__ == "__main__":
    main()
