//
//  chi_squaredBRM.c
//
//  Created by Tjark Miener and Bruce Allen on 02.02.18.
//  Copyright (c) 2018 Tjark Miener and Bruce Allen. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct Data{
  int index;       //0,...,255
  float intensityI;
  float linearQ;
  float linearU;
  float linearL;
  float circularV;
  float psi;
  float psi_error;
  float e1;
};
struct Data list[256];

//Degrees to Radians conversion.
double deg2rad(float val){
  return val*M_PI/180.0;
}

//Radians to Degrees conversion.
float rad2deg(double val){
  return val*180.0/M_PI;
}

//All arguments in the head of the function are in radians, except the variable called 'index'. The returning position angle is declared in radians.
double positionAngle(float alpha, float zeta, int index, float phi0, float psi0) {
  double phi_rad=index*2*M_PI/256.0;
  //Calculating the position angle with the rotating vector model (RVM). The result is declared in radians.
  double x=-(sin(zeta)*cos(alpha)-cos(zeta)*sin(alpha)*cos(phi_rad-phi0));
  double y=sin(phi_rad-phi0)*sin(alpha);
  return atan2(y,x) + psi0;
}

//main function
int main (int argc, char** argv) {
  //Initialization of helpful variables.
  int i,counter,counterA,counterB,counterE2,counterP2;
  float alpha_val[3],zeta_val[3],psi0_val[3],e2_val[3],psi2_val[3];
  int phi0_val[3];
  
  //Bin size for the contourplot.
  int bins=21;
  int e2_bins=29;
  int psi2_bins=31;
    
  double chi[bins][bins];
  float alpha_Array[bins][bins];
  float zeta_Array[bins][bins];
  int phi0_Array[bins][bins];
  float psi0_Array[bins][bins];
  double chi_min[e2_bins][psi2_bins];
  float alpha_minArray[e2_bins][psi2_bins];
  float zeta_minArray[e2_bins][psi2_bins];
  int phi0_minArray[e2_bins][psi2_bins];
  float psi0_minArray[e2_bins][psi2_bins];
    
  FILE *input_file=fopen(argv[1],"r");
  char str[35];
    
  float alpha_theo,zeta_theo,phi0_theo,psi0_theo;
  //Detecting if the input file is the profile of 1500MHz or 2000MHz.
  if (argv[1][6]=='1') {
    sprintf(str, "J2007_BRM1500Chi.txt");
    alpha_theo=deg2rad(111.6);
    zeta_theo=deg2rad(104.1);
    phi0_theo=deg2rad(193.0);
    psi0_theo=deg2rad(-12.9);
    e2_val[0]=0.1;
    e2_val[1]=1.5;
    e2_val[2]=0.05;
    psi2_val[0]=40.0;
    psi2_val[1]=70.0;
    psi2_val[2]=1.0;
  } else if (argv[1][6]=='2') {
    sprintf(str, "J2007_BRM2000Chi.txt");
    alpha_theo=deg2rad(115.1);
    zeta_theo=deg2rad(109.7);
    phi0_theo=deg2rad(202.0);
    psi0_theo=deg2rad(-5.5);
    e2_val[0]=0.1;
    e2_val[1]=1.5;
    e2_val[2]=0.05;
    psi2_val[0]=45.0;
    psi2_val[1]=75.0;
    psi2_val[2]=1.0;
  } else {
    printf("Error! Inputfile %s is wrong!\n", argv[1]);
    exit(1);
  }
  //This ranges are the same for both observing frequencies.
  alpha_val[0]=100.0;
  alpha_val[1]=120.0;
  alpha_val[2]=1.0;
  zeta_val[0]=100.0;
  zeta_val[1]=120.0;
  zeta_val[2]=1.0;
  phi0_val[0]=100;
  phi0_val[1]=200;
  phi0_val[2]=1;
  psi0_val[0]=-20.0;
  psi0_val[1]=20.0;
  psi0_val[2]=1.0;
    
  //Open the output file.
  FILE *output_file=fopen(str,"w");
  if (!input_file) {
    printf("Error! Inputfile %s not found\n", argv[1]);
    exit(2);
  }
  if (!output_file) {
    printf("Error! Outputfile %s not found\n", str);
    exit(3);
  }
    
  //Reading the input data and writing it into the structure (struct Data).
  counter=0;
  while (1) {
    int foobar1,foobar2;
    if (9==(i=fscanf(input_file,"%d %d %d %f %f %f %f %f %f\n",&foobar1,&foobar2,&list[counter].index,&list[counter].intensityI,&list[counter].linearQ,&list[counter].linearU,&list[counter].circularV,&list[counter].psi,&list[counter].psi_error))) {
        
      list[counter].linearL=sqrt(pow(list[counter].linearQ,2)+pow(list[counter].linearU,2));
      counter++;
    } else break;
  }
    
  if (i!=EOF) {
    fprintf(stderr, "Error! Problem reading %s at line %d . fscanf() returned %d\n", argv[1], counter+1,i);
    exit(4);
  }

  //Calculating the prefactor of the Chi^2 Equation.
  float n=0.0;
  for (int i=0; i<=255;i++) {
    if (list[i].psi!=0 && list[i].psi_error!=0) {
      n+=1.0;
    }
  }
  double prefactor=1.0/(n-6.0);
  //Minimizing Chi^2 over the parameter space (alpha,zeta,phi0,psi0) plus E2 and psi2.
  counterE2=-1;
  for (float e2=e2_val[0]; e2<=e2_val[1]; e2+=e2_val[2]) {
    counterE2+=1;
    counterP2=-1;
    for (float psi2=psi2_val[0]; psi2<=psi2_val[1]; psi2+=psi2_val[2]) {
      printf("Iteration:%.2f %.1f\n",e2,psi2);
      counterP2+=1;
      double e3;
      double y[256],psi_rvm[256];
      for (int i=0; i<=255; i++) {
        e3=list[i].linearL;
        //If there is no measured psi, the theoretical psi from the RVM by Allen and Knispel is taken. This is an acceptable approximation!
        if (list[i].psi==0.0) {
          psi_rvm[i]=positionAngle(alpha_theo, zeta_theo, i, phi0_theo, psi0_theo);
          if (psi_rvm[i]>0.0) {
            y[i]=e3+e2-2*sqrt(e3)*sqrt(e2)*cos(fabs(deg2rad(psi2)-psi_rvm[i]));
          } else {
            y[i]=e3+e2-2*sqrt(e3)*sqrt(e2)*cos(fabs(psi_rvm[i])+deg2rad(psi2));
          }
          list[i].e1=y[i];
        } else if (list[i].psi<0.0) {
          list[i].e1=e3+e2-2*sqrt(e3)*sqrt(e2)*cos(fabs(deg2rad(list[i].psi))+deg2rad(psi2));
        } else {
          list[i].e1=e3+e2-2*sqrt(e3)*sqrt(e2)*cos(fabs(deg2rad(psi2)-deg2rad(list[i].psi)));
        }
      }
      counterA=-1;
      for (float alpha=alpha_val[0]; alpha<=alpha_val[1]; alpha+=alpha_val[2]) {
        double alpha_rad=deg2rad(alpha);
        counterA+=1;
        counterB=-1;
        for (float zeta=zeta_val[0]; zeta<=zeta_val[1]; zeta+=zeta_val[2]) {
          double zeta_rad=deg2rad(zeta);
          counterB+=1;
          for (int phi0=phi0_val[0]; phi0<=phi0_val[1]; phi0+=phi0_val[2]) {
            double phi0_rad=deg2rad(phi0);
            for (float psi0=psi0_val[0]; psi0<=psi0_val[1]; psi0+=psi0_val[2]) {
              double psi0_rad=deg2rad(psi0);
              //Calculating the Chi^2 for one point of the parameter space.
              double psi_theo[256];
              double sum=0.0;
              for (int i=0; i<=255; i++) {
                if (list[i].psi!=0 && list[i].psi_error!=0) {
                  psi_rvm[i]=positionAngle(alpha_rad, zeta_rad, i, phi0_rad, psi0_rad);
                  psi_theo[i]=rad2deg(atan2(list[i].e1*sin(psi_rvm[i])+e2*sin(deg2rad(psi2)),list[i].e1*cos(psi_rvm[i])+e2*cos(deg2rad(psi2))));
                  sum+=(pow((list[i].psi-psi_theo[i])/list[i].psi_error,2));
                }
              }
              //Initializing Chi^2.
              if (phi0==phi0_val[0] && psi0==psi0_val[0]) {
                chi[counterA][counterB]=prefactor*sum;
                alpha_Array[counterA][counterB]=alpha;
                zeta_Array[counterA][counterB]=zeta;
                phi0_Array[counterA][counterB]=phi0;
                psi0_Array[counterA][counterB]=psi0;
              } else {
                //Minimizing Chi^2.
                if (chi[counterA][counterB] > prefactor*sum) {
                  chi[counterA][counterB]=prefactor*sum;
                  alpha_Array[counterA][counterB]=alpha;
                  zeta_Array[counterA][counterB]=zeta;
                  phi0_Array[counterA][counterB]=phi0;
                  psi0_Array[counterA][counterB]=psi0;
                }
              }
            }
          }
        }
      }
      //Detecting the minimal Chi^2 for each e2 and psi2.
      for (int i=0; i<bins; i++) {
        for (int j=0; j<bins; j++) {
          if (i==0 && j==0) {
              chi_min[counterE2][counterP2]=chi[i][j];
              alpha_minArray[counterE2][counterP2]=alpha_Array[i][j];
              zeta_minArray[counterE2][counterP2]=zeta_Array[i][j];
              phi0_minArray[counterE2][counterP2]=phi0_Array[i][j];
              psi0_minArray[counterE2][counterP2]=psi0_Array[i][j];
          } else {
            if (chi_min[counterE2][counterP2] > chi[i][j]) {
              chi_min[counterE2][counterP2]=chi[i][j];
              alpha_minArray[counterE2][counterP2]=alpha_Array[i][j];
              zeta_minArray[counterE2][counterP2]=zeta_Array[i][j];
              phi0_minArray[counterE2][counterP2]=phi0_Array[i][j];
              psi0_minArray[counterE2][counterP2]=psi0_Array[i][j];
            }
          }
        }
      }
    }
  }
  //Printing the minimized Chi_min^2 into the output file. The python code 'J2007ContourplotBRM.py' will read in this file and visualize the Chi_min^2 using a heatmap.
  counterE2=-1;
  for (float e2=e2_val[0]; e2<=e2_val[1]; e2+=e2_val[2]) {
    counterE2+=1;
    counterP2=-1;
    for (float psi2=psi2_val[0]; psi2<=psi2_val[1]; psi2+=psi2_val[2]) {
      counterP2+=1;
      fprintf(output_file,"%.2f %.1f %lf %.1f %.1f %d %.1f\n",e2, psi2,chi_min[counterE2][counterP2],alpha_minArray[counterE2][counterP2],zeta_minArray[counterE2][counterP2],phi0_minArray[counterE2][counterP2],psi0_minArray[counterE2][counterP2]);
    }
  }
  fclose(input_file);
  fclose(output_file);
  return 0;
}
