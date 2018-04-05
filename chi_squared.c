//
//  chi_squared.c
//
//  Created by Tjark Miener and Bruce Allen on 17.08.17.
//  Copyright (c) 2017 Tjark Miener and Bruce Allen. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//To compile this code on a Linux machine, the constant M_PI (actually a constant from math.h) has to be defined manually.
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct Data{
  int index;       //0,...,255
  float intensityI;
  float linearQ;
  float linearU;
  float linearL;
  float circularV;
  float psi;
  float psi_error;
};
struct Data list[256];

//Degrees to radians conversion.
double deg2rad(float val){
  return val*M_PI/180.0;
}
//Radians to degrees conversion.
float rad2deg(double val){
  return val*180.0/M_PI;
}

//All arguments in the head of the function are in radians, except the variable called 'index'. The returning position angle is declared in degrees.
double positionAngle(float alpha, float zeta, int index, float phi0, float psi0) {
  double phi_rad=index*2*M_PI/256.0;
  //Calculating the position angle with the rotating vector model (RVM).
  double x=-(sin(zeta)*cos(alpha)-cos(zeta)*sin(alpha)*cos(phi_rad-phi0));
  double y=sin(phi_rad-phi0)*sin(alpha);
  return rad2deg(atan2(y,x)) + psi0;
}
//main function
int main (int argc, char** argv) {
    
  //Initialization of helpful variables.
  int i,counter,counterA,counterB;
  float alpha_val[3],zeta_val[3],psi0_val[3];
  int phi0_val[3];
    
  //For reasons of parallelizing the code, two variables from the console mark the start and end of the alpha values. The third value of the array is the iteration size.
  alpha_val[0]=atof(argv[2]);
  alpha_val[1]=atof(argv[3]);
  alpha_val[2]=1.0;   //0.1 for the finer grid
    
  //Bin size for the heatmap, which the python code 'J2007heatmap.py' will visualize.
  int alpha_bins=(int)(alpha_val[1]-alpha_val[0]);    //(int)((alpha_val[1]-alpha_val[0])*10) for the finer grid
  int zeta_bins=181;  //201 for the finer grid
  double chi[alpha_bins][zeta_bins];
  int phi0_minArray[alpha_bins][zeta_bins];
  float psi0_minArray[alpha_bins][zeta_bins];
    
  FILE *input_file=fopen(argv[1],"r");
  int file_index=atoi(argv[2]);
  char str[35];
  //Detecting if the input file is the profile of 1500MHz or 2000MHz.
  if (argv[1][6]=='1') {
    sprintf(str, "J2007_1500MHz_Chi%d.txt", file_index);
    zeta_val[0]=0.0;      //95.0 for the finer grid
    zeta_val[1]=180.0;    //115.0 for the finer grid
    zeta_val[2]=1.0;      //0.1 for the finer grid
    phi0_val[0]=0;        //150 for the finer grid
    phi0_val[1]=360;      //250 for the finer grid
    phi0_val[2]=1;        //1 for the finer grid
    psi0_val[0]=-90.0;   //-30.0 for the finer grid
    psi0_val[1]=90.0;    //10.0 for the finer grid
    psi0_val[2]=1.0;      //0.1 for the finer grid
  } else if (argv[1][6]=='2') {
    sprintf(str, "J2007_2000MHz_Chi%d.txt", file_index);
    zeta_val[0]=0.0;      //100.0 for the finer grid
    zeta_val[1]=180.0;    //120.0 for the finer grid
    zeta_val[2]=1.0;      //0.1 for the finer grid
    phi0_val[0]=0;        //150 for the finer grid
    phi0_val[1]=360;      //250 for the finer grid
    phi0_val[2]=1;        //1 for the finer grid
    psi0_val[0]=-90.0;   //-20.0 for the finer grid
    psi0_val[1]=90.0;    //20.0 for the finer grid
    psi0_val[2]=1.0;      //0.1 for the finer grid
  } else {
    printf("Error! Inputfile %s is wrong!\n", argv[1]);
    exit(1);
  }
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
  double prefactor=1.0/(n-4.0);
    
  //Minimizing Chi^2 over the parameter space (alpha,zeta,phi0,psi0).
  alpha_val[1]-=1.0;      //0.1 for the finer grid
  counterA=-1;
  for (float alpha=alpha_val[0]; alpha<=alpha_val[1]; alpha+=alpha_val[2]) {
    double alpha_rad=deg2rad(alpha);
    counterA+=1;
    counterB=-1;
    for (float zeta=zeta_val[0]; zeta<=zeta_val[1]; zeta+=zeta_val[2]) {
      printf("Iteration: alpha=%.1f (%.1f-%.1f); zeta=%.1f\n",alpha,alpha_val[0],alpha_val[1],zeta);
      double zeta_rad=deg2rad(zeta);
      counterB+=1;
      for (int phi0=phi0_val[0]; phi0<=phi0_val[1]; phi0+=phi0_val[2]) {
        double phi0_rad=deg2rad(phi0);
        for (float psi0=psi0_val[0]; psi0<=psi0_val[1]; psi0+=psi0_val[2]) {
                    
          //Calculating the Chi^2 for one point of the parameter space.
          double psi_theo[256];
          double sum=0.0;
          for (int i=0; i<=255;i++) {
            if (list[i].psi!=0 && list[i].psi_error!=0) {
              psi_theo[i]=positionAngle(alpha_rad, zeta_rad, i, phi0_rad, psi0);
              sum+=(pow((list[i].psi-psi_theo[i])/list[i].psi_error,2));
            }
          }
            
          //Initializing Chi^2.
          if (phi0==phi0_val[0] && psi0==psi0_val[0]) {
            chi[counterA][counterB]=prefactor*sum;
            phi0_minArray[counterA][counterB]=phi0;
            psi0_minArray[counterA][counterB]=psi0;
              
          } else {
            //Minimizing Chi^2.
            if (chi[counterA][counterB] > prefactor*sum) {
              chi[counterA][counterB]=prefactor*sum;
              phi0_minArray[counterA][counterB]=phi0;
              psi0_minArray[counterA][counterB]=psi0;
            }
          }
        }
      }
    }
  }
    
  //Printing the minimized Chi^2 into the output file. The python code 'J2007heatmap.py' will read in this file and visualize the Chi^2 using a heatmap.
  counterA=-1;
  for (float alpha=alpha_val[0]; alpha<=alpha_val[1]; alpha+=alpha_val[2]) {
    counterA+=1;
    counterB=-1;
    for (float zeta=zeta_val[0]; zeta<=zeta_val[1]; zeta+=zeta_val[2]) {
      counterB+=1;
      fprintf(output_file,"%.1f %.1f %d %.1f %lf\n",alpha, zeta, phi0_minArray[counterA][counterB], psi0_minArray[counterA][counterB],chi[counterA][counterB]);
    }
  }
    
  fclose(input_file);
  fclose(output_file);
    
  return 0;
}
