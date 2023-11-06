#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// 12 2 2
#define Offset_a 2.8309
#define Offset_b -2.1341
#define Offset_c 4.4310

void matrix_mult_2d(float array_1[8], float array[4][8], float array_final[4]){

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 8; j++){
    array_final[i] = array_final[i]+(array[i][j]*array_1[j]);

    }

  }
  
}

int FH_ODE(float y[4], float Inp_curr,float EXPVNA,float EXPEK, float EXPMEL, float EXPVSAT, float EXPMVGK, float TMAT[4][8], float dydt[4]){
  float I_Vec[8] = {Inp_curr, 
                    1-(EXPMEL*(pow(2, y[0]))),
                    pow(2,(-0.75*y[3]))*(1-EXPVNA*pow(2,y[0])),
                    pow(2,(-0.75*y[2])),
                    1-(EXPVSAT*(pow(2, -y[3]-23.253))),
                    pow(2,(y[0]-(0.75*y[1])))- EXPEK*pow(2,-0.75*(y[1])),
                    1-(EXPMVGK*(pow(2, y[1]))),
                    pow(2,(y[2]))-pow(2,(y[3]))
                    };

  matrix_mult_2d(I_Vec,TMAT,dydt);
}

int main() {
  //////// Initialization (Only Needs to Run Once) ///////
  // Define Model Constants
  float y[4] = {-Offset_a,-Offset_b,-Offset_c,-Offset_c}; // This is -offoff in the MATLAB code
  float Inp_curr = 0;
  //float EXPMEL = 1.016374731485013*pow(2,Offset_a-12.0);
  //float EXPVSAT = 2.320369728666329*pow(2,-Offset_c);
  //float EXPMVGK = 0.947133465162177*pow(2,Offset_b);

  float EXPMEL = 7.3306;
  float EXPVSAT= 6.6116;
  float EXPMVGK = 0.2279;
  float EXPVNA = 0.0182;
  float EXPEK = 0.1364;


 /* float T_MAT[4][8] = {{72.1347520444482/1000.0,	0.417041544961956/1000.0,	3.25462226097323*pow(2,-0.75*Offset_c)/1000.0,	9453.65999053692*pow(2,-0.75*Offset_c)/1000.0,	-9116.17833656249/1000.0,	-0.4284458685*pow(2,-0.75*Offset_b+Offset_a)/1000.0,	141.943709361294/1000.0,	0.0},
  {72.1347520444482/1000.0,	0.417041544961956/1000.0,	3.25462226097323*pow(2,-0.75*Offset_c)/1000.0,	9453.65999053692*pow(2,-0.75*Offset_c)/1000.0,	-9116.17833656249/1000.0,	-0.4284458685*pow(2,-0.75*Offset_b+Offset_a)/1000.0,	425.831128083882/1000.0,	0.0},
  {72.1347520444482/1000.0,	0.417041544961956/1000.0,  3.25462226097323*pow(2,-0.75*Offset_c)/1000.0,	17331.7099826510*pow(2,-0.75*Offset_c)/1000.0,	-16712.9936170312/1000.0, -0.4284458685*pow(2,-0.75*Offset_b+Offset_a)/1000.0,	141.943709361294/1000.0,	0.0},
  {72.1347520444482/1000.0,	0.417041544961956/1000.0, 	3.25462226097323*pow(2,-0.75*Offset_c)/1000.0,	80356.1099195639*pow(2,-0.75*Offset_c)/1000.0,	-77487.5158607812/1000.0,	-0.4284458685*pow(2,-0.75*Offset_b+Offset_a)/1000.0,	141.943709361294/1000.0,	1837.85855789555*pow(2,Offset_c)/1000.0}};
*/

  /*float T_MAT[4][8] = {
    //{1, 0, 0, 0, 0, 0, 0, 0},
  {0.0721, 0.0004, 0.0012, 3.3424, -9.1162, -4.9636, 0.1419, 0},
  {0.0721, 0.0004, 0.0012, 3.3424, -9.1162, -4.9636, 0.4258,  0},
  {0.0721, 0.0004, 0.0012, 6.1277, -16.7130, -4.9636, 0.1419, 0},
  {0.0721, 0.0004, 0.0012, 28.4102, -77.4875, -4.9636, 0.1419 ,7.3514}};*/

  float T_MAT[4][8] = {{70.785082606739/1000, 20.2173858708657/1000,	1.32921811478727/1000,	2387.77324751470/1000,	-23893.5281204093/1000,	-9026.24819333978/1000,	135.171307831437/1000,	-232.108910489866/1000},
  {70.7850826067390/1000,	20.2173858708657/1000,	1.32921811478727/1000,	2387.77324751470/1000,	-23893.5281204093/1000,	-9026.24819333978/1000,	410.668585361385/1000,	-232.108910489866/1000},
  {69.6603580753146/1000,	19.8961460133584/1000,	1.30809778595160/1000,	4377.58428711028/1000,	-43804.8015540837/1000, -8882.82754034463/1000,	133.023532053490/1000,	425.533002564754/1000},
  {65.3065856956075/1000,	18.6526368875235/1000,	1.22634167432962/1000,	20587.9115563489/1000,	-206015.309127084/1000,	-8327.65081907309/1000,	124.709561300147/1000,	23636.4240515513/1000}};


  // Define simulation parameters
  float dt     = 0.001; // Units (ms)
  int   NSteps = 1e5;  // Total Number of steps to take in the simulation
  int PulseStartIdx = 20000;
  float PulseHeight = 350.0; // In pA


  // Open a file for writing
  FILE *fptr;
  fptr = fopen("Output.txt","w");
  if(fptr == NULL){
      printf("File Error!");   
      exit(1);             
   }


  for(int i = 0; i < NSteps; i++){
    // Compute the current input current
    if(i >= PulseStartIdx){
      Inp_curr = PulseHeight;
    } else{
      Inp_curr = 0;
    }

    // Run the actual ODE (Explicit Euler)
    float dydt[4] = {0,0,0,0};
    FH_ODE(y, Inp_curr,  EXPVNA, EXPEK, EXPMEL, EXPVSAT, EXPMVGK, T_MAT, dydt);
    

    // Take Euler Steps (y_next = y_current + dydt*dt)
    for(int j = 0; j < 4; j++){
      y[j] += dydt[j]*dt;
    }

    // Print Vmem to file
    //Print_1dArray(fixed_point_dydt,4);
    fprintf(fptr,"%f\n",y[0]);

  }
  

  // Close file
  fclose(fptr);

  return 0;
}