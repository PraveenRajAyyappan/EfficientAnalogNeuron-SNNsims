#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Macros
#define fi(x,FracBits) ((int)((x) * (1<<FracBits)))
#define mul_24_24_24(x,y) (((long)(x) * (long)(y)) >> 24)
#define mul_24_8_8(x,y) (((long)(x) * (long)(y)) >> 24)
#define mul_24_8_24(x,y) (((long)(x) * (long)(y)) >> 8)
#define mul_8_8_8(x,y) (((long)(x) * (long)(y)) >> 8)
#define mul_16_8_8(x,y) (((long)(x) * (long)(y)) >> 16)
#define mul_23_8_8(x,y) (((long)(x) * (long)(y)) >> 23)


#define Offset_a 2.8309
#define Offset_b -2.1341
#define Offset_c 4.4310




// Functions
void float_vector_to_fix(float vec[4], int fix_vec[4]) {          // Converts a floating point vector to a fixed point vector
  for (int i = 0; i < 4; i++) {
    fix_vec[i] = fi(vec[i],24);
  }
}

void float_array_to_fix(float array[4][8], int fix_array[4][8]) { // Converts a floating point array to a fixed point vector
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 8; j++)
      fix_array[i][j] = fi(array[i][j],23);
  }
}

float two_power_of_appx(int x){       // Fixed point radix 2 exponentiation (linear version)
/*
  int intr_x = ((-((x & 0x80000000)>>31)) & (0xFFFFFF80)) | ((x & 0x7f000000) >> 24);
  if(intr_x >= -6){
  int shift  = intr_x-16; 
  int frac_x = x & 0x00ffffff; // Input in Q8.24
  int Output;
  if(shift >= 0)
    Output = ((0x1000000+frac_x)<<shift); // Output in Q24.8
  else
    Output = ((0x1000000+frac_x)>>abs(shift)); // Output in Q24.8

  // Convert
    return(Output);
  } else{
    return(0);
  }
  */

  // Convert input (Q8.24) to floating point
  float float_x = (float) x / (1<<24);

  // Compute 2^x
  float float_out = pow(2.0,float_x);

  // Convert back to fix (Q24.8) 
  return(fi(float_out,8));
}

float two_power_of(int x){       // Fixed point radix 2 exponentiation (linear version)

  int intr_x = ((-((x & 0x80000000)>>31)) & (0xFFFFFF80)) | ((x & 0x7f000000) >> 24);
  if(intr_x >= -6){
  int shift  = intr_x-16; 
  int frac_x = x & 0x00ffffff; // Input in Q8.24
  int Output;
  if(shift >= 0)
    Output = ((0x1000000+frac_x)<<shift); // Output in Q24.8
  else
    Output = ((0x1000000+frac_x)>>abs(shift)); // Output in Q24.8

  // Convert
    return(Output);
  } else{
    return(0);


    
  }

}

void matrix_mult_2d(int array_1[8], int array[4][8], int array_final[4]){
  int Result;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 8; j++){
    array_final[i] += mul_23_8_8(array[i][j],array_1[j]);
    }
  }
}

float two_power_of_cubic(int x){       // Fixed point radix 2 exponentiation (linear version)

  int intr_x = ((-((x & 0x80000000)>>31)) & (0xFFFFFF80)) | ((x & 0x7f000000) >> 24);
  if(intr_x >= -6){
    int shift  = intr_x-16; 
    int frac_x = x & 0x00ffffff; // Input in Q8.24
    int Output;
    if(shift >= 0)
      Output = mul_24_8_8((0x1000000+frac_x)<<shift,(((mul_24_24_24(frac_x,frac_x) - frac_x)>>2) + 16777216)); // Output in Q24.8
    else
      Output = mul_24_8_8((0x1000000+frac_x)>>abs(shift),(((mul_24_24_24(frac_x,frac_x) - frac_x)>>2) + 16777216)); // Output in Q24.8

    // Convert
      return(Output);
  } else{
    return(0);
  }
}




int FH_ODE(int fix_y[4], int fix_Inp_curr,int fix_EXPVNA, int fix_EXPEK, int fix_EXPMEL, int fix_EXPVSAT, int fix_EXPMVGK, int Kp, int fix_TMAT[4][8], int fix_dydt[4]){
  // Compute IVec
  int fix_I_Vec[8] = {fix_Inp_curr, 
                    256-mul_24_8_8(fix_EXPMEL,two_power_of_cubic(fix_y[0])),//12,201326592
                    mul_8_8_8( two_power_of(mul_24_24_24(Kp,fix_y[3])),(256- mul_24_8_8(fix_EXPVNA,two_power_of(fix_y[0]))) ),
                    two_power_of(mul_24_24_24(Kp,fix_y[2])),
                    256-mul_24_8_8(fix_EXPVSAT,two_power_of(-fix_y[3]-390120603)),//24,402653184
                    two_power_of_cubic(fix_y[0]+mul_24_24_24(Kp,fix_y[1])) - mul_24_8_8( fix_EXPEK,two_power_of(mul_24_24_24(Kp,fix_y[1])) ),
                    256-mul_24_8_8(fix_EXPMVGK,two_power_of(fix_y[1])),
                    two_power_of(fix_y[2])-two_power_of(fix_y[3])
                    };



  matrix_mult_2d(fix_I_Vec,fix_TMAT,fix_dydt); // dydt in Q24.8

}

int main() {
  //////// Initialization (Only Needs to Run Once) ///////
  // Define Model Constants
  float y[4] = {-Offset_a,-Offset_b,-Offset_c,-Offset_c}; // This is -offoff in the MATLAB code
  int fix_Inp_curr = 0;
  //float EXPMEL = 1.016374731485013*pow(2,Offset_a-12.0);
  //float EXPVSAT = 2.320369728666329*pow(2,-Offset_c);
  //float EXPMVGK = 0.947133465162177*pow(2,Offset_b);

  float EXPMEL = 7.3306;
  float EXPVSAT= 6.6116;
  float EXPMVGK = 0.2279;
  float EXPVNA = 0.0182;
  float EXPEK = 0.1364;

  float T_MAT[4][8] = {{70.785082606739/1000, 20.2173858708657/1000,	1.32921811478727/1000,	2387.77324751470/1000,	-23893.5281204093/1000,	-9026.24819333978/1000,	135.171307831437/1000,	-232.108910489866/1000},
  {70.7850826067390/1000,	20.2173858708657/1000,	1.32921811478727/1000,	2387.77324751470/1000,	-23893.5281204093/1000,	-9026.24819333978/1000,	410.668585361385/1000,	-232.108910489866/1000},
  {69.6603580753146/1000,	19.8961460133584/1000,	1.30809778595160/1000,	4377.58428711028/1000,	-43804.8015540837/1000, -8882.82754034463/1000,	133.023532053490/1000,	425.533002564754/1000},
  {65.3065856956075/1000,	18.6526368875235/1000,	1.22634167432962/1000,	20587.9115563489/1000,	-206015.309127084/1000,	-8327.65081907309/1000,	124.709561300147/1000,	23636.4240515513/1000}};



  // Define simulation parameters
  float dt     = 0.001; // Units (ms)
  int   NSteps = 1e5;  // Total Number of steps to take in the simulation
  int PulseStartIdx = 20000;
  int PulseHeight = fi(100.0,8); // In pA

  // Convert to fix point
  int fix_TMAT[4][8] = {};
  float_array_to_fix(T_MAT, fix_TMAT);
  int fix_dt = fi(dt,24);
  int fix_EXPMEL = fi(EXPMEL,24);
  int fix_EXPVSAT= fi(EXPVSAT,24);
  int fix_EXPMVGK= fi(EXPMVGK,24);
  int fix_EXPVNA = fi(EXPVNA,24);
  int fix_EXPEK = fi(EXPEK,24);
  int Kp = fi(-0.75,24);
  int fix_y[4] = {};
  float_vector_to_fix(y, fix_y);

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
      fix_Inp_curr = PulseHeight;
    } else{
      fix_Inp_curr = 0;
    }

    // Run the actual ODE (Explicit Euler)
    int fix_dydt[4] = {0,0,0,0};
    FH_ODE(fix_y, fix_Inp_curr, fix_EXPVNA, fix_EXPEK, fix_EXPMEL, fix_EXPVSAT, fix_EXPMVGK, Kp, fix_TMAT, fix_dydt);
    
    // Take Euler Steps (y_next = y_current + dydt*dt)
    for(int j = 0; j < 4; j++){
      fix_y[j] += mul_24_8_24(fix_dydt[j],fix_dt);
    }

    // Print Vmem to file
    fprintf(fptr,"%f\n",(float) fix_y[0] / (1<<24));
  }
  
  // Close file
  fclose(fptr);

  return 0;
}