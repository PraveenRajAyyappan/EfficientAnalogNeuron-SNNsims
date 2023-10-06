#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Macros
#define fi(x,FracBits) ((int)((x) * (1<<FracBits)))
#define mul_24_24_24(x,y) (((long)(x) * (long)(y)) >> 24)
#define mul_24_8_24(x,y) (((long)(x) * (long)(y)) >> 8)
#define mul_24_8_8(x,y) (((long)(x) * (long)(y)) >> 24)
#define mul_16_24_24(x,y) (((long)(x) * (long)(y)) >> 16)
#define mul_16_24_16(x,y) (((long)(x) * (long)(y)) >> 24)
#define mul_8_8_8(x,y) (((long)(x) * (long)(y)) >> 8)
#define mul_22_8_8(x,y) (((long)(x) * (long)(y)) >> 22)
#define mul_23_8_8(x,y) (((long)(x) * (long)(y)) >> 23)


#define VDD_Q16 0x009044FD             // VDD in Q16.16
#define MinusSpikeThresh -24397427  // 3.2 in Q8.24
#define Six 0x06000000             // 6.0 in Q8.24
#define VDD 144.2695               // In floating point
#define NeuronPopulation 100000
#define SimSteps 240e6
#define Kp 0.75                    // Kappa for p-channel
#define UT 0.025                   // Thermal voltage
#define EffCoupling -10066329     // Combination of both the CG-gate and gate-channel coupling (negative here).
#define Tri_Top_Off 8509364       // Resting offset of the control gate of the excitatory synapse. Q16.16
#define Tri_Bot_Off 945487         // Resting offset of the control gate of the inhibitory synapse. Q16.16
#define TriScaling  0xCCCCC         // Scaling of the triangle generator outputs. Q8.24

#define Offset_a 2.8309
#define Offset_b -2.1341
#define Offset_c 4.4310
#define OFF 21


//////////////////// TRIANGLE GENERATOR ////////////////////
void float_array2_to_fix(float array[2][2], int fix_array[2][2]) { // Converts a floating point array to a fixed point vector
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++)
      fix_array[i][j] = fi(array[i][j],16);
  }
}

int Tri_gen(int fix_dydt[2], int fix_y[2],  int fix_AMAT[2][2], int fix_Inp){
  // Compute Derivative
  int Gamma = ((~(fix_Inp+MinusSpikeThresh))>>31);
  int NGamma = ~Gamma;
  fix_dydt[0] = (fix_AMAT[0][0]&Gamma)|(fix_AMAT[0][1]&NGamma);

  // Saturation Logic
  fix_dydt[0] &= ~((((~(fix_y[0]-VDD_Q16))>>31) & ((~fix_dydt[0])>>31)) | (((~(-fix_y[0]))>>31) & ((~(-fix_dydt[0]))>>31)));
}

/////////////////// NEURON PART /////////////////////////////
void float_vector_to_fix(float vec[4], int fix_vec[4]) {          // Converts a floating point vector to a fixed point vector
  for (int i = 0; i < 4; i++) {
    fix_vec[i] = fi(vec[i],24);
  }
}

void matrix_mult_2d(int array_1[8], int array[4][8], int array_final[4]){
  for (int i = 0; i < 4; i++) {
    array_final[i] = 0;
    for (int j = 0; j < 8; j++){
    array_final[i] += mul_23_8_8(array[i][j],array_1[j]);
    }
  }
}

void float_array_to_fix(float array[4][8], int fix_array[4][8]) { // Converts a floating point array to a fixed point vector
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 8; j++)
      fix_array[i][j] = fi(array[i][j],23);
  }
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

int FH_ODE(int fix_y[4], int fix_Inp_curr, int fix_EXPVNA, int fix_EXPEK, int fix_EXPMEL, int fix_EXPVSAT, int fix_EXPMVGK, int fix_Kp, int fix_TMAT[4][8], int fix_dydt[4]){
  // Compute IVec  
  int fix_I_Vec[8] = {fix_Inp_curr, 
                    256-mul_24_8_8(fix_EXPMEL,two_power_of_cubic(fix_y[0])),//12,201326592
                    mul_8_8_8( two_power_of(mul_24_24_24(fix_Kp,fix_y[3])),(256- mul_24_8_8(fix_EXPVNA,two_power_of(fix_y[0]))) ),
                    two_power_of(mul_24_24_24(fix_Kp,fix_y[2])),
                    256-mul_24_8_8(fix_EXPVSAT,two_power_of(-fix_y[3]-390120603)),//24,402653184
                    two_power_of_cubic(fix_y[0]+mul_24_24_24(fix_Kp,fix_y[1])) - mul_24_8_8( fix_EXPEK,two_power_of(mul_24_24_24(fix_Kp,fix_y[1])) ),
                    256-mul_24_8_8(fix_EXPMVGK,two_power_of(fix_y[1])),
                    two_power_of(fix_y[2])-two_power_of(fix_y[3])
                    };

  // Do the matrix multiplication
  matrix_mult_2d(fix_I_Vec,fix_TMAT,fix_dydt); // dydt in Q24.8
}

///////////////////////// END OF NEURON PART /////////////////////////////////
void initialize_y(float *y){

  int temp=0;
  for(int i=0; i<(NeuronPopulation*6);i++){

    if(i<NeuronPopulation*2){
      y[i] = VDD;
    }

    else if (i >= (NeuronPopulation*2))
    {
      switch ((i-(NeuronPopulation*2)) % 4) {
        case 0:
            y[i] = -Offset_a;  // if i mod 4 equals 0
            break;
        case 1:
            y[i] = -Offset_b;   // if i mod 4 equals 1
            break;
        case 2:
            y[i] = -Offset_c;   // if i mod 4 equals 2
            break;
        case 3:
            y[i] = -Offset_c;   // if i mod 4 equals 3
            break;
        default:
            break;
    }
    }
  }
}

void update_Vmem(int *y, int Vmem[NeuronPopulation]){
  int temp=0;
  for(int i=0; i<(NeuronPopulation*6);i++){
    if (i == (NeuronPopulation*2 +temp))
    {
      Vmem[temp >> 2]= y[i];
      temp += 4;
    }
  }
}

void InputCurrent_compute(int CurrPrefix, int Triangle_scaled[2*NeuronPopulation], int fix_Inp_current[NeuronPopulation], int InputCurrent[NeuronPopulation]){
  int fix_temp1;
  fix_temp1 = two_power_of(mul_16_24_24(EffCoupling,Triangle_scaled[2*(NeuronPopulation-1)]) + 0x56000000);
  InputCurrent[0] = fix_Inp_current[0] + mul_24_8_8(CurrPrefix,fix_temp1);
  for(int i=1; i< NeuronPopulation; i ++){
    fix_temp1 = two_power_of(mul_16_24_24(EffCoupling,Triangle_scaled[2*(i-1)]) + 0x56000000);
    InputCurrent[i] = fix_Inp_current[i] + mul_24_8_8(CurrPrefix,fix_temp1);
  }
}

int Network_ODE(int *fix_dydt, int *y, int fix_Vmem[NeuronPopulation], int CurrPrefix, int fix_Inp_current[NeuronPopulation], int fix_Kp, int fix_EXPMEL, int fix_EXPVSAT, int fix_EXPMVGK, int fix_EXPVNA , int fix_EXPEK , int fix_AMAT[2][2], int fix_T_MAT_neu[4][8] ){
  ///////// Triangle Pulse Scaling ///////
  int fix_Triangle_scaled[2*NeuronPopulation] = {0};
  float FloatingTri;
  for (int i = 0; i < 2*NeuronPopulation; i +=2)
  {
    // Scale triangle
    fix_Triangle_scaled[i] = mul_16_24_16(y[i],TriScaling) + Tri_Top_Off;
  }

  ///////// Input Current Computation ///////
  int fix_InputCurrent[NeuronPopulation];
  InputCurrent_compute(CurrPrefix, fix_Triangle_scaled, fix_Inp_current, fix_InputCurrent);

  ///////// Triangle Gen Part ///////
  int fix_dydt_tri[2],fix_y_tri[2];
  for(int i=0; i<NeuronPopulation; i++){
    // Extract triangle gen y components
    for (int a = 0; a < 2; a++) 
        fix_y_tri[a] = y[a+(2*i)];

    // Compute triangle gen derivative
    Tri_gen(fix_dydt_tri,fix_y_tri,fix_AMAT,fix_Vmem[i]);

    // Update dydt for triangle gen dydt components
    for (int a = 0; a < 2; a++) 
        fix_dydt[a+(2*i)] = fix_dydt_tri[a];
  }

  ///////// Neuron Part ///////
  int fix_y_neur[4];
  for(int i=0; i<NeuronPopulation; i++){
    // Extract neuron y components
    for (int a = 0; a < 4; a++)
        fix_y_neur[a] = y[a+(2*NeuronPopulation)+(4*i)];
    int fix_dydt_neur[4];

    // Update triangle gen derivative 
    FH_ODE(fix_y_neur, fix_InputCurrent[i], fix_EXPVNA, fix_EXPEK, fix_EXPMEL, fix_EXPVSAT, fix_EXPMVGK, fix_Kp, fix_T_MAT_neu, fix_dydt_neur);

    // Update neuron dydt components
    for (int a = 0; a < 4; a++) 
        fix_dydt[a+(2*NeuronPopulation)+(4*i)] = fix_dydt_neur[a];
  }
}

int main() {
  // Define simulation parameters
  float dt     = 0.001; // Units (ms)
  int   NSteps = SimSteps;  // Total Number of steps to take in the simulation
  // For 1 neuron 20000-80000  - 0.06s/60ms
  // For 3 neuron 20000-120000 - 0.1s/100ms
  int PulseStartIdx = 70000; 
  int PulseendIdx = 160000;
  int fix_Inp_current[NeuronPopulation] = {0};
  int fix_PulseHeight = fi(80.0,8); 
  int PrintThresh = fi(4.0,24);
  
  // Define Network Connectivity Via the Synapse Weight Matrix
  // float weights[NeuronPopulation][NeuronPopulation*2]={2.5,2.5};
  float Syn = 0.26;
  float Synapse_Mat_fixed_term = 2.1425;

  // Intialize values for y, and the Vmem array for all neurons


  float *y = (float *)malloc(NeuronPopulation*6 * sizeof(int));
  initialize_y(y);

  int fix_Vmem[NeuronPopulation];
  int *fix_dydt = (int *)malloc(6*NeuronPopulation * sizeof(int));
 

  // Define other static variables and transition matrices
  float EXPMEL = 7.3306;
  float EXPVSAT= 6.6116;
  float EXPMVGK = 0.2279;
  float EXPVNA = 0.0182;
  float EXPEK = 0.1364;

  float T_MAT_neu[4][8] = {{70.785082606739/1000, 20.2173858708657/1000,	1.32921811478727/1000,	2387.77324751470/1000,	-23893.5281204093/1000,	-9026.24819333978/1000,	135.171307831437/1000,	-232.108910489866/1000},
  {70.7850826067390/1000,	20.2173858708657/1000,	1.32921811478727/1000,	2387.77324751470/1000,	-23893.5281204093/1000,	-9026.24819333978/1000,	410.668585361385/1000,	-232.108910489866/1000},
  {69.6603580753146/1000,	19.8961460133584/1000,	1.30809778595160/1000,	4377.58428711028/1000,	-43804.8015540837/1000, -8882.82754034463/1000,	133.023532053490/1000,	425.533002564754/1000},
  {65.3065856956075/1000,	18.6526368875235/1000,	1.22634167432962/1000,	20587.9115563489/1000,	-206015.309127084/1000,	-8327.65081907309/1000,	124.709561300147/1000,	23636.4240515513/1000}};

  float T_MAT[2][2] = {{-71.2378,13.4597},{-272.0602,3.4893}};

  // Convert static constants to fixed point
  int fix_dt = fi(dt,24);
  int fix_Kp = fi(-Kp,24);
  int fix_EXPMEL = fi(EXPMEL,24);
  int fix_EXPVSAT = fi(EXPVSAT,24);
  int fix_EXPMVGK = fi(EXPMVGK,24);
  int fix_EXPVNA = fi(EXPVNA,24);
  int fix_EXPEK = fi(EXPEK,24);


  int fix_AMAT[2][2];
  float_array2_to_fix(T_MAT, fix_AMAT);
  int fix_T_MAT_neu[4][8];
  float_array_to_fix(T_MAT_neu, fix_T_MAT_neu);
  int CurrPrefix = mul_24_8_24(fi(Synapse_Mat_fixed_term * pow(2,(((-Kp*Syn)/(UT*log(2)))+9.9658)),24),0x58);

  int *fix_y = (int *)malloc(6*NeuronPopulation * sizeof(int));
  for(int j = 0; j < 2*NeuronPopulation; j++)
    fix_y[j] = fi(y[j],16);
  for(int j = 2*NeuronPopulation; j < 6*NeuronPopulation; j++)
    fix_y[j] = fi(y[j],24);
  update_Vmem(fix_y,fix_Vmem);

  // Open a file for writing
  char CurrStr[50];
  char Bigstr[100000] = "";
  int Offset;

  FILE *fptr,*fptr1;
  fptr = fopen("Synfire_100K.txt","w");
  if(fptr == NULL){
      printf("File Error!");   
      exit(1);             
   }

  clock_t begin = clock();
  for(int i = 0; i < NSteps; i++){
    // Compute the external input current
    if(i >= PulseStartIdx && i <= PulseendIdx){
      fix_Inp_current[0] = fix_PulseHeight;
    }
    else{
      fix_Inp_current[0] = 0;
    }

    // Run the actual ODE (Explicit Euler)
    Network_ODE(fix_dydt, fix_y, fix_Vmem, CurrPrefix, fix_Inp_current, fix_Kp, fix_EXPMEL, fix_EXPVSAT, fix_EXPMVGK, fix_EXPVNA, fix_EXPEK, fix_AMAT, fix_T_MAT_neu);
    
    // Take Euler Steps (y_next = y_current + dydt*dt)
    for(int j = 0; j < 6*NeuronPopulation; j++){
        if(j<2*NeuronPopulation){ // Triangle gen state
            fix_y[j] += mul_16_24_16(fix_dydt[j],fix_dt);
        }
        else{                     // Neuron state
            fix_y[j] += mul_24_8_24(fix_dydt[j],fix_dt);
        }
    }
    update_Vmem(fix_y,fix_Vmem);
    
    // Choose the y accordingly to print to the file.
    if (!(i%1000)){ // Print once every ms
    printf("%d ms\n",(int) i/1000);
    
    Offset = ((int)(0.0004*i))<<2;
    //for(int j=((2*NeuronPopulation)+Offset-1000); j<((2*NeuronPopulation)+Offset); j +=4){
    for(int j=(2*NeuronPopulation); j<(NeuronPopulation*6); j +=4){
      if(fix_y[j] >= PrintThresh){
        //fprintf(fptr,"%d ",(j-2*NeuronPopulation)>>2);
        sprintf(CurrStr,"%d ",(j-2*NeuronPopulation)>>2);
        //printf("%s",CurrStr);
        strcat(Bigstr,CurrStr);
      }
    }
    //fprintf(fptr,"\n");
    strcat(Bigstr,"\n");
    //printf("\n");

    // Print rarely
    if(!(i%100000)){
        printf("Now Printing!\n");
        fprintf(fptr,"%s",Bigstr);
        //fwrite(Bigstr , 1 , sizeof(Bigstr) , fptr);
        memset(Bigstr, 0, sizeof(Bigstr));
    }
    }
  }
  
  // Close file
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("%f\n",time_spent);
  fclose(fptr);
  printf("Done.\n");
  free(y);
  free(fix_y);
  free(fix_dydt);
  return 0;
}


