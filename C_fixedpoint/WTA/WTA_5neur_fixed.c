#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Macros
#define fi(x,FracBits) ((int)((x) * (1<<FracBits)))
#define mul_24_24_24(x,y) (((long)(x) * (long)(y)) >> 24)
#define mul_24_8_24(x,y) (((long)(x) * (long)(y)) >> 8)
#define mul_24_8_8(x,y) (((long)(x) * (long)(y)) >> 24)
#define mul_16_24_24(x,y) (((long)(x) * (long)(y)) >> 16)
#define mul_16_24_16(x,y) (((long)(x) * (long)(y)) >> 24)
#define mul_8_8_8(x,y) (((long)(x) * (long)(y)) >> 8)
#define mul_22_8_8(x,y) (((long)(x) * (long)(y)) >> 22)
#define mul_23_8_8(x,y) (((long)(x) * (long)(y)) >> 22)


#define VDD_Q16 0x009044FD             // VDD in Q16.16
#define MinusSpikeThresh fi(-1.4,24)  // 3.2 in Q8.24
#define Six 0x06000000             // 6.0 in Q8.24
#define VDD 144.2695               // In floating point
#define NeuronPopulation 5
#define Kp 0.75                    // Kappa for p-channel
#define UT 0.025                   // Thermal voltage
#define EffCoupling -10066329     // Combination of both the CG-gate and gate-channel coupling (negative here).
#define Tri_Top_Off 8509364       // Resting offset of the control gate of the excitatory synapse. Q16.16
#define Tri_Bot_Off 945487         // Resting offset of the control gate of the inhibitory synapse. Q16.16
#define TriScaling  0xCCCCC         // Scaling of the triangle generator outputs. Q8.24
/*
#define Offset_a 2.8366
#define Offset_b -2.0972
#define Offset_c 4.4344
*/
#define Offset_a 1
#define Offset_b -1
#define Offset_c 4

#define OFF 21

//////////////////// TRIANGLE GENERATOR ////////////////////
void float_array2_to_fix(float array[2][2], int fix_array[2][2]) { // Converts a floating point array to a fixed point vector
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++)
      fix_array[i][j] = fi(array[i][j],16);
  }
}

void float_synarray_to_fix(float array[NeuronPopulation][NeuronPopulation*2], int fix_array[NeuronPopulation][NeuronPopulation*2]) { // Converts a floating point array to a fixed point vector
  for (int i = 0; i < NeuronPopulation; i++) {
    for (int j = 0; j < NeuronPopulation*2; j++)
      fix_array[i][j] = fi(array[i][j],24);
  }
}

int Tri_gen(int fix_dydt[2], int fix_y[2],  int fix_AMAT[2][2], int fix_Inp){
  // Compute Derivative
  int Gamma = ((~(fix_Inp+MinusSpikeThresh))>>31);
  int NGamma = ~Gamma;
  fix_dydt[0] = (fix_AMAT[0][0]&Gamma)|(fix_AMAT[0][1]&NGamma);
  fix_dydt[1] = (fix_AMAT[1][0]&Gamma)|(fix_AMAT[1][1]&NGamma);

  // Saturation Logic
  fix_dydt[0] &= ~((((~(fix_y[0]-VDD_Q16))>>31) & ((~fix_dydt[0])>>31)) | (((~(-fix_y[0]))>>31) & ((~(-fix_dydt[0]))>>31)));
  fix_dydt[1] &= ~((((~(fix_y[1]-VDD_Q16))>>31) & ((~fix_dydt[1])>>31)) | (((~(-fix_y[1]))>>31) & ((~(-fix_dydt[1]))>>31)));
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

  // Convert input (Q8.24) to floating point
  float float_x = (float) x / (1<<24);

  // Compute 2^x
  float float_out = pow(2.0,float_x);

  // Convert back to fix (Q24.8) 
  return(fi(float_out,8));
}


float two_power_of_appx(int x){       // Fixed point radix 2 exponentiation (linear version)
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
                    mul_8_8_8( two_power_of_appx(mul_24_24_24(fix_Kp,fix_y[3])),(256- mul_24_8_8(fix_EXPVNA,two_power_of_appx(fix_y[0]))) ),
                    two_power_of_appx(mul_24_24_24(fix_Kp,fix_y[2])),
                    256-mul_24_8_8(fix_EXPVSAT,two_power_of_appx(-fix_y[3]-390120603)),//24,402653184
                    two_power_of_cubic(fix_y[0]+mul_24_24_24(fix_Kp,fix_y[1])) - mul_24_8_8( fix_EXPEK,two_power_of_appx(mul_24_24_24(fix_Kp,fix_y[1])) ),
                    256-mul_24_8_8(fix_EXPMVGK,two_power_of_appx(fix_y[1])),
                    two_power_of_appx(fix_y[2])-two_power_of_appx(fix_y[3])
                    };

  // Do the matrix multiplication
  matrix_mult_2d(fix_I_Vec,fix_TMAT,fix_dydt); // dydt in Q24.8
}

///////////////////////// END OF NEURON PART /////////////////////////////////
void compute_synapse_mat(float weights[NeuronPopulation][NeuronPopulation*2], float Synapse_Mat_fixed_term, float synapse_mat[NeuronPopulation][NeuronPopulation*2]){

  for (int i = 0; i < NeuronPopulation; i++) {
    for (int j = 0; j < (NeuronPopulation*2); j++){
      if(weights[i][j]<2.5){
        synapse_mat[i][j] = Synapse_Mat_fixed_term * pow(2,(((-Kp*weights[i][j])/(UT*log(2)))+9.9658));
      }
      else
        synapse_mat[i][j] = 0;
    }
  }
}

void initialize_y(float y[NeuronPopulation*6]){

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

void update_Vmem(int y[NeuronPopulation*6], int Vmem[NeuronPopulation]){
  int temp=0;
  for(int i=0; i<(NeuronPopulation*6);i++){
    if (i == (NeuronPopulation*2 +temp))
    {
      Vmem[temp >> 2]= y[i];
      temp += 4;
    }
  }
}

void intialize_top_bot_idx(int top_idx[2*NeuronPopulation], int bot_idx[2*NeuronPopulation]){
  for(int i=0; i<(NeuronPopulation*2);i++){
    if(i%2==0){
      top_idx[i]=0xFFFFFFFF;
      bot_idx[i]=0x0;
    }
    else{
      top_idx[i]=0x0;
      bot_idx[i]=0xFFFFFFFF;
    }
  }
}

void InputCurrent_compute1(int array_1[2*NeuronPopulation], int array[NeuronPopulation][NeuronPopulation*2], int array_final[NeuronPopulation], int Inp_current[NeuronPopulation]){
  for (int i = 0; i < NeuronPopulation; i++) {
    array_final[i] = 0;
    for (int j = 0; j < 2*NeuronPopulation; j++){
      array_final[i] += mul_24_8_8(array[i][j],array_1[j]);
    }
    array_final[i] += Inp_current[i];
  }
}

void InputCurrent_compute(int fix_Synapse_Mat[NeuronPopulation][NeuronPopulation*2], int fix_Scale_Mat[NeuronPopulation][NeuronPopulation*2], int Triangle_scaled[2*NeuronPopulation], int top_idx[2*NeuronPopulation], int bot_idx[2*NeuronPopulation], int fix_Inp_current[NeuronPopulation], int InputCurrent[NeuronPopulation]){

  int fix_temp[NeuronPopulation][NeuronPopulation*2], fix_temp1[NeuronPopulation*2];
  for (int i = 0; i < NeuronPopulation; i++) {
    for (int j = 0; j < 2*NeuronPopulation; j++){
      fix_temp[i][j] = mul_24_8_24(fix_Synapse_Mat[i][j],fix_Scale_Mat[i][j]); // Synapse in Q8.24, Scale in Q24.8
    }
  }

  int TopidxTemp, BotidxTemp;
  for(int i=0; i< 2*NeuronPopulation; i++){
    //fix_temp1[i] = two_power_of(mul_16_24_24(EffCoupling,Triangle_scaled[i])  + (0x56000000&top_idx[i]) + (0xFE000000&bot_idx[i]));
    //problem
    fix_temp1[i] = two_power_of_appx(mul_16_24_24(EffCoupling,Triangle_scaled[i])  + (0x56000000&top_idx[i]) + (fi(OFF-2.0,24)&bot_idx[i]));
    
  }

  // Do the main computation
  InputCurrent_compute1(fix_temp1, fix_temp, InputCurrent,fix_Inp_current);
}

int Network_ODE(int fix_dydt[NeuronPopulation*6], int y[NeuronPopulation*6], int fix_Vmem[NeuronPopulation], int fix_Synapse_Mat[NeuronPopulation][NeuronPopulation*2], int fix_Inp_current[NeuronPopulation], int fix_Kp, int fix_EXPMEL, int fix_EXPVSAT, int fix_EXPMVGK, int fix_EXPVNA , int fix_EXPEK , int fix_AMAT[2][2], int fix_T_MAT_neu[4][8] ){
  ///////// Set Up Some Indexing ///////
  int top_idx[2*NeuronPopulation], bot_idx[2*NeuronPopulation];
  intialize_top_bot_idx(top_idx,bot_idx);

  ///////// Triangle Pulse Scaling ///////
  int fix_Triangle_scaled[2*NeuronPopulation];
  int TopMask, BotBask;
  float FloatingTri;
  for (int i = 0; i < 2*NeuronPopulation; i++)
  {
    // Scale triangle
    fix_Triangle_scaled[i] = mul_16_24_16(y[i],TriScaling) + (top_idx[i]&Tri_Top_Off) + (bot_idx[i]&Tri_Bot_Off);
  }
  
  ///////// Scaling Matrix Computation ///////
  int fix_Scale_Mat[NeuronPopulation][NeuronPopulation*2];
  for (int i = 0; i < NeuronPopulation; i++)
  {
    for (int j = 0; j < 2*NeuronPopulation; j++)
    {
      if(top_idx[j])
        fix_Scale_Mat[i][j] = 0x58;
      else{
        //fix_Scale_Mat[i][j] = -two_power_of(fix_Vmem[i]+0xD348E8A);
        //problem
        fix_Scale_Mat[i][j] = -two_power_of_appx(fix_Vmem[i]+fi(13.2053-OFF,24));
      }
    }
  }

  ///////// Input Current Computation ///////
  int fix_InputCurrent[NeuronPopulation];
  InputCurrent_compute(fix_Synapse_Mat, fix_Scale_Mat, fix_Triangle_scaled, top_idx, bot_idx, fix_Inp_current, fix_InputCurrent);

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
  int   NSteps = 6e4;  // Total Number of steps to take in the simulation
  // For 1 neuron 20000-80000  - 0.06s/60ms
  // For 3 neuron 20000-120000 - 0.1s/100ms
  int PulseStartIdx = 1e4; 
  int PulseendIdx = 5e4;
  int fix_Inp_current[NeuronPopulation] = {0};

  int fix_PulseHeight1 = fi(85.0,8); 
  int fix_PulseHeight2 = fi(55.0,8); 
  int fix_PulseHeight3 = fi(45.0,8); 
  int fix_PulseHeight4 = fi(35.0,8); 


  
  // Define Network Connectivity Via the Synapse Weight Matrix
  // float weights[NeuronPopulation][NeuronPopulation*2]={2.5,2.5};
  float Syn = 0.26, Ih=0.23;
  float weights[NeuronPopulation][NeuronPopulation*2]={{2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, Ih},
         {2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, Ih},
         {2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, Ih},
         {2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, Ih},
         {Syn, 2.5, Syn, 2.5, Syn, 2.5, Syn, 2.5, 2.5, 2.5}};
  float Synapse_Mat_fixed_term = 2.1425, Synapse_Mat[NeuronPopulation][NeuronPopulation*2]={0};
  compute_synapse_mat(weights,Synapse_Mat_fixed_term,Synapse_Mat);

  // Intialize values for y, and the Vmem array for all neurons
  float y[NeuronPopulation*6];
  initialize_y(y);
  int fix_Vmem[NeuronPopulation];
  int fix_dydt[6*NeuronPopulation];

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
  int fix_Synapse_Mat[NeuronPopulation][NeuronPopulation*2];
  float_synarray_to_fix(Synapse_Mat, fix_Synapse_Mat);
  int fix_y[6*NeuronPopulation];
  for(int j = 0; j < 2*NeuronPopulation; j++)
    fix_y[j] = fi(y[j],16);
  for(int j = 2*NeuronPopulation; j < 6*NeuronPopulation; j++)
    fix_y[j] = fi(y[j],24);
  update_Vmem(fix_y,fix_Vmem);

  // Open a file for writing
  FILE *fptr,*fptr1;
  fptr = fopen("WTA.txt","w");
  if(fptr == NULL){
      printf("File Error!");   
      exit(1);             
   }

  for(int i = 0; i < NSteps; i++){
    // Compute the external input current
    if(i >= PulseStartIdx && i <= PulseendIdx){
      fix_Inp_current[0] = fix_PulseHeight1;
      fix_Inp_current[1] = fix_PulseHeight2;
      fix_Inp_current[2] = fix_PulseHeight3;
      fix_Inp_current[3] = fix_PulseHeight4;


    }
    else{
      fix_Inp_current[0] = 0;
      fix_Inp_current[1] = 0;
      fix_Inp_current[2] = 0;
      fix_Inp_current[3] = 0;


    }

    // Run the actual ODE (Explicit Euler)
    Network_ODE(fix_dydt, fix_y, fix_Vmem, fix_Synapse_Mat, fix_Inp_current, fix_Kp, fix_EXPMEL, fix_EXPVSAT, fix_EXPMVGK, fix_EXPVNA, fix_EXPEK, fix_AMAT, fix_T_MAT_neu);
    
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
    fprintf(fptr,"%f %f %f %f %f\n",(float) fix_y[10]/(1<<24),(float) fix_y[14]/(1<<24),(float) fix_y[18]/(1<<24),(float) fix_y[22]/(1<<24),(float) fix_y[26]/(1<<24));
  }
  
  // Close file
  fclose(fptr);
  printf("Done.\n");
  return 0;
}