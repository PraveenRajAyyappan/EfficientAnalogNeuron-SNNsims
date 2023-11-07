#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#define VDD 144.2695
#define NeuronPopulation 4
#define Kp 0.75
#define UT 0.025


// 12 2 2
#define Offset_a 2.8309
#define Offset_b -2.1341
#define Offset_c 4.4310


//////////////////// TRIANGLE GENERATOR ////////////////////

void matrix_mult_2d(float array_1[2], float array[2][2], float array_final[2]){

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++){
    array_final[i] = array_final[i]+(array[i][j]*array_1[j]);

    }

  }
  
}


int Tri_gen(float dydt[2], float y[2],  float TMAT[2][2], float Inp_voltage){

  float zeta_mat[2] = {(Inp_voltage>=3),((Inp_voltage>=3)-1)};


  matrix_mult_2d(zeta_mat,TMAT,dydt);



  for(int j=0; j<2; j++){

    if ((y[j] >= VDD) && (dydt[j] > 0)){
      dydt[j] = 0;

    }
      
    else if ((y[j] <= 0) && (dydt[j] < 0)){
      dydt[j] = 0;

    }

  }

}


void slicing_mat(float dydt[NeuronPopulation*6], float dydt_tri[2], float y[NeuronPopulation*6], float y_tri[2], int i){

  for (int a = 0; a < 2; a++) {
    dydt_tri[a] = dydt[a+(2*i)];
    y_tri[a] = y[a+(2*i)];

  }
  
}


void update_trimat(float dydt[NeuronPopulation*6], float dydt_tri[2], float y[NeuronPopulation*6], float y_tri[2], int i){

  for (int a = 0; a < 2; a++) {
    dydt[a+(2*i)] = dydt_tri[a];
    y[a+(2*i)]= y_tri[a] ;

  }
  
}



//////////////////// END OF TRIANGLE GENERATOR //////////////////



/////////////////// NEURON PART /////////////////////////////


void matrix_mult_neu_2d(float array_1[8], float array[4][8], float array_final[4]){

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

  matrix_mult_neu_2d(I_Vec,TMAT,dydt);
}



void slicing_mat_neur(float dydt[NeuronPopulation*6], float dydt_neur[4], float y[NeuronPopulation*6], float y_neur[4], int i){

  for (int a = 0; a < 4; a++) {
    dydt_neur[a] = dydt[a+(2*NeuronPopulation)+(4*i)];
    y_neur[a] = y[a+(2*NeuronPopulation)+(4*i)];

  }
  
}


void update_neurmat(float dydt[NeuronPopulation*6], float dydt_neur[4], float y[NeuronPopulation*6], float y_neur[4], int i){

  for (int a = 0; a < 4; a++) {
    dydt[a+(2*NeuronPopulation)+(4*i)] = dydt_neur[a];
    y[a+(2*NeuronPopulation)+(4*i)]= y_neur[a] ;

  }
  
}



///////////////////////// END OF NEURON PART /////////////////////////////////




void compute_synapse_mat(float weights[NeuronPopulation][NeuronPopulation*2], float Synapse_Mat_fixed_term, float synapse_mat[NeuronPopulation][NeuronPopulation*2]){

  for (int i = 0; i < NeuronPopulation; i++) {
    for (int j = 0; j < (NeuronPopulation*2); j++){

      if(weights[i][j]<2.5){
        synapse_mat[i][j] = Synapse_Mat_fixed_term * pow(2,(((-Kp*weights[i][j])/(UT*log(2)))+9.9658));


      }

      else{
        synapse_mat[i][j] = 0;

      }

    }

  }
  
}


void intialize_y(float y[NeuronPopulation*6], float Vmem[NeuronPopulation]){

  int temp=0;
  for(int i=0; i<(NeuronPopulation*6);i++){

    if(i<NeuronPopulation*2){
      y[i] = VDD;
    }

    else if (i == (NeuronPopulation*2 +temp))
    {
      y[i] = -Offset_a;
      Vmem[temp/4]= -Offset_a;
      temp = temp+4;
    }

    else if (i == (NeuronPopulation*2 +temp-3))
    {
      y[i] = -Offset_b;
      
    }

    else if (i == (NeuronPopulation*2 +temp-2))
    {
      y[i] = -Offset_c;
      
    }

    else{

      y[i] = -Offset_c;

    }

  }
  
}

void update_Vmem(float y[NeuronPopulation*6], float Vmem[NeuronPopulation]){

  int temp=0;
  for(int i=0; i<(NeuronPopulation*6);i++){

    if (i == (NeuronPopulation*2 +temp))
    {
      Vmem[temp/4]= y[i];
      temp = temp+4;
    }

  }
  
}

void intialize_top_bot_idx(int top_idx[2*NeuronPopulation], int bot_idx[2*NeuronPopulation]){

  for(int i=0; i<(NeuronPopulation*2);i++){

    if(i%2==0){
      top_idx[i]=1;
      bot_idx[i]=0;

    }

    else{
      top_idx[i]=0;
      bot_idx[i]=1;
    }

  }
  
}


void InputCurrent_compute1(float array_1[2*NeuronPopulation], float array[NeuronPopulation][NeuronPopulation*2], float array_final[NeuronPopulation], float Inp_current[NeuronPopulation]){

  for (int i = 0; i < NeuronPopulation; i++) {
    for (int j = 0; j < 2*NeuronPopulation; j++){
      array_final[i] = array_final[i]+(array[i][j]*array_1[j]);

    }
    array_final[i]+=Inp_current[i];
  }
  
}


void InputCurrent_compute(float Synapse_Mat[NeuronPopulation][NeuronPopulation*2], float Scale_Mat[NeuronPopulation][NeuronPopulation*2], float EffCoupling, float Triangle_scaled[2*NeuronPopulation], int top_idx[2*NeuronPopulation], int bot_idx[2*NeuronPopulation], float Inp_current[NeuronPopulation], float InputCurrent[NeuronPopulation]){

  float temp[NeuronPopulation][NeuronPopulation*2] = {}, temp1[NeuronPopulation*2] = {};
  for (int i = 0; i < NeuronPopulation; i++) {
    for (int j = 0; j < 2*NeuronPopulation; j++){
      temp[i][j] = Synapse_Mat[i][j]*Scale_Mat[i][j];
    }
  }

  for(int i=0; i< 2*NeuronPopulation; i++){
    temp1[i] = pow(2, ( (-EffCoupling*Triangle_scaled[i]) + (86*top_idx[i]) - (2*bot_idx[i]) ) );
  }

  InputCurrent_compute1(temp1, temp, InputCurrent, Inp_current);

}


int Network_ODE(float dydt[NeuronPopulation*6], float y[NeuronPopulation*6], float Vmem[NeuronPopulation], float Synapse_Mat[NeuronPopulation][NeuronPopulation*2], float Inp_current[NeuronPopulation]){

  float CG_Coupling = 0.8;               //Control Gate Coupling Factor
  float EffCoupling = CG_Coupling*Kp;    //Combination of both the CG-gate and gate-channel coupling
  float Tri_Top_Off = 129.8426;          // Resting offset of the control gate of the excitatory synapse
  float Tri_Bot_Off = 14.427;            // Resting offset of the control gate of the inhibitory synapse
  float TriScaling  = 0.05;              // Scaling of the triangle generator outputs

  //COMPUTE INPUT CURRENT
  int top_idx[2*NeuronPopulation]={}, bot_idx[2*NeuronPopulation] ={};
  intialize_top_bot_idx(top_idx,bot_idx);

  float Triangle_scaled[2*NeuronPopulation]={};

  for (int i = 0; i < 2*NeuronPopulation; i++)
  {
    Triangle_scaled[i] = y[i]*TriScaling + top_idx[i]*Tri_Top_Off + bot_idx[i]*Tri_Bot_Off;
  }
  
  float Scale_Mat[NeuronPopulation][NeuronPopulation*2]={};

  for (int i = 0; i < NeuronPopulation; i++)
  {
    for (int j = 0; j < 2*NeuronPopulation; j++)
    {
      if(top_idx[j]==1){
        Scale_Mat[i][j] = pow(2,VDD-145.7947);
        
      }

      else{
        Scale_Mat[i][j] = -pow(2,Vmem[i]+13.2053);

      }
    }
  }
  
  float InputCurrent[NeuronPopulation] ={};

  InputCurrent_compute(Synapse_Mat, Scale_Mat, EffCoupling, Triangle_scaled, top_idx, bot_idx, Inp_current, InputCurrent);
   
  float T_MAT[2][2] = {{-71.2378,-13.4597},{-272.0602,-3.4893}};
  float dydt_tri[2]={}, y_tri[2]={}; 

  for(int i=0; i<NeuronPopulation; i++){
    slicing_mat(dydt,dydt_tri,y,y_tri,i);
    Tri_gen(dydt_tri,y_tri,T_MAT,Vmem[i]);
    update_trimat(dydt,dydt_tri,y,y_tri,i);
  }



///////// Neuron Part///////
 
 
  float EXPMEL = 7.3306;
  float EXPVSAT= 6.6116;
  float EXPMVGK = 0.2279;
  float EXPVNA = 0.0182;
  float EXPEK = 0.1364;

  float T_MAT_neu[4][8] = {{70.785082606739/1000, 20.2173858708657/1000,	1.32921811478727/1000,	2387.77324751470/1000,	-23893.5281204093/1000,	-9026.24819333978/1000,	135.171307831437/1000,	-232.108910489866/1000},
  {70.7850826067390/1000,	20.2173858708657/1000,	1.32921811478727/1000,	2387.77324751470/1000,	-23893.5281204093/1000,	-9026.24819333978/1000,	410.668585361385/1000,	-232.108910489866/1000},
  {69.6603580753146/1000,	19.8961460133584/1000,	1.30809778595160/1000,	4377.58428711028/1000,	-43804.8015540837/1000, -8882.82754034463/1000,	133.023532053490/1000,	425.533002564754/1000},
  {65.3065856956075/1000,	18.6526368875235/1000,	1.22634167432962/1000,	20587.9115563489/1000,	-206015.309127084/1000,	-8327.65081907309/1000,	124.709561300147/1000,	23636.4240515513/1000}};



  float dydt_neur[4]={}, y_neur[4]={}; 

  for(int i=0; i<NeuronPopulation; i++){
    slicing_mat_neur(dydt,dydt_neur,y,y_neur,i);
    FH_ODE(y_neur, InputCurrent[i], EXPVNA, EXPEK, EXPMEL, EXPVSAT, EXPMVGK, T_MAT_neu, dydt_neur);
    update_neurmat(dydt,dydt_neur,y,y_neur,i);
  }

}



int main() {

  // Define Network Connectivity Via the Synapse Weight Matrix
  //float weights[NeuronPopulation][NeuronPopulation*2]={2.5,2.5};
  float Syn = 0.26, Ih=0.22;
  float weights[NeuronPopulation][NeuronPopulation*2]={{2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, Ih},
         {2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, Ih},
         {2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, Ih},
         {Syn, 2.5, Syn, 2.5, Syn, 2.5, 2.5, 2.5}};

  //float ExternalCurrentIdx = 1, M_syn = 1.6667, Ithp = 1.25e-7, Kp = 0.75, VTp = 40.3955, UT= 0.025;
  //(Params.M_Syn*Ithp*2^(Kp*(VDD-VT0p)-VDD).*1e30 ./2^10) .*  2.^((-Kp*Weights/(UT*log(2)))+10) .* (+(Weights < 2.5));

  float Synapse_Mat_fixed_term = 2.1425, Synapse_Mat[NeuronPopulation][NeuronPopulation*2]={}; //2.1939

  compute_synapse_mat(weights,Synapse_Mat_fixed_term,Synapse_Mat);

  //Intial values for y, and Vmem array for all neurons
  float y[NeuronPopulation*6] ={}, Vmem[NeuronPopulation]={};
  intialize_y(y,Vmem);




  // Define simulation parameters
  float dt     = 0.001; // Units (ms)
  int   NSteps = 1.2e5;  // Total Number of steps to take in the simulation
  int PulseStartIdx = 50000;
  int PulseendIdx = 90000;
  float Inp_current[NeuronPopulation] = {};
  
  // For 1 neuron 20000-80000  - 0.06s/60ms
  // For 3 neuron 20000-120000 - 0.1s/100ms 


  float PulseHeight_1 = 50.0;
  float PulseHeight_2 = 60.0; 
  float PulseHeight_3 = 70.0;  
 

  // Open a file for writing
  FILE *fptr,*fptr1;
  fptr = fopen("WTA.txt","w");
  if(fptr == NULL){
      printf("File Error!");   
      exit(1);             
   }


  for(int i = 0; i < NSteps; i++){
    // Compute the current input current
    if(i >= PulseStartIdx && i <= PulseendIdx){
      Inp_current[0] = PulseHeight_1;
      Inp_current[1] = PulseHeight_2;
      Inp_current[2] = PulseHeight_3;
    }
    else{
      Inp_current[0] = 0;
      Inp_current[1] = 0;
      Inp_current[2] = 0;
    }

    // Run the actual ODE (Explicit Euler)
    float dydt[6*NeuronPopulation] = {};
    
    Network_ODE(dydt, y, Vmem, Synapse_Mat,Inp_current);
    

    // Take Euler Steps (y_next = y_current + dydt*dt)
    for(int j = 0; j < 6*NeuronPopulation; j++){
      y[j] += dydt[j]*dt;
    }

    update_Vmem(y,Vmem);

    // Print Vmem to file
    //Print_1dArray(fixed_point_dydt,4);


    //choose the y accordingly to print to the file.
    fprintf(fptr,"%f %f %f %f\n",y[8],y[12],y[16],y[20]);
    //fprintf(fptr,"%f %f %f %f\n",y[0],y[2],y[4],y[6]);



//    fprintf(fptr1,"%f\n",y[1]);
  }
  

  // Close file
  fclose(fptr);
//  fclose(fptr1);

  printf("success!!");

  return 0;

}