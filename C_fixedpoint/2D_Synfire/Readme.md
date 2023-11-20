# Running 2D_Synfire
Compile the C file using: gcc Synfire_2D_100Kneur_fixed.c -o Syn -lm

Run the C file using: ./Syn   ->  This will generate Synfire_num.txt (raster of each neuron). This code might take upto a day to run.

Run the Matlab plotting script (MakeRasterPlot.m) to read-in and plot the output.

The main parameters that one can tune are: 

INPUTNEURONS -> The number of neurons in the first layer to which a stimulus will be applied during a brief excitatory period 

INPUTNOISE -> Standard deviation of the voltages on the control gates of the FG synapses (in V)

SimTime -> How many timesteps to run the simulation (max value will be proportional to the number of layers, but typically the spike packet dies much faster)

Syn -> Mean synaptic weights

NeuronPopulation -> Total number of neurons in the network

Group -> Number of neurons per layer of the network

The idea is that spike trains can only survive and propagate in such a heterogeneous 2D synfire chain if the standard deviation of the synaptic weights (i.e., INPUTNOISE) is low enough given the number of neurons in the input pool to which the stimulus is being applied (INPUTNEURONS), or if INPUTNEURONS is high given INPUTNOISE. For more information, please see the corresponding journal paper.