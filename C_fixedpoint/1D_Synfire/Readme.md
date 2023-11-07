# Running 1D_Synfire
Compile the C file using: gcc Synfire_1D_3neur_fixed.c -o Syn -lm
Run the C file using: ./Syn   ->  This will generate Synfire.txt (membrane voltages of each neuron)
Run the Matlab plotting script (Plot_Synfire_1D_3neur.m) to read-in and plot the output.

# Running 2D_Synfire
Compile the C file using: gcc Synfire_1D_100kneur_fixed.c -o Syn2 -lm
Run the C file using: ./Syn2   ->  This will generate Synfire_100K.txt (spike times/raster... this is done because saving to file is extremely expensive). This code might take on the order of days to run.
Run the Matlab plotting script (Synfire_100k.m) to read-in and plot the output.