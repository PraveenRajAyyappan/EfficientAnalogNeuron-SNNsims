# Simulating a Single Neuron
Compile the C file using: gcc neuron_fixed.c -o neu -lm
Run the C file using: ./neu   ->  This will generate Output.txt (membrane voltages of each neuron)
Run the Matlab plotting script (Plot_Neuron.m) to read-in and plot the output (membrane voltage).

# Generating IF curve
Compile the C file using: gcc IF_curve.c -o IF -lm
Run the C file using: ./IF   ->  This will generate multiple Vmem_Iin_%d.txt files (membrane voltages of each neuron) where %d represents the range of IF curve input current.
Run the Matlab plotting script (Plot_IFcurve.m) to read-in and plot the fixed point IF curve.
