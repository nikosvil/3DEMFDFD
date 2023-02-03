# 3DEMFDFD

GAIA 3DEMFDFD is a 4th-order finite difference solver used for the approximation of the residual electric field of a three-dimensional half-space. The repository contains 2 directories with the necessary files for the solver to provide results and an empty temp folder as well as a config file and the shell script. The directories are named as

1. halfspace
2. repo
3. temp

The first directory contains the Fortran files and the Makefile as well as the subdirectory Results where the approximation of each E component exist. The second directory contains the file that computes the primary electric field intensity in Matlab. The last directory is being used as a temporary files folder.

**Instructions**

The config file contains the necessary parameters of the physical problem in the following manner:

Length-of-x-side-in-the-domain Length-of-y-side-in-the-domain Length-of-z-side-in-the-domain  
Position-of-transmitter-in-x-direction Position-of-transmitter-in-y-direction Position-of-transmitter-in-z-direction*  
Ground-in-z-direction  
Background-resistivity  
Foreground-resistivity  
Frequency  

Distance is given in m, resistivity in Ohm and frequency in Hz.

***Important Notice**: The z component is assumed to be zero at the ground surface and positive towards the Earth. Thus, z means depth, which is positive below the ground surface and negative above it. The Position-of-transmitter-in-z-direction in the config file should be entered as a negative number when referring to the height of the transmitter above the ground surface. For instance, (the Position-of-transmitter-in-z-direction)= -22 value when the ground surface (Ground-in-z-direction) is set to 208 indicates that the transmitter is at 230.

The solver's parameters (tolerance and maximum number of steps for the iterative method) may be altered from within the main.f file located in halfspace folder.

To use the solver, provide the necessary parameters in the config.dat file as explained above and then run the command gaia.sh < config.dat**
