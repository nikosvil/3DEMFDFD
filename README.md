# 3DEMFDFD

GAIA 3DEMFDFD is a 4th-order finite difference solver used for the approximation of the residual electric field of a three-dimensional half-space. The repository contains 2 directories with the necessary files for the solver to provide results and an empty temp folder as well as a config file and the shell script. The directories are named as

1. halfspace
2. repo
3. temp

The first directory contains the Fortran files and the Makefile as well as the subdirectory Results where the approximation of each E component and the Magnetic Field z-component at the end of the computation. The second directory contains key files that compute the primary electric field intensity in Matlab as well as scripts that convert the output E files to Matlab variables. The last directory is being used in as a temporary files folder.

**Instructions**

The config file contains the necessary parameters of the physical problem in the following manner:

Length-of-x-side-in-the-domain Length-of-y-side-in-the-domain Length-of-z-side-in-the-domain  
Position-of-transmitter-in-x-direction Position-of-transmitter-in-y-direction Position-of-transmitter-in-z-direction*  
Ground-in-z-direction  
Background-resistivity  
Foreground-resistivity  
Frequency  

Distance is given in m, resistivity in Ohm and frequency in Hz.

***Important Notices**:

The Position-of-transmitter-in-z-direction in the config file should be entered as a negative number which refers to the height of the transmitter above the ground. For instance, (the Position-of-transmitter-in-z-direction=) -22 value when the ground (Ground-in-z-direction) is set to 208 indicates the transmitter to be at 230.

The solver's parameters (tolerance and maximum number of steps for the iterative method) may be altered from within the main.f file located in halfspace folder.

The gaia.sh script launches Matlab without the desktop and runs the necessary .m files using the matlab command from the shell. The user may need to modify the path in lines 106 and 123 of the gaia.sh script depending on the installation path and the Matlab version available.

The Fortran compiler used is pgf90. Python3 should also be installed in the system for the post-processing conversion.
 
To use the solver, provide the necessary parameters in the config.dat file as explained above and then run the command gaia.sh < config.dat**
