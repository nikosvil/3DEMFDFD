# 3DEMFDFD

GAIA 3DEMFDFD is a 4th-order finite difference solver used for the approximation of the residual electric field of a three-dimensional half-space. The repository contains 2 directories with the necessary files for the solver to provide results, an empty temp folder and the test folder as well as a config file and the shell script. The directories are named as

1. halfspace
2. repo
3. temp
4. test

The first directory contains the Fortran files and the Makefile as well as the subdirectory Results where the approximation of each E component is saved on .txt and .mat files at the end of the computation. The second directory contains key files that compute the primary electric field intensity in Matlab as well as scripts that convert the output E files to Matlab variables. The temp directory is being used in as a temporary files folder. The test folder contains the results of the program's execution (.txt and .mat files) when using the specific parameters of the uploaded config.dat file - to be used for testing purposes.

**Instructions**

The config file contains the necessary parameters of the physical problem, line after line, in the following manner:

Length-of-x-side-in-the-domain Length-of-y-side-in-the-domain Length-of-z-side-in-the-domain  
Position-of-transmitter-in-x-direction Position-of-transmitter-in-y-direction Position-of-transmitter-in-z-direction*  
Ground-in-z-direction  
Background-resistivity  
Foreground-resistivity  
Frequency

Example:

640 640 640
320 320 -22
208
0.001
0.01
3000

Distance is given in m, resistivity in Ohm and frequency in Hz.

***Important Notices**:

The Position-of-transmitter-in-z-direction in the config file should be entered as a negative number which refers to the height of the transmitter above the ground. For instance, (the Position-of-transmitter-in-z-direction=) -22 value when the ground (Ground-in-z-direction) is set to 208 indicates the transmitter to be at 230.

The solver's parameters (tolerance and maximum number of steps for the iterative method) may be altered from within the main.f file located in halfspace folder.

The gaia.sh script launches Matlab without the desktop and runs the necessary .m files using the matlab command from the shell. The user may need to modify the path in lines 106 and 123 of the gaia.sh script depending on the installation path and the Matlab version available.

The Fortran compiler used is pgf90 (PGI compilers are free to download and use). Python3 should also be installed in the system for the post-processing conversion.
 
To use the solver, provide the necessary parameters in the config.dat file as explained above and then run the command **gaia.sh < config.dat**

At the end of the program's execution, the results are being written in .txt files (EtotalX.txt, EtotalY.txt, EX.txt, EY.txt, EZ.txt) and .mat files (ETX.mat, ETY.mat, Ex.mat, Ey.mat, TotalEX.mat, TotalEY.mat in the f2m folder) for further process. EX, EY and EZ (Ex, Ey respectively) refer to the values of the secondary electric field components, EtotalX and EtotalY (TotalEX and TotalEY respectively) refer to the total electric field (background + secondary) and ETX, ETY are the analytical solution values of the components. File info.txt contains information on the problem's parameters and total execution time.

***Parameters of the implementation**:

The main.f contains the parameters of the solver (maximum number of steps of the iterative method, tolerance) and fixed parameters like magnetic permeability and air resistivity which may be modified if desired.
