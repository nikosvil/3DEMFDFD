# ===================================================================
# Title: convert_files.py
# Authors: N. D. Vilanakis
# Details: Applied Mathematics and Computers Lab, Technical University of Crete
# * A 3D frequency-domain electromagnetic solver employing a high order compact finite-difference scheme
# * N. D. Vilanakis, N. Economou, E. N. Mathioudakis, A. Vafidis
# * Computers and Geosciences
#====================================================================
# Script that reads .txt files that contain the total electric field 
# component values produced by the Fortran program and writes new ones
# complex-like strings to be processed by Matlab
#====================================================================
import math
import numpy

'''******************************************'''
file1=open('../temp/EtotalX.txt','r')
EtotalX=file1.read()
file1.close()

EtotalX=EtotalX.replace('(','')
EtotalX=EtotalX.replace(')','j;')
EtotalX=EtotalX.replace(',','+')
EtotalX=EtotalX.replace('+-','-')

file1=open('../temp/mEtotalX.txt','w')
file1.write(EtotalX)
file1.close()

file1=open('../temp/EtotalY.txt','r')
EtotalY=file1.read()
file1.close()

EtotalY=EtotalY.replace('(','')
EtotalY=EtotalY.replace(')','j;')
EtotalY=EtotalY.replace(',','+')
EtotalY=EtotalY.replace('+-','-')

file1=open('../temp/mEtotalY.txt','w')
file1.write(EtotalY)
file1.close()

file1=open('../temp/EX.txt','r')
EX=file1.read()
file1.close()

EX=EX.replace('(','')
EX=EX.replace(')','j;')
EX=EX.replace(',','+')
EX=EX.replace('+-','-')

file1=open('../temp/mEX.txt','w')
file1.write(EX)
file1.close()

file1=open('../temp/EY.txt','r')
EY=file1.read()
file1.close()

EY=EY.replace('(','')
EY=EY.replace(')','j;')
EY=EY.replace(',','+')
EY=EY.replace('+-','-')

file1=open('../temp/mEY.txt','w')
file1.write(EY)
file1.close()

file1.close()
