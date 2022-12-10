import math
import sys
import os
from scipy.io import savemat
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

file1=open('../temp/HZ.txt','r')
HZ=file1.read()
file1.close()

HZ=HZ.replace('(','')
HZ=HZ.replace(')','j;')
HZ=HZ.replace(',','+')
HZ=HZ.replace('+-','-')

file1=open('../temp/mHZ.txt','w')
file1.write(HZ)
file1.close()








