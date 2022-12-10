import matplotlib 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as pltick
import math
import sys
import os
from near_node_E import near_node_E
from error import error
from E_Val_node import E_Val_node

n=64
Lx=640
Ly=640
Lz=640
nx=ny=nz=n
hx=Lx/nx
hy=Ly/ny
hz=Lz/nz
gl=208.0
lenR=2.0

PosTx=318.0
PosTy=318.0
PosTz=-42.0
coor='C11'
print('Creating Node Values File')

PosRx=PosTx+lenR
PosRy=PosTy
PosRz=gl-PosTz

sPrim='0001'
sSec='001'
freq=5000


'''******************************************'''
file1=open('EtotalX.txt','r')
EtotalX=file1.read()
file1.close()

EtotalX=EtotalX.replace('(','')
EtotalX=EtotalX.replace(')','j')
EtotalX=EtotalX.replace(',','+')
EtotalX=EtotalX.replace('+-','-')

file1=open('EtotalX.txt','w')
file1.write(EtotalX)
file1.close()

file1=open('EtotalY.txt','r')
EtotalY=file1.read()
file1.close()

EtotalY=EtotalY.replace('(','')
EtotalY=EtotalY.replace(')','j')
EtotalY=EtotalY.replace(',','+')
EtotalY=EtotalY.replace('+-','-')

file1=open('EtotalY.txt','w')
file1.write(EtotalY)
file1.close()

'''******************************************'''
file1=open('../../../EbackFiles/%s/E_Prim_%s_Sec_%s_f_%s/ETX%s.txt' %(coor,sPrim,sSec,freq,nx),'r')
ETX=file1.read()
file1.close()

ETX=ETX.replace('   -','-')
ETX=ETX.replace('    -','-')
ETX=ETX.replace('   +','+')
ETX=ETX.replace('    +','+')
ETX=ETX.replace('    ','+')

file1=open('../../../EbackFiles/%s/E_Prim_%s_Sec_%s_f_%s/ETX%s.txt' %(coor,sPrim,sSec,freq,nx),'w')
file1.write(ETX)
file1.close()

file1=open('../../../EbackFiles/%s/E_Prim_%s_Sec_%s_f_%s/ETY%s.txt' %(coor,sPrim,sSec,freq,nx),'r')
ETY=file1.read()
file1.close()

ETY=ETY.replace('   -','-')
ETY=ETY.replace('    -','-')
ETY=ETY.replace('   +','+')
ETY=ETY.replace('    +','+')
ETY=ETY.replace('    ','+')

file1=open('../../../EbackFiles/%s/E_Prim_%s_Sec_%s_f_%s/ETY%s.txt' %(coor,sPrim,sSec,freq,nx),'w')
file1.write(ETY)
file1.close()
'''******************************************'''
dataEtotalX=np.loadtxt('EtotalX.txt' ,dtype=complex)
dataEtotalY=np.loadtxt('EtotalY.txt' ,dtype=complex)
dataETX=np.loadtxt('../../../EbackFiles/%s/E_Prim_%s_Sec_%s_f_%s/ETX%s.txt' %(coor,sPrim,sSec,freq,nx),dtype=complex)
dataETY=np.loadtxt('../../../EbackFiles/%s/E_Prim_%s_Sec_%s_f_%s/ETY%s.txt' %(coor,sPrim,sSec,freq,nx),dtype=complex)
#print(type(datac),datac.shape)
#print(type(dataa),dataa.shape)

dataEtotalXreal=dataEtotalX.real
dataEtotalYreal=dataEtotalY.real
dataETXreal=dataETX.real
dataETYreal=dataETY.real

i=0
file1=open('All_Nodes_EtotalX_Values_%s_%s_%s.txt' %(sPrim,sSec,freq), 'w')
file3=open('Error_Below_6_EX.txt', 'w')
for zc in range (10, int(Lz), int(hz)):
    for yc in range (10, int(Ly), int(hy)):
        for xc in range (5, int(Lx), int(hx)):
            dfs=round(math.sqrt((PosTx-xc)**2+(PosTy-yc)**2+(gl-PosTz-zc)**2))
            abs_e=abs(dataEtotalXreal[i]-dataETXreal[i])
            rel_e=100*(abs_e/abs(dataETXreal[i]))
            if rel_e<6.0:
             file3.write('%.2f %d %d %d %s\n' % (rel_e,xc,yc,zc,dfs))
            file1.write('%s %s %s %.5e %.2f %d %d %d %s\n' % (i+1,dataEtotalXreal[i],dataETXreal[i],abs_e,rel_e,xc,yc,zc,dfs))
            i=i+1
file1.close()
file3.close()

i=0
file2=open('All_Nodes_EtotalY_Values_%s_%s_%s.txt' %(sPrim,sSec,freq), 'w')
file3=open('Error_Below_6_EY.txt', 'w')
for zc in range (10, int(Lz), int(hz)):
    for yc in range (5, int(Ly), int(hy)):
        for xc in range (10, int(Lx), int(hx)):
            dfs=round(math.sqrt((PosTx-xc)**2+(PosTy-yc)**2+(gl-PosTz-zc)**2))
            abs_e=abs(dataEtotalYreal[i]-dataETYreal[i])
            rel_e=100*(abs_e/abs(dataETYreal[i]))
            if rel_e<6.0:
             file3.write('%.2f %d %d %d %s\n' % (rel_e,xc,yc,zc,dfs))
            file2.write('%s %s %s %.5e %.2f %d %d %d %s\n' % (i+1,dataEtotalYreal[i],dataETYreal[i],abs_e,rel_e,xc,yc,zc,dfs))
            i=i+1
file2.close()
file3.close()

