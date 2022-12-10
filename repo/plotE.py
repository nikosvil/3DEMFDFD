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

file1=open('info.txt','r')
filec=file1.readlines()
print(filec)
file1.close()
sys.exit()
'''
nx=
ny=
nz=
Lx=
Ly=
Lz=
PosTx=
PosTy=
PosTz=
gl=
sigma0=
sigma=
sigmaAIR=

hx=Lx/nx
hy=Ly/ny
hz=Lz/nz
lenR=2.0
'''
coor='C8'
print('Transmitter Coordinates')
print('(',PosTx,PosTy,gl-PosTz,')')

PosRx=PosTx+lenR
PosRy=PosTy
PosRz=gl-PosTz

sPrim='0001'
sSec='001'
freq=3000


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


[zcoor,Ex_ycoor,Ey_ycoor]=near_node_E(nx,ny,nz,Lx,Ly,Lz,gl,PosRx,PosRy,PosRz)


EtotalX=dataEtotalX[zcoor*nx*(ny-1)+Ex_ycoor*nx+0:zcoor*nx*(ny-1)+Ex_ycoor*nx+nx]
EtotalY=dataEtotalY[zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)+0:zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)+nx-1]

EtotalXreal=dataEtotalX.real[zcoor*nx*(ny-1)+Ex_ycoor*nx+0:zcoor*nx*(ny-1)+Ex_ycoor*nx+nx]
ETXreal=dataETX.real[zcoor*nx*(ny-1)+Ex_ycoor*nx+0:zcoor*nx*(ny-1)+Ex_ycoor*nx+nx]
EtotalXimag=dataEtotalX.imag[zcoor*nx*(ny-1)+Ex_ycoor*nx+0:zcoor*nx*(ny-1)+Ex_ycoor*nx+nx]
ETXimag=dataETX.imag[zcoor*nx*(ny-1)+Ex_ycoor*nx+0:zcoor*nx*(ny-1)+Ex_ycoor*nx+nx]
#print(zcoor,zcoor*nx*(ny-1),Ex_ycoor,Ex_ycoor*nx,zcoor*nx*(ny-1)+Ex_ycoor*nx,zcoor*nx*(ny-1)+Ex_ycoor*nx+16)

EtotalYreal=dataEtotalY.real[zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)+0:zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)+nx-1]
ETYreal=dataETY.real[zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)+0:zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)+nx-1]
EtotalYimag=dataEtotalY.imag[zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)+0:zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)+nx-1]
ETYimag=dataETY.imag[zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)+0:zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)+nx-1]
xpp=int(nx/2)-1+zcoor*nx*(ny-1)+Ey_ycoor*nx 
xp=int(nx/2)-1

error(EtotalXreal,EtotalXimag,EtotalYreal,EtotalYimag,ETXreal,ETXimag,ETYreal,ETYimag)

E_Val_node(nx,ny,nz,hx,hy,hz,zcoor,Ex_ycoor,Ey_ycoor,EtotalXreal, EtotalYreal, ETXreal, ETYreal, EtotalXimag, EtotalYimag, ETXimag, ETYimag,sPrim, sSec, freq)

print('Transmitter position (',PosTx,PosTy,gl-PosTz,')')
#print('Receiver position (',PosRx,PosRy,PosRz,')')
print('****************************************')
print('Closest Ex Node at (',Ex_ycoor*hy+hy/2, Ex_ycoor*hy+hy,   (zcoor+1)*hz,')')
print('Closest Ey Node at (',Ey_ycoor*hy+hy,   Ey_ycoor*hy+hy/2, (zcoor+1)*hz,')')
print('****************************************')
print('Ex at \t (',Ex_ycoor*hy+7*hy/2, Ex_ycoor*hy+hy, (zcoor+1)*hz,')')
print('1D',ETXreal[xp+3])
print('3D',EtotalXreal[xp+3])
re=abs(ETXreal[xp+3]-EtotalXreal[xp+3])/abs(ETXreal[xp+3])
print('Deviation EX Real Central Node',abs(EtotalXreal[xp+3]-ETXreal[xp+3]),'%',100*re)
print('1D',ETXimag[xp+3])
print('3D',EtotalXimag[xp+3])
re=abs(ETXimag[xp+3]-EtotalXimag[xp+3])/abs(ETXimag[xp+3])
print('Deviation EX Imag Central Node',abs(EtotalXimag[xp+3]-ETXimag[xp+3]),'%',100*re)
print('*************************')
print('Ex at \t (',Ex_ycoor*hy+9*hy/2, Ex_ycoor*hy+hy,(zcoor+1)*hz,')')
print('1D',ETXreal[xp+4])
print('3D',EtotalXreal[xp+4])
re=abs(ETXreal[xp+4]-EtotalXreal[xp+4])/abs(ETXreal[xp+4])
print('Deviation',abs(EtotalXreal[xp+4]-ETXreal[xp+4]),'%',100*re)
print('1D',ETXimag[xp+4])
print('3D',EtotalXimag[xp+4])
re=abs(ETXimag[xp+4]-EtotalXimag[xp+4])/abs(ETXimag[xp+4])
print('Deviation',abs(EtotalXimag[xp+4]-ETXimag[xp+4]),'%',100*re)
print('*************************')
print('Ex at \t (',Ex_ycoor*hy-5*hy/2,Ex_ycoor*hy+hy,(zcoor+1)*hz,')')
print('1D',ETXreal[xp-3])
print('3D',EtotalXreal[xp-3])
re=abs(ETXreal[xp-3]-EtotalXreal[xp-3])/abs(ETXreal[xp-3])
print('Deviation',abs(EtotalXreal[xp-3]-ETXreal[xp-3]),'%',100*re)
print('1D',ETXimag[xp-3])
print('3D',EtotalXimag[xp-3])
re=abs(ETXimag[xp-3]-EtotalXimag[xp-3])/abs(ETXimag[xp-3])
print('Deviation',abs(EtotalXimag[xp-3]-ETXimag[xp-3]),'%',100*re)
print('*************************')
print('Ex at \t (',Ex_ycoor*hy-7*hy/2,Ex_ycoor*hy+hy,(zcoor+1)*hz,')')
print('1D',ETXreal[xp-4])
print('3D',EtotalXreal[xp-4])
re=abs(ETXreal[xp-4]-EtotalXreal[xp-4])/abs(ETXreal[xp-4])
print('Deviation',abs(EtotalXreal[xp-4]-ETXreal[xp-4]),100*re)
print('1D',ETXimag[xp-4])
print('3D',EtotalXimag[xp-4])
re=abs(ETXimag[xp-4]-EtotalXimag[xp-4])/abs(ETXimag[xp-4])
print('Deviation',abs(EtotalXimag[xp-4]-ETXimag[xp-4]),'%',100*re)
print('****************************************')

#sys.exit()
x = np.arange(hx/2,Lx,hx)
print('****************************************')
'''print('Node for Ey with respect to y and z axes (',Ey_ycoor*hy+hy, Ey_ycoor*hy+hy/2, (zcoor+1)*hz,')')
print('1D',ETYreal[xp],'3D',EtotalYreal[xp])
print('Deviation EY Real Central Node',abs(EtotalYreal[xp]-ETYreal[xp]))
print('1D',ETYimag[xp],'3D',EtotalYimag[xp])
print('Deviation EYImag Central Node',abs(EtotalYimag[xp]-ETYimag[xp]))
print('*************************')
print('Previous node for Ey at \t (',Ey_ycoor*hy,Ey_ycoor*hy+hy/2,(zcoor+1)*hz,')')
print('1D',ETYreal[xp-1],'3D',EtotalYreal[xp-1])
print('Deviation EY Real Previous Node',abs(EtotalYreal[xp-1]-ETYreal[xp-1]))
print('1D',ETYimag[xp-1],'3D',EtotalYimag[xp-1])
print('Deviation EY Imag Previous Node',abs(EtotalYimag[xp-1]-ETYimag[xp-1]))
print('*************************')
print('Next node for Ey at \t (',Ey_ycoor*hy+2*hy, Ey_ycoor*hy+hy/2,(zcoor+1)*hz,')')
print('1D',ETYreal[xp+1],'3D',EtotalYreal[xp+1])
print('Deviation EY Real Next Node',abs(EtotalYreal[xp+1]-ETYreal[xp+1]))
print('1D',ETYimag[xp+1],'3D',EtotalYimag[xp+1])
print('Deviation EY Imag Next Node',abs(EtotalYimag[xp+1]-ETYimag[xp+1]))'''
print('Ey at (',Ey_ycoor*hy+3*hy, Ey_ycoor*hy+hy/2, (zcoor+1)*hz,')')
print('1D',ETYreal[xp+2])
print('3D',EtotalYreal[xp+2])
re=abs(ETYreal[xp+2]-EtotalYreal[xp+2])/abs(ETYreal[xp+2])
print('Deviation',abs(EtotalYreal[xp+2]-ETYreal[xp+2]),'%',100*re)
print('1D',ETYimag[xp+2])
print('3D',EtotalYimag[xp+2])
re=abs(ETYimag[xp+2]-EtotalYimag[xp+2])/abs(ETYimag[xp+2])
print('Deviation',abs(EtotalYimag[xp+2]-ETYimag[xp+2]),'%',100*re)
print('*************************')
print('Ey at \t (',Ey_ycoor*hy+4*hy, Ey_ycoor*hy+hy/2,(zcoor+1)*hz,')')
print('1D',ETYreal[xp+3])
print('3D',EtotalYreal[xp+3])
re=abs(ETYreal[xp+3]-EtotalYreal[xp+3])/abs(ETYreal[xp+3])
print('Deviation',abs(EtotalYreal[xp+3]-ETYreal[xp+3]),'%',100*re)
print('1D',ETYimag[xp+3])
print('3D',EtotalYimag[xp+3])
re=abs(ETYimag[xp+3]-EtotalYimag[xp+3])/abs(ETYimag[xp+3])
print('Deviation',abs(EtotalYimag[xp+3]-ETYimag[xp+3]),'%',100*re)
print('*************************')
print('Ey at \t (',Ey_ycoor*hy-3*hy,Ey_ycoor*hy+hy/2,(zcoor+1)*hz,')')
print('1D',ETYreal[xp-4])
print('3D',EtotalYreal[xp-4])
re=abs(ETYreal[xp-4]-EtotalYreal[xp-4])/abs(ETYreal[xp-4])
print('Deviation EY',abs(EtotalYreal[xp-4]-ETYreal[xp-4]),'%',100*re )
print('1D',ETYimag[xp-4])
print('3D',EtotalYimag[xp-4])
re=abs(ETYimag[xp-4]-EtotalYimag[xp-4])/abs(ETYimag[xp-4])
print('Deviation EY',abs(EtotalYimag[xp-4]-ETYimag[xp-4]),'%',100*re)
print('*************************')
print('Ey at \t (',Ey_ycoor*hy-2*hy,Ey_ycoor*hy+hy/2,(zcoor+1)*hz,')')
print('1D',ETYreal[xp-3])
print('3D',EtotalYreal[xp-3])
re=abs(ETYreal[xp-3]-EtotalYreal[xp-3])/abs(ETYreal[xp-3])
print('Deviation',abs(EtotalYreal[xp-3]-ETYreal[xp-3]),'%',100*re)
print('1D',ETYimag[xp-3])
print('3D',EtotalYimag[xp-3])
re=abs(ETYimag[xp-3]-EtotalYimag[xp-3])/abs(ETYimag[xp-3])
print('Deviation',abs(EtotalYimag[xp-3]-ETYimag[xp-3]),'%',100*re)
print('*************************')
#sys.exit()
y = np.arange(hy,Ly,hy)
#print(type(x),x.shape,x)
#print(type(y),y.shape,y)
#print(type(EtotalXreal),EtotalXreal.shape)
#print(type(EtotalXimag),EtotalXimag.shape)
#print(type(ETXreal),ETXreal.shape)
#print(type(ETXimag),ETXimag.shape)
#print(type(EtotalYreal),EtotalYreal.shape)
#print(type(EtotalYimag),EtotalYimag.shape)
#print(type(ETYreal),ETYreal.shape)
#print(type(ETYimag),ETYimag.shape)
#sys.exit()
plt.figure(1)
plt.subplot(211)
plt.plot(x,EtotalXreal,'r*-',x,ETXreal,'b')
#plt.ticklabel_format(axis='x',style='sci',scilimits=(0,640))
#plt.ticklabel_format(axis='y',style='sci',scilimits=(-1,-11))
plt.legend(['EtotalX Real','ETX Real'])
plt.subplot(212)
plt.plot(x,EtotalXimag,'r*-',x,ETXimag,'b')
plt.legend(['EtotalX Imag','ETX Imag'])
#sys.exit()
plt.figure(2)
plt.subplot(211)
plt.plot(y,EtotalYreal,'r*-',y,ETYreal,'b')
#plt.ticklabel_format(axis='x',style='sci',scilimits=(0,640))
#plt.ticklabel_format(axis='y',style='sci',scilimits=(-10,-10))
plt.legend(['EtotalY Real','ETY Real'])
plt.subplot(212)
plt.plot(y,EtotalYimag,'r*-',y,ETYimag,'b')
plt.legend(['EtotalY Imag','ETY Imag'])
plt.show()
