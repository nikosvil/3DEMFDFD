import math
import numpy as np

def near_node_E(nx,ny,nz,Lx,Ly,Lz,gl,PosRx,PosRy,PosRz):
 PosTx=PosRx-2
 PosTy=PosRy
 PosTz=PosRz
 hx=Lx/nx
 hy=Ly/ny
 hz=Lz/nz

#x-z plane
 zcoor=round(PosRz/hz)-1
 Ex_ycoor=round(PosRy/hy)
 Ey_ycoor=round(PosRy/hy)

 if (abs((Ey_ycoor-1)*hy+hy/2-PosRy)<=abs(Ey_ycoor*hy+hy/2-PosRy)):
  Ey_ycoor=Ey_ycoor-1
 if (abs((Ex_ycoor-1)*hy+hy-PosRy)<=abs(Ex_ycoor*hy+hy-PosRy)):
  Ex_ycoor=Ex_ycoor-1

 sel=str(input('Select Different Height?\n'))
 if (sel=='y'or sel=='Y'):
  chnode_z=int(input('Move towards the z-axis\n')) 
  zcoor=zcoor+chnode_z

 sel=str(input('Select Different Node for Ex?\n'))
 if (sel=='y'or sel=='Y'):
  chnodeEx_y=int(input('Move Ex towards the y-axis\n'))
  Ex_ycoor=Ex_ycoor+chnodeEx_y
 
 sel=str(input('Select Different Node for Ey?\n'))
 if (sel=='y'or sel=='Y'):
  chnodeEy_y=int(input('Move Ex towards the y-axis\n'))
  Ey_ycoor=Ey_ycoor+chnodeEy_y

 return zcoor,Ex_ycoor,Ey_ycoor


