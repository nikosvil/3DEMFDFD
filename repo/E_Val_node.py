import math
import numpy as np

def E_Val_node(nx,ny,nz,hx,hy,hz,zcoor,Ex_ycoor,Ey_ycoor,EtotalXreal, EtotalYreal, ETXreal, ETYreal, EtotalXimag, EtotalYimag, ETXimag, ETYimag,sPrim, sSec, freq):

    EX_Index_First_in_Line=zcoor*nx*(ny-1)+Ex_ycoor*nx
    EX_Index_Last_in_Line=zcoor*nx*(ny-1)+Ex_ycoor*nx+nx

    EY_Index_First_in_Line=zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)
    EY_Index_Last_in_Line=zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)+nx-1

    xc=hx/2
    yc=Ex_ycoor*hy+hy
    zc=(zcoor+1)*hz

    file1=open('EtotalX_Node_Values_%s_%s_%s.txt' %(sPrim,sSec,freq), 'w')
    for i in range (0,len(EtotalXreal)):
        abs_e=abs(EtotalXreal[i]-ETXreal[i])
        rel_e=100*(abs_e/abs(ETXreal[i]))
        file1.write('%s %s %s %s %.2f %s %s %s\n' % (i+1,EtotalXreal[i],ETXreal[i],abs_e,rel_e,xc,yc,zc))
        xc=xc+hx
                
    file1.close()

    xc=hx
    yc=Ey_ycoor*hy+hy/2
    zc=(zcoor+1)*hz

    file2=open('EtotalY_Node_Values_%s_%s_%s.txt' %(sPrim,sSec,freq), 'w')
    for i in range (0,len(EtotalYreal)):
        abs_e=abs(EtotalYreal[i]-ETYreal[i])
        rel_e=100*(abs_e/abs(ETYreal[i]))
        file2.write('%s %s %s %s %.2f %s %s %s\n' % (i+1,EtotalYreal[i],ETYreal[i],abs_e,rel_e,xc,yc,zc))
        xc=xc+hx
                
    file2.close()
    return


