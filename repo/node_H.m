function [xcoor,ycoor, zcoor]=node_H(nx,ny,nz,Lx,Ly,Lz,PosTx,PosTy,PosTz)


hx=Lx/nx;hy=Ly/ny;hz=Lz/nz;
PosRx=PosTx;PosRy=PosTy;PosRz=PosTz;

zcoor=round(PosRz/hz)-1;
ycoor=round(PosRy/hy);
xcoor=round(PosRx/hx);

if (abs((ycoor-1)*hy+hy-PosRy)<=abs(ycoor*hy+hy-PosRy));
  ycoor=ycoor-1;
end  
sel=input('Select Different Height?\n','s');
if (sel=='y'|| sel=='Y')
  chnode_z=input('Move towards the z-axis\n');
  zcoor=zcoor+chnode_z;
end
sel=input('Move towards the y-axis?\n','s');
if (sel=='y'|| sel=='Y')
  chnode_y=input('How many nodes?\n');
  ycoor=ycoor+chnode_y;
end

end