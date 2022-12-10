function [Ex_xcoor,Ey_xcoor,Ex_ycoor, Ey_ycoor,zcoor]=node_E(hx,hy,hz,PosTx,PosTy,PosTz)

PosRx=PosTx;PosRy=PosTy;PosRz=PosTz;


Ex_xcoor=round(PosRx/hx);Ey_xcoor=round(PosRx/hx);
Ex_ycoor=round(PosRy/hy);Ey_ycoor=round(PosRy/hy);

zcoor=round(PosRz/hz)-1;

if (abs((Ey_ycoor-1)*hy+hy/2-PosRy)<=abs(Ey_ycoor*hy+hy/2-PosRy));
  Ey_ycoor=Ey_ycoor-1;
end  
if (abs((Ex_ycoor-1)*hy+hy-PosRy)<=abs(Ex_ycoor*hy+hy-PosRy));
  Ex_ycoor=Ex_ycoor-1;
end
fprintf('Transmitter placed at %d %d %d\n',PosTx,PosTy,PosTz);
fprintf('Distance between nodes on y-axis:%d, on z-axis:%d\n', hy,hz);
fprintf('Displaying 4 nodes on x-axis on each side of transmitter, starting 20m away\n');
sel=input('Select Different Height?\n','s');
if (sel=='y'|| sel=='Y')
  chnode_z=input('Move towards the z-axis\n');
  zcoor=zcoor+chnode_z;
end

sel=input('Select Different Node for Ex?\n','s');
if (sel=='y'|| sel=='Y')
 chnodeEx_y=input('Move Ex towards the y-axis\n');
 Ex_ycoor=Ex_ycoor+chnodeEx_y;
end
sel=input('Select Different Node for Ey?\n','s');
if (sel=='y'|| sel=='Y')
 chnodeEy_y=input('Move Ex towards the y-axis\n');
 Ey_ycoor=Ey_ycoor+chnodeEy_y;
end

end