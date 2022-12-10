function plot_H(nx,ny,nz,Lx,Ly,Lz,PosTx,PosTy,PosTz)

[xcoor,ycoor,zcoor]=node_H(nx,ny,nz,Lx,Ly,Lz,PosTx,PosTy,PosTz);

load Hz.mat HZ
p_HZ=HZ(zcoor*nx*ny+ycoor*nx+1:zcoor*nx*ny+ycoor*nx+nx);

p_rHZ=real(p_HZ);
p_iHZ=imag(p_HZ);
hx=Lx/nx;hy=Ly/ny;hz=Lz/nz;

x=[hx/2:hx:Lx-hx/2]';

figure;
subplot(2,1,1)
plot(x,p_rHZ,'*')
legend('real Hz')
subplot(2,1,2)
plot(x,p_iHZ,'*')
legend('imag Hz')
return

v1=real(HZ);
v2=imag(HZ);


for iz=1:nz-1 
    q((iz-1)*nx+1:(iz-1)*nx+nx)=v1((ny/2-1)*nx+(iz-1)*ny*nx+1:(ny/2-1)*(nx)+(iz-1)*ny*(nx)+nx);
end
 
Qe = vec2mat(q,ny);
%save realHZ.mat Qe
r1=hx/2:hx:nx*hx-hx/2;
r2=hz:hz:(nz-1)*hz;
[t1,t2]=meshgrid(r1,r2);
figure
surf(t1,t2,Qe);
set(gca,'FontSize',30)
colorbar
legend('real HZ')

clear v1 Qe q r1 r2 t1 t2 h
 
for iz=1:nz-1 
    q((iz-1)*nx+1:(iz-1)*nx+nx)=v2((ny/2-1)*nx+(iz-1)*ny*nx+1:(ny/2-1)*(nx)+(iz-1)*ny*(nx)+nx);
end
 
Qe = vec2mat(q,ny);
%save imagHZ.mat Qe
r1=hx/2:hx:(nx)*hx-hx/2;
r2=hz:hz:(nz-1)*hz;
[t1,t2]=meshgrid(r1,r2);
figure
surf(t1,t2,Qe);
set(gca,'FontSize',30)
colorbar
legend('imag(HZtotal)')