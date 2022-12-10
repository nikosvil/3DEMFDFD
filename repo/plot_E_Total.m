function plot_E_Total(nx,ny,nz,hx,hy,hz,rETx,iETx,rETy,iETy,rEx_Total,iEx_Total,rEy_Total,iEy_Total)
 
%##################################### 
v1=rEy_Total;
v2=iEy_Total;
 
for iz=1:nz-1 
    q((iz-1)*(nx-1)+1:(iz-1)*(nx-1)+nx-1)=v1((ny/2-1)*(nx-1)+(iz-1)*ny*(nx-1)+1:(ny/2-1)*(nx-1)+(iz-1)*ny*(nx-1)+nx-1);
end
Qe = vec2mat(q,ny-1);
r1=hx:hx:(nx-1)*hx;
r2=hz:hz:(nz-1)*hz;
[t1,t2]=meshgrid(r1,r2);
h=figure;
surf(t1,t2,Qe);
set(gca,'FontSize',30)
colorbar
legend('real(EtotalY)')
% savefig(h,(sprintf('%sx%sx%s_real(Etotal_y).fig',int2str(nx),int2str(ny),int2str(nz))));
%saveas(h,(sprintf('%sx%sx%s_real(Ey-Back).png',int2str(nx),int2str(ny),int2str(nz))));
clear v1 Qe q r1 r2 t1 t2 h
 
for iz=1:nz-1 
    q((iz-1)*(nx-1)+1:(iz-1)*(nx-1)+nx-1)=v2((ny/2-1)*(nx-1)+(iz-1)*ny*(nx-1)+1:(ny/2-1)*(nx-1)+(iz-1)*ny*(nx-1)+nx-1);
end
Qe = vec2mat(q,ny-1);
r1=hx:hx:(nx-1)*hx;
r2=hz:hz:(nz-1)*hz;
[t1,t2]=meshgrid(r1,r2);
h=figure;
surf(t1,t2,Qe);
set(gca,'FontSize',30)
colorbar
legend('imag(EtotalY)')

% savefig(h,(sprintf('%sx%sx%s_imag(Etotal_y).fig',int2str(nx),int2str(ny),int2str(nz))));
%saveas(h,(sprintf('%sx%sx%s_imag(Ey-Back).png',int2str(nx),int2str(ny),int2str(nz))));
clear v1 v2 Qe q r1 r2 t1 t2 h
%#####################################     
v1=rEx_Total;
v2=iEx_Total;
 
for iz=1:nz-1
    for k=1:ny-1
        q((iz-1)*(ny-1)+k)=v1((k-1)*(nx)+(iz-1)*nx*(ny-1)+(nx/2));
    end
end
Qe = vec2mat(q,ny-1);
r1=hy:hy:(ny-1)*hy;
r2=hz:hz:(nz-1)*hz;
[t1,t2]=meshgrid(r1,r2);
h=figure;
surf(t1,t2,Qe);
set(gca,'FontSize',30)
colorbar
legend('real(EtotalX)')
% savefig(h,(sprintf('%sx%sx%s_real(Etotal_x).fig',int2str(nx),int2str(ny),int2str(nz))));
%saveas(h,(sprintf('%sx%sx%s_real(Ex-Back).png',int2str(nx),int2str(ny),int2str(nz))));
clear v1 Qe q r1 r2 t1 t2 h
   
for iz=1:nz-1
    for k=1:ny-1
        q((iz-1)*(ny-1)+k)=v2((k-1)*(nx)+(iz-1)*nx*(ny-1)+(nx/2));
    end
end
Qe = vec2mat(q,ny-1);
r1=hy:hy:(ny-1)*hy;
r2=hz:hz:(nz-1)*hz;
[t1,t2]=meshgrid(r1,r2);
h=figure;
surf(t1,t2,Qe);
set(gca,'FontSize',30)
colorbar
legend('imag(EtotalX)')
% savefig(h,(sprintf('%sx%sx%s_imag(Etotal_x).fig',int2str(nx),int2str(ny),int2str(nz))));
% saveas(h,(sprintf('%sx%sx%s_imag(Ex-Back).png',int2str(nx),int2str(ny),int2str(nz))));
clear v1 v2 Qe q r1 r2 t1 t2 h
%#####################################     
    
v1=rETy;
v2=iETy;
 
for iz=1:nz-1 
     q((iz-1)*(nx-1)+1:(iz-1)*(nx-1)+nx-1)=v1((ny/2-1)*(nx-1)+(iz-1)*ny*(nx-1)+1:(ny/2-1)*(nx-1)+(iz-1)*ny*(nx-1)+nx-1);
end
Qe = vec2mat(q,ny-1);
r1=hx:hx:(nx-1)*hx;
r2=hz:hz:(nz-1)*hz;
[t1,t2]=meshgrid(r1,r2);
h=figure;
surf(t1,t2,Qe);
set(gca,'FontSize',30)
colorbar
legend('real(ETY)')
% savefig(h,(sprintf('%sx%sx%s_real(EYY=Ety-ETY).fig',int2str(nx),int2str(ny),int2str(nz))));
% % saveas(h,(sprintf('%sx%sx%s_real(EYY=Ety-ETY).png',int2str(nx),int2str(ny),int2str(nz))));
clear v1 Qe q r1 r2 t1 t2 h
 
for iz=1:nz-1 
    q((iz-1)*(nx-1)+1:(iz-1)*(nx-1)+nx-1)=v2((ny/2-1)*(nx-1)+(iz-1)*ny*(nx-1)+1:(ny/2-1)*(nx-1)+(iz-1)*ny*(nx-1)+nx-1);
end
Qe = vec2mat(q,ny-1);
r1=hx:hx:(nx-1)*hx;
r2=hz:hz:(nz-1)*hz;
[t1,t2]=meshgrid(r1,r2);
h=figure;
surf(t1,t2,Qe);
set(gca,'FontSize',30)
colorbar
legend('imag(ETY)')

% savefig(h,(sprintf('%sx%sx%s_imag(EYY=Ety-ETY).fig',int2str(nx),int2str(ny),int2str(nz))));
% % % saveas(h,(sprintf('%sx%sx%s_imag(EYY=Ety-ETY).png',int2str(nx),int2str(ny),int2str(nz))));
% % clear v1 v2 Qe q r1 r2 t1 t2 h
%    
%     
v1=rETx;
v2=iETx;
 
for iz=1:nz-1
    for k=1:ny-1
        q((iz-1)*(ny-1)+k)=v1((k-1)*(nx)+(iz-1)*nx*(ny-1)+(nx/2));
    end
end
Qe = vec2mat(q,ny-1);
r1=hy:hy:(ny-1)*hy;
r2=hz:hz:(nz-1)*hz;
[t1,t2]=meshgrid(r1,r2);
h=figure;
surf(t1,t2,Qe);
set(gca,'FontSize',30)
colorbar
legend('real(ETX)')
% savefig(h,(sprintf('%sx%sx%s_real(EXX=Etx-ETX).fig',int2str(nx),int2str(ny),int2str(nz))));
% % saveas(h,(sprintf('%sx%sx%s_real(EXX=Etx-ETX).png',int2str(nx),int2str(ny),int2str(nz))));
clear v1 Qe q r1 r2 t1 t2 h
 
for iz=1:nz-1
    for k=1:ny-1
        q((iz-1)*(ny-1)+k)=v2((k-1)*(nx)+(iz-1)*nx*(ny-1)+(nx/2));
    end
end
Qe = vec2mat(q,ny-1);
r1=hy:hy:(ny-1)*hy;
r2=hz:hz:(nz-1)*hz;
[t1,t2]=meshgrid(r1,r2);
h=figure;
surf(t1,t2,Qe);
set(gca,'FontSize',30)
colorbar
legend('imag(ETX)')
%savefig(h,(sprintf('%sx%sx%s_imag(EXX=Etx-ETX).fig',int2str(nx),int2str(ny),int2str(nz))));
% % saveas(h,(sprintf('%sx%sx%s_imag(EXX=Etx-ETX).png',int2str(nx),int2str(ny),int2str(nz))));
% clear v1 v2 Qe q r1 r2 t1 t2 h

return








v1=real(ESX);
v2=imag(ESX);
 
for iz=1:nz-1
    for k=1:ny-1
        q((iz-1)*(ny-1)+k)=v1((k-1)*(nx)+(iz-1)*nx*(ny-1)+(nx/2));
    end
end
Qe = vec2mat(q,ny-1);
r1=hy:hy:(ny-1)*hy;
r2=hz:hz:(nz-1)*hz;
[t1,t2]=meshgrid(r1,r2);
h=figure;
surf(t1,t2,Qe);
set(gca,'FontSize',30)
colorbar
legend('real(ESX)')

% savefig(h,(sprintf('%sx%sx%s_real(ESX).fig',int2str(nx),int2str(ny),int2str(nz))));
% saveas(h,(sprintf('%sx%sx%s_real(Etx=Ex+EX).png',int2str(nx),int2str(ny),int2str(nz))));
clear v1 Qe q r1 r2 t1 t2 h
 
for iz=1:nz-1
    for k=1:ny-1
        q((iz-1)*(ny-1)+k)=v2((k-1)*(nx)+(iz-1)*nx*(ny-1)+(nx/2));
    end
end
Qe = vec2mat(q,ny-1);
r1=hy:hy:(ny-1)*hy;
r2=hz:hz:(nz-1)*hz;
[t1,t2]=meshgrid(r1,r2);
h=figure;
surf(t1,t2,Qe);
set(gca,'FontSize',30)
colorbar
legend('imag(ESX)')
% savefig(h,(sprintf('%sx%sx%s_imag(ESX).fig',int2str(nx),int2str(ny),int2str(nz))));
% saveas(h,(sprintf('%sx%sx%s_imag(Etx=Ex+EX).png',int2str(nx),int2str(ny),int2str(nz))));
clear v1 v2 Qe q r1 r2 t1 t2
    
    
v1=real(ESY);
v2=imag(ESY);
 
for iz=1:nz-1 
     q((iz-1)*(nx-1)+1:(iz-1)*(nx-1)+nx-1)=v1((ny/2-1)*(nx-1)+(iz-1)*ny*(nx-1)+1:(ny/2-1)*(nx-1)+(iz-1)*ny*(nx-1)+nx-1);
end
Qe = vec2mat(q,ny-1);
r1=hx:hx:(nx-1)*hx;
r2=hz:hz:(nz-1)*hz;
[t1,t2]=meshgrid(r1,r2);
h=figure;
surf(t1,t2,Qe);
set(gca,'FontSize',30)
colorbar
legend('real(ESY)')

% savefig(h,(sprintf('%sx%sx%s_real(ESY).fig',int2str(nx),int2str(ny),int2str(nz))));
% saveas(h,(sprintf('%sx%sx%s_real(Ety=Ey+EY).png',int2str(nx),int2str(ny),int2str(nz))));
clear v1 Qe q r1 r2 t1 t2
 
for iz=1:nz-1
    q((iz-1)*(nx-1)+1:(iz-1)*(nx-1)+nx-1)=v2((ny/2-1)*(nx-1)+(iz-1)*ny*(nx-1)+1:(ny/2-1)*(nx-1)+(iz-1)*ny*(nx-1)+nx-1);
end
Qe = vec2mat(q,ny-1);
r1=hx:hx:(nx-1)*hx;
r2=hz:hz:(nz-1)*hz;
[t1,t2]=meshgrid(r1,r2);
h=figure;
surf(t1,t2,Qe);
set(gca,'FontSize',30)
colorbar
legend('imag(ESY)')

% savefig(h,(sprintf('%sx%sx%s_imag(ESY).fig',int2str(nx),int2str(ny),int2str(nz))));
% saveas(h,(sprintf('%sx%sx%s_imag(Ety=Ey+EY).png',int2str(nx),int2str(ny),int2str(nz))));
clear v1 v2 Qe q r1 r2 t1 t2 h
%#####################################
v1=real(ErX);
v2=imag(ErX);
 
for iz=1:nz-1
    for k=1:ny-1
        q((iz-1)*(ny-1)+k)=v1((k-1)*(nx)+(iz-1)*nx*(ny-1)+(nx/2));
    end
end
Qe = vec2mat(q,ny-1);
r1=hy:hy:(ny-1)*hy;
r2=hz:hz:(nz-1)*hz;
[t1,t2]=meshgrid(r1,r2);
h=figure;
surf(t1,t2,Qe);
set(gca,'FontSize',30)
colorbar
legend('real(ErrorX)')
% savefig(h,(sprintf('%sx%sx%s_real(ErrorX).fig',int2str(nx),int2str(ny),int2str(nz))));
% % saveas(h,(sprintf('%sx%sx%s_real(ETX).png',int2str(nx),int2str(ny),int2str(nz))));
clear v1 Qe q r1 r2 t1 t2 h
   
for iz=1:nz-1
    for k=1:ny-1
        q((iz-1)*(ny-1)+k)=v2((k-1)*(nx)+(iz-1)*nx*(ny-1)+(nx/2));
    end
end
Qe = vec2mat(q,ny-1);
r1=hy:hy:(ny-1)*hy;
r2=hz:hz:(nz-1)*hz;
[t1,t2]=meshgrid(r1,r2);
h=figure;
surf(t1,t2,Qe);
set(gca,'FontSize',30)
colorbar
legend('imag(ErrorX)')
% savefig(h,(sprintf('%sx%sx%s_imag(ErrorX).fig',int2str(nx),int2str(ny),int2str(nz))));
% % saveas(h,(sprintf('%sx%sx%s_imag(ETX).png',int2str(nx),int2str(ny),int2str(nz))));
% clear v1 v2 Qe q r1 r2 t1 t2 h
% %#####################################     
    
v1=real(ErY);
v2=imag(ErY);
 
for iz=1:nz-1 
     q((iz-1)*(nx-1)+1:(iz-1)*(nx-1)+nx-1)=v1((ny/2-1)*(nx-1)+(iz-1)*ny*(nx-1)+1:(ny/2-1)*(nx-1)+(iz-1)*ny*(nx-1)+nx-1);
end
Qe = vec2mat(q,ny-1);
r1=hx:hx:(nx-1)*hx;
r2=hz:hz:(nz-1)*hz;
[t1,t2]=meshgrid(r1,r2);
h=figure;
surf(t1,t2,Qe);
set(gca,'FontSize',30)
colorbar
legend('real(ErrorY)')

% savefig(h,(sprintf('%sx%sx%s_real(ErrorY).fig',int2str(nx),int2str(ny),int2str(nz))));
% % saveas(h,(sprintf('%sx%sx%s_real(ETY).png',int2str(nx),int2str(ny),int2str(nz))));
clear v1 Qe q r1 r2 t1 t2 h
 
for iz=1:nz-1
    q((iz-1)*(nx-1)+1:(iz-1)*(nx-1)+nx-1)=v2((ny/2-1)*(nx-1)+(iz-1)*ny*(nx-1)+1:(ny/2-1)*(nx-1)+(iz-1)*ny*(nx-1)+nx-1);
end
Qe = vec2mat(q,ny-1);
r1=hx:hx:(nx-1)*hx;
r2=hz:hz:(nz-1)*hz;
[t1,t2]=meshgrid(r1,r2);
h=figure;
surf(t1,t2,Qe);
set(gca,'FontSize',30)
colorbar
legend('imag(ErrorY)')

% savefig(h,(sprintf('%sx%sx%s_imag(ErrorY).fig',int2str(nx),int2str(ny),int2str(nz))));
% % saveas(h,(sprintf('%sx%sx%s_imag(ETY).png',int2str(nx),int2str(ny),int2str(nz))));

%#######################################################
v1=real(ETY);
v2=imag(ETY);
 
for iz=1:nz-1 
     q((iz-1)*(nx-1)+1:(iz-1)*(nx-1)+nx-1)=v1((ny/2-1)*(nx-1)+(iz-1)*ny*(nx-1)+1:(ny/2-1)*(nx-1)+(iz-1)*ny*(nx-1)+nx-1);
end
Qe = vec2mat(q,ny-1);
r1=hx:hx:(nx-1)*hx;
r2=hz:hz:(nz-1)*hz;
[t1,t2]=meshgrid(r1,r2);
h=figure;
surf(t1,t2,Qe);
set(gca,'FontSize',30)
colorbar
legend('real(ETY)')

% savefig(h,(sprintf('%sx%sx%s_real(ETY).fig',int2str(nx),int2str(ny),int2str(nz))));
% % saveas(h,(sprintf('%sx%sx%s_real(ETY).png',int2str(nx),int2str(ny),int2str(nz))));
clear v1 Qe q r1 r2 t1 t2 h
 
for iz=1:nz-1
    q((iz-1)*(nx-1)+1:(iz-1)*(nx-1)+nx-1)=v2((ny/2-1)*(nx-1)+(iz-1)*ny*(nx-1)+1:(ny/2-1)*(nx-1)+(iz-1)*ny*(nx-1)+nx-1);
end
Qe = vec2mat(q,ny-1);
r1=hx:hx:(nx-1)*hx;
r2=hz:hz:(nz-1)*hz;
[t1,t2]=meshgrid(r1,r2);
h=figure;
surf(t1,t2,Qe);
set(gca,'FontSize',30)
colorbar
legend('imag(ETY)')
% savefig(h,(sprintf('%sx%sx%s_imag(ETY).fig',int2str(nx),int2str(ny),int2str(nz))));
% % saveas(h,(sprintf('%sx%sx%s_imag(ETY).png',int2str(nx),int2str(ny),int2str(nz))));

clear v1 Qe q r1 r2 t1 t2 h
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
v1=real(ETX);
v2=imag(ETX);
 
for iz=1:nz-1
    for k=1:ny-1
        q((iz-1)*(ny-1)+k)=v1((k-1)*(nx)+(iz-1)*nx*(ny-1)+(nx/2));
    end
end
Qe = vec2mat(q,ny-1);
r1=hy:hy:(ny-1)*hy;
r2=hz:hz:(nz-1)*hz;
[t1,t2]=meshgrid(r1,r2);
h=figure;
surf(t1,t2,Qe);
set(gca,'FontSize',30)
colorbar
legend('real(ETX)')
% savefig(h,(sprintf('%sx%sx%s_real(ETX).fig',int2str(nx),int2str(ny),int2str(nz))));
% % saveas(h,(sprintf('%sx%sx%s_real(ETX).png',int2str(nx),int2str(ny),int2str(nz))));
clear v1 Qe q r1 r2 t1 t2 h
   
for iz=1:nz-1
    for k=1:ny-1
        q((iz-1)*(ny-1)+k)=v2((k-1)*(nx)+(iz-1)*nx*(ny-1)+(nx/2));
    end
end
Qe = vec2mat(q,ny-1);
r1=hy:hy:(ny-1)*hy;
r2=hz:hz:(nz-1)*hz;
[t1,t2]=meshgrid(r1,r2);
h=figure;
surf(t1,t2,Qe);
set(gca,'FontSize',30)
colorbar
savefig(h,(sprintf('%sx%sx%s_imag(ETX).fig',int2str(nx),int2str(ny),int2str(nz))));
% % saveas(h,(sprintf('%sx%sx%s_imag(ETX).png',int2str(nx),int2str(ny),int2str(nz))));
% clear v1 v2 Qe q r1 r2 t1 t2 h