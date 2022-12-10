function plot_err_all_nodes(nx,ny,nz,Lx,Ly,Lz,rEx_Total,rEy_Total,rETx,rETy)

hx=Lx/nx;hy=Ly/ny;hz=Lz/nz;
f1=fopen('Nodes_Error_EX_below_6.txt','w');
f2=fopen('Nodes_Error_EY_below_6.txt','w');
i=1;
k=1;
for zc=hz:hz:(nz-1)*hz
    for yc=hy:hy:(ny-1)*hy
        for xc=hx/2:hx:nx*hx
            %dfs=round(sqrt((posTx-xc)^2+(posTy-yc)^2+(posTz-zc)^2));
            rAbs_EX(i)=abs(rEx_Total(i)-rETx(i));
            rRel_EX(i)=rAbs_EX(i)/abs(rETx(i));
            if rRel_EX(i)<0.06
                fprintf(f1,'%d %d %d %4.2f%%\n',xc,yc,zc,100*rRel_EX(i));
            end
%             if 10*hx<xc<(nx-10)*hx && 10*hy<yc<(ny-10)*hy && 10*hz<zc<(nz-10)*hz
%                 re_mc_Ex(k,1)=xc;re_mc_Ex(k,2)=yc;re_mc_Ex(k,3)=zc;re_mc_Ex(k,4)=rRel_EX(i);
%                 k=k+1;
%             end
            i=i+1;
        end
    end
end

i=1;
k=1;
for zc=hz:hz:(nz-1)*hz
    for yc=hy/2:hy:ny*hy
        for xc=hx:hx:(nx-1)*hx
            %dfs=round(sqrt((posTx-xc)^2+(posTy-yc)^2+(posTz-zc)^2));
            rAbs_EY(i)=abs(rEy_Total(i)-rETy(i));
            rRel_EY(i)=rAbs_EY(i)/abs(rETy(i));
            if rRel_EY(i)<0.06
                fprintf(f2,'%d %d %d %4.2f%%\n',xc,yc,zc,100*rRel_EY(i));
            end
%             if 10*hx<xc<(nx-10)*hx && 10*hy<yc<(ny-10)*hy && 10*hz<zc<(nz-10)*hz
%                 re_mc_Ey(k,1)=xc;re_mc_Ey(k,2)=yc;re_mc_Ey(k,3)=zc;re_mc_Ey(k,4)=rRel_EY(i);
%                 k=k+1;
%             end
            i=i+1;
        end
    end
end

% EX
v1=rAbs_EX;
v2=rRel_EX;
 
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
axis([10*hx Lx-10*hx 10*hy Ly-10*hy 0 100]);
set(gca,'FontSize',30)
colorbar
legend('Absolute real-EX Error')

%savefig(h,(sprintf('%sx%sx%s_real(EX-Sec).fig',int2str(nx),int2str(ny),int2str(nz))));
%saveas(h,(sprintf('%sx%sx%s_real(EX-Sec).png',int2str(nx),int2str(ny),int2str(nz))));
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
axis([10*hx Lx-10*hx 10*hy Ly-10*hy 0 100]);
set(gca,'FontSize',30)
colorbar
legend('Relative real-EX Error')
%savefig(h,(sprintf('%sx%sx%s_imag(EX-Sec).fig',int2str(nx),int2str(ny),int2str(nz))));
%saveas(h,(sprintf('%sx%sx%s_imag(EX-Sec).png',int2str(nx),int2str(ny),int2str(nz))));


% EY
v1=rAbs_EY;
v2=rRel_EY;
 
for iz=1:nz-1 
    q((iz-1)*(nx-1)+1:(iz-1)*(nx-1)+nx-1)=v1((ny/2-1)*(nx-1)+(iz-1)*ny*(nx-1)+1:(ny/2-1)*(nx-1)+(iz-1)*ny*(nx-1)+nx-1);
end
 
Qe = vec2mat(q,ny-1);
r1=hx:hx:(nx-1)*hx;
r2=hz:hz:(nz-1)*hz;
[t1,t2]=meshgrid(r1,r2);
h=figure;
surf(t1,t2,Qe);
axis([10*hx Lx-10*hx 10*hy Ly-10*hy 0 100]);
set(gca,'FontSize',30)
colorbar
legend('Absolute real-EY Error')

%savefig(h,(sprintf('%sx%sx%s_real(EY-Sec).fig',int2str(nx),int2str(ny),int2str(nz))));
%saveas(h,(sprintf('%sx%sx%s_real(EY-Sec).png',int2str(nx),int2str(ny),int2str(nz))));
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
axis([10*hx Lx-10*hx 10*hy Ly-10*hy 0 100]);
set(gca,'FontSize',30)
colorbar
legend('Relative real-EY Error')