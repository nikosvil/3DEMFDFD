format shortEng
[data,ETX,ETY,EtotalX,EtotalY]=case_data_validate_error;

nx=round(data(1,1));ny=round(data(1,2));nz=round(data(1,3));
Lx=round(data(2,1));Ly=round(data(2,2));Lz=round(data(2,3));
posTx=round(data(3,1));posTy=round(data(3,2));posTz=round(data(3,3));
gl=round(data(4,1));sigma0=data(5,1);sigma=data(6,1);
f=data(7,1);sigmaAIR=data(8,1);tolU=data(9,1);maxU=data(10,1);
clear data

hx=Lx/nx;hy=Ly/ny;hz=Lz/nz;
posTz=gl-posTz;

rETx=real(ETX);iETx=imag(ETX);clear ETX;
rETy=real(ETY);iETy=imag(ETY);clear ETY;

rEx_Total=real(EtotalX);iEx_Total=imag(EtotalX);clear EtotalX;
rEy_Total=real(EtotalY);iEy_Total=imag(EtotalY);clear EtotalY;

[Ex_xcoor,Ey_xcoor, Ex_ycoor, Ey_ycoor,zcoor]=node_E(hx,hy,hz,posTx,posTy,posTz);

%p_rEX=rEx(zcoor*nx*(ny-1)+Ex_ycoor*nx+1:zcoor*nx*(ny-1)+Ex_ycoor*nx+nx);
p_rETx=rETx(zcoor*nx*(ny-1)+Ex_ycoor*nx+1:zcoor*nx*(ny-1)+Ex_ycoor*nx+nx);
p_rEx_Total=rEx_Total(zcoor*nx*(ny-1)+Ex_ycoor*nx+1:zcoor*nx*(ny-1)+Ex_ycoor*nx+nx);
%p_iEX=iEx(zcoor*nx*(ny-1)+Ex_ycoor*nx+1:zcoor*nx*(ny-1)+Ex_ycoor*nx+nx);
p_iETx=iETx(zcoor*nx*(ny-1)+Ex_ycoor*nx+1:zcoor*nx*(ny-1)+Ex_ycoor*nx+nx);
p_iEx_Total=iEx_Total(zcoor*nx*(ny-1)+Ex_ycoor*nx+1:zcoor*nx*(ny-1)+Ex_ycoor*nx+nx);

%p_rEY=rEy(zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)+1:zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)+nx-1);
p_rETy=rETy(zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)+1:zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)+nx-1);
p_rEy_Total=rEy_Total(zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)+1:zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)+nx-1);
%p_iEY=iEy(zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)+1:zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)+nx-1);
p_iETy=iETy(zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)+1:zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)+nx-1);
p_iEy_Total=iEy_Total(zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)+1:zcoor*(nx-1)*ny+Ey_ycoor*(nx-1)+nx-1);

[p_dev_er_rEx,p_dev_er_rEy,p_dev_er_iEx,p_dev_er_iEy,p_rel_er_rEx,p_rel_er_iEx,p_rel_er_rEy,p_rel_er_iEy]=...
    comp_dev_rel_err(p_rEx_Total,p_rETx,p_iEx_Total,p_iETx,p_rEy_Total,p_rETy,p_iEy_Total,p_iETy);

x=[hx/2:hx:Lx-hx/2]';
y=[hy:hy:Ly-hy]';
xy=[hx:hx:Lx-hx]';

figEx=figure;
subplot(2,1,1)
plot(x,p_rETx,'b',x,p_rEx_Total,'r*')
legend('real EtotalX vs ETX')
subplot(2,1,2)
plot(x,p_iETx,'b',x,p_iEx_Total,'r*')
legend('imag EtotalX vs ETX')

figEy=figure;
subplot(2,1,1)
plot(xy,p_rETy,'b',xy,p_rEy_Total,'r*')
legend('real EtotalY vs ETY')
subplot(2,1,2)
plot(xy,p_iETy,'b',xy,p_iEy_Total,'r*')
legend('imag EtotalY vs ETY')

fig_dev=figure;
subplot(2,1,1)
plot(x,p_dev_er_rEx,'*',x,p_dev_er_iEx,'+')
legend('Deviation real EtotalX vs ETX','Deviation imag EtotalX vs ETX')
subplot(2,1,2)
plot(xy,p_dev_er_rEy,'*',xy,p_dev_er_iEy,'+')
legend('Deviation real EtotalY vs ETY','Deviation imag EtotalY vs ETY')

fig_rel=figure;
error_line_1=5*ones(length(p_rel_er_rEx));
error_line_2=10*ones(length(p_rel_er_rEx));
subplot(2,1,1)
plot(x,p_rel_er_rEx,'*',x,p_rel_er_iEx,'+',x,error_line_1,x,error_line_2)
legend('Relative Error real EtotalX vs ETX','Relative Error imag EtotalX vs ETX')
subplot(2,1,2)
plot(xy,p_rel_er_rEy,'*',xy,p_rel_er_iEy,'+',x,error_line_1,x,error_line_2)
legend('Relative Error real EtotalY vs ETY','Relative Error imag EtotalY vs ETY')




dfs=round(posTx/hx);

rel_er_10_nodes_ba(dfs,hx,hy,hz,Ex_xcoor,Ex_ycoor,Ey_xcoor,Ey_ycoor,zcoor,p_rETx,p_rETy,p_rEx_Total,p_rEy_Total)

quest1='Need the surface plots for electric field intensity and errors?';
answer = questdlg(quest1,'Surface Plots','Yes');
switch answer
    case 'Yes'
        plot_E_Total(nx,ny,nz,hx,hy,hz,rETx,iETx,rETy,iETy,rEx_Total,iEx_Total,rEy_Total,iEy_Total);
        plot_err_all_nodes(nx,ny,nz,Lx,Ly,Lz,rEx_Total,rEy_Total,rETx,rETy);
    otherwise
        disp('Done')
end

quest2='Plot the magnetic field on the same height?';
answer = questdlg(quest2,'Magnetic Field Plot','Yes');
switch answer
    case 'Yes'
        plot_H(nx,ny,nz,Lx,Ly,Lz,posTx,posTy,posTz);
    otherwise
end