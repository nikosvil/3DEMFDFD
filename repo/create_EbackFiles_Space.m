close all;clc;format long;
%Initial Parameters
ns=64;
f=1000;

sigma0=0.001;        %Background - Primary ground conductivity
sigma=0.01;          %Foreground - Secondary ground conductivity

nx=ns; ny=ns; nz=ns;    %Discretization

% Model assumptions

%Lx=320; Ly=320; Lz=320; %Domain
Lx=640; Ly=640; Lz=640; %Domain 

%gl=94;                 %Ground line with respect to Lz
gl=208;                 %Ground line with respect to Lz

%Transmitter Positioning with respect to Lx,Ly,Lz
%posTx = 158; posTy = 158; posTz = -26; %Transmitter Positioning with respect to Lx,Ly,Lz
posTx = 318; posTy = 318;

%posTz = -32; % for the E_to_one_point_Henkel function
posTz = 250; % for the magnetic_dipole_N function
%**********************************************

hx = Lx/nx; hy = Ly/ny; hz = Lz/nz;
Ex = zeros(nx*(ny-1)*(nz-1),1);
Ey = zeros((nx-1)*ny*(nz-1),1);
ETX = zeros(nx*(ny-1)*(nz-1),1);
ETY = zeros((nx-1)*ny*(nz-1),1);
SIGMAEX = zeros(nx*(ny-1)*(nz-1),1);
SIGMAEY = zeros((nx-1)*ny*(nz-1),1);
SIGMADiffEX = zeros(nx*(ny-1)*(nz-1),1);
SIGMADiffEY = zeros((nx-1)*ny*(nz-1),1);
SIGMAEZ = zeros((nx-1)*nz*(ny-1),1);
fileID2=fopen(sprintf('../../../Desktop/EbackFiles_Space/C2/ExBack%s_Prim_%s_Sec_%s_f_%s.txt',int2str(nx),num2str(sigma0),num2str(sigma),int2str(f)),'w');
fileID3=fopen(sprintf('../../../Desktop/EbackFiles_Space/C2/EyBack%s_Prim_%s_Sec_%s_f_%s.txt',int2str(nx),num2str(sigma0),num2str(sigma),int2str(f)),'w');
fileID6=fopen(sprintf('../../../Desktop/EbackFiles_Space/C2/ETX%s_Prim_%s_Sec_%s_f_%s.txt',int2str(nx),num2str(sigma0),num2str(sigma),int2str(f)),'w');
fileID7=fopen(sprintf('../../../Desktop/EbackFiles_Space/C2/ETY%s_Prim_%s_Sec_%s_f_%s.txt',int2str(nx),num2str(sigma0),num2str(sigma),int2str(f)),'w');

% Create Ex Ey RHS of Main Diff Equations using given function
format long
k=0;
for iz=1:nz-1
    zzz=iz*hz;
    zz=gl-zzz;
    for iy=1:ny-1
        for ix=1:nx
            k=k+1;
            %[~, Ex(k,1)]=E_to_onepoint_hankel_transform([ix*hx-hx/2 iy*hy zz ], posTx, posTy, posTz, sigma0, f);
            [Ex(k,1),~]=magnetic_dipole_N(f, [posTx;posTy;posTz], [ix*hx-hx/2;iy*hy;iz*hz], sigma0);
            %[Ex(k,1),~]=electric_dipole_N(f, [posTx;posTy;posTz], [ix*hx-hx/2;iy*hy;iz*hz], sigma0);
            %[Ex2(k,1), ~]=backgroundE_onepointS3([ix*hx-hx/2 iy*hy zz ],posTx, posTy, posTz, sigma0, f); %Ey(k)
        end
    end
end

k=0;
for iz=1:nz-1
    zzz=iz*hz;
    zz=gl-zzz;
    for iy=1:ny
        for ix=1:nx-1
            k=k+1;
             %[Ey(k,1), ~]=E_to_onepoint_hankel_transform([ix*hx iy*hy-hy/2 zz], posTx, posTy, posTz, sigma0, f);
             [~,Ey(k,1)]=magnetic_dipole_N(f, [posTx;posTy;posTz],[ix*hx;iy*hy-hy/2;iz*hz], sigma0);
             %[~,Ey(k,1)]=electric_dipole_N(f, [posTx;posTy;posTz], [ix*hx-hx/2;iy*hy;iz*hz], sigma0);
             %[~, Ey2(k,1)]=backgroundE_onepointS3([ix*hx iy*hy-hy/2 zz ],posTx, posTy, posTz, sigma0,f); %Ex1(k)
        end
    end
end



for i=1:length(Ex)
   fprintf(fileID2,'%25.16e %25.16e\n',real(Ex(i)),imag(Ex(i)));
   fprintf(fileID3,'%25.16e %25.16e\n',real(Ey(i)),imag(Ey(i)));
end


if (1==1)
%Create ETX ETY for Python Script
% fileID6=fopen(sprintf('ETX%s_01_30000_310.txt',int2str(nx)),'w');
k=0;
for iz=1:nz-1
    zzz=iz*hz;
    zz=gl-zzz;
    for iy=1:ny-1
        for ix=1:nx
             k=k+1;
            %[ETX(k,1), ~]=backgroundE_onepointS3([ix*hx-hx/2 iy*hy zz ],posTx, posTy, posTz,sigma,f); %ETY(k)
            %[~, ETX(k,1)]=E_to_onepoint_hankel_transform([ix*hx-hx/2 iy*hy zz ], posTx, posTy, posTz, sigma, f);
            [ETX(k,1),~]=magnetic_dipole_N(f, [posTx;posTy;posTz], [ix*hx-hx/2;iy*hy;iz*hz], sigma);
            %[ETX(k,1),~]=electric_dipole_N(f, [posTx;posTy;posTz], [ix*hx-hx/2;iy*hy;iz*hz], sigma);
        end
    end
end

k=0;
for iz=1:nz-1
    zzz=iz*hz;
    zz=gl-zzz;
    for iy=1:ny
        for ix=1:nx-1
            k=k+1;
            %[~, ETY(k,1)]=backgroundE_onepointS3([ix*hx iy*hy-hy/2 zz ],posTx, posTy, posTz, sigma,f); %Ett(k)
            %[ETY(k,1), ~]=E_to_onepoint_hankel_transform([ix*hx iy*hy-hy/2 zz], posTx, posTy, posTz, sigma, f);
            [~,ETY(k,1)]=magnetic_dipole_N(f, [posTx;posTy;posTz], [ix*hx;iy*hy-hy/2;iz*hz], sigma);
            %[~,ETY(k,1)]=electric_dipole_N(f, [posTx;posTy;posTz], [ix*hx;iy*hy-hy/2;iz*hz], sigma);
            
        end
    end
end


for i=1:length(ETX)
	fprintf(fileID6,'%25.16e %25.16e%s\n',real(ETX(i)),imag(ETX(i)),'j');
    fprintf(fileID7,'%25.16e %25.16e%s\n',real(ETY(i)),imag(ETY(i)),'j');
end
end
