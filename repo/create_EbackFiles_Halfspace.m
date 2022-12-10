close all;clc;format long;clear all;
%Path
%arg_list=argv();wpath=arg_list{1};

%Initial Parameters
ns=64;
%f=5000; omega=2.0*pi*f; m0=4.0*pi*10^-7;
%sigma0=0.001;      %Background - Primary ground conductivity
%sigma=0.01;      %Foreground - Secondary ground conductivity

A = dlmread('../config.dat');
sigma0 = A(4,1) ;
sigma = A(5,1);
f = A(6,1);

%nx=ns; ny=ns; nz=ns;    %Discretization
nx=128; ny=32; nz=128;

% Model assumptions
%Lx=640; Ly=640; Lz=640;  %Domain
Lx = A(1,1);
Ly = A(1,2);
Lz = A(1,3);

%gl=208;                 %Ground line with respect to Lz
gl = A(3,1);

%posTx = 318; posTy = 318; posTz = -42; %Transmitter Positioning with respect to Lx,Ly,Lz
posTx = A(2,1);
posTy = A(2,2);
posTz = A(2,3);
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
%fileID2=fopen(sprintf('ExBack%s_Prim_%s_Sec_%s_f_%s.txt',int2str(nx),num2str(sigma0),num2str(sigma),int2str(f)),'w');
%fileID3=fopen(sprintf('EyBack%s_Prim_%s_Sec_%s_f_%s.txt',int2str(nx),num2str(sigma0),num2str(sigma),int2str(f)),'w');
%fileID6=fopen(sprintf('ETX%s_Prim_%s_Sec_%s_f_%s.txt',int2str(nx),num2str(sigma0),num2str(sigma),int2str(f)),'w');
%fileID7=fopen(sprintf('ETY%s_Prim_%s_Sec_%s_f_%s.txt',int2str(nx),num2str(sigma0),num2str(sigma),int2str(f)),'w');
fileID2=fopen('../temp/ExBack.txt','w');
fileID3=fopen('../temp/EyBack.txt','w');
fileID6=fopen('../temp/ETX.txt','w');
fileID7=fopen('../temp/ETY.txt','w');

if (0==0)
% Create Ex Ey RHS of Main Diff Equations using given function
format long
k=0;
for iz=1:nz-1
    zzz=iz*hz;
    zz=gl-zzz;
    for iy=1:ny-1
        for ix=1:nx
            k=k+1;
            [~, Ex(k,1)]=E_to_onepoint_hankel_transform([ix*hx-hx/2 iy*hy zz ], posTx, posTy, posTz, sigma0, f);
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
             [Ey(k,1), ~]=E_to_onepoint_hankel_transform([ix*hx iy*hy-hy/2 zz], posTx, posTy, posTz, sigma0, f);
        end
    end
end



for i=1:length(Ex)
   fprintf(fileID2,'%25.16e %25.16e\n',real(Ex(i)),imag(Ex(i)));
end
for i=1:length(Ey)
   fprintf(fileID3,'%25.16e %25.16e\n',real(Ey(i)),imag(Ey(i)));
end

end

if (0==0)
%Create ETX ETY for Python Script
% fileID6=fopen(sprintf('ETX%s_01_30000_310.txt',int2str(nx)),'w');
k=0;
for iz=1:nz-1
    zzz=iz*hz;
    zz=gl-zzz;
    for iy=1:ny-1
        for ix=1:nx
             k=k+1;
            [~, ETX(k,1)]=E_to_onepoint_hankel_transform([ix*hx-hx/2 iy*hy zz ], posTx, posTy, posTz, sigma, f);
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
            [ETY(k,1), ~]=E_to_onepoint_hankel_transform([ix*hx iy*hy-hy/2 zz], posTx, posTy, posTz, sigma, f);
        end
    end
end

 for i=1:length(ETX)
     fprintf(fileID6,'%25.16e %25.16e%s\n',real(ETX(i)),imag(ETX(i)),'j');
 end
 for i=1:length(ETY)
     fprintf(fileID7,'%25.16e %25.16e%s\n',real(ETY(i)),imag(ETY(i)),'j');
 end

save '../temp/ETX.mat' ETX
save '../temp/ETY.mat' ETY

disp('Background Electric Field Intensity Parsed')

% for i=1:length(ETX)
%    fprintf(fileID6,'%25.16e %25.16e\n',real(ETX(i)),imag(ETX(i)));
%    fprintf(fileID7,'%25.16e %25.16e\n',real(ETY(i)),imag(ETY(i)));
%end
end
