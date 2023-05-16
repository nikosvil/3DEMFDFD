% ===================================================================
% Title: convert_file_type.m
% Authors: N. D. Vilanakis
% Details: Applied Mathematics and Computers Lab, Technical University of Crete
% Cite: 
% * A 3D frequency-domain electromagnetic solver employing a high order compact finite-difference scheme
% * N. D. Vilanakis, N. Economou, E. N. Mathioudakis, A. Vafidis
% * Computers and Geosciences
%====================================================================
% Script that saves .txt files that contain the total electric field 
% component values as Matlab variables
%====================================================================

fid1 = fopen('../temp/mEtotalX.txt');
EtotalX=cell2mat(textscan(fid1,'%f %*[^\n]'));
save ../temp/TotalEX.mat EtotalX
fclose(fid1);

fid2 = fopen('../temp/mEtotalY.txt');
EtotalY=cell2mat(textscan(fid2,'%f %*[^\n]'));
save ../temp/TotalEY.mat EtotalY
fclose(fid2);

fid3=fopen('../temp/mEX.txt');
EX=cell2mat(textscan(fid3,'%f %*[^\n]'));
save ../temp/Ex.mat EX
fclose(fid3);

fid4=fopen('../temp/mEY.txt');
EY=cell2mat(textscan(fid4,'%f %*[^\n]'));
save ../temp/Ey.mat EY
fclose(fid4);
