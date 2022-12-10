function [data,EX,EY]=case_data_E

load EX.mat; load EY.mat;

fid1 = fopen('../info.txt');
data=cell2mat(textscan(fid1,'%f %f %f%*[^\n]'));
fclose(fid1);

rEx=real(EX);iEx=imag(EX);
rEy=real(EY);iEy=imag(EY);

plot_E_Sec(nx,ny,nz,hx,hy,hz,EX,EY)

end