function [data,ETX,ETY,EtotalX,EtotalY]=case_data_validate_error

load ETX.mat; load ETY.mat;
load TotalEX.mat; load TotalEY.mat;

fid1 = fopen('../info.txt');
data=cell2mat(textscan(fid1,'%f %f %f%*[^\n]'));
fclose(fid1);

end