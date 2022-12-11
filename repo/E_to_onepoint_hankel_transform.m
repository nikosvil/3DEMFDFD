function [Ey, Ex]=E_to_onepoint_hankel_transform(AB, xT, yT, zT, cond, freq)

% Info
% Calculates electric  fields either in the Earth or above it,
% with a transmitter loop placed either on the Earth's surface of above it
% (this is why zT<=0 always).

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------- User-input -------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AB is a vector containing [xR, yR, zR], which are the coordinates of the transmitter in meters. zR is positive below the Earth's surface and negative above it.
% xT is the x-coordinate of the transmitter (m)
% yT is the y-coordinate of the transmitter (m)
% zT is the z-coordinate of the transmitter (m)-Use a negative number only,which means you always have the transmitter above the Earth
% cond is the conductivity in S/m
% freq is the frequency in Herz

% -------------------------------------------------------Output-----------%
% E is the total field response. 



orient='ZZE';
[Et, En] = inifile_econ1(AB, xT, yT, zT, cond, freq, orient) ;

%Transform to cartesian coordinates for Ef 
r= sqrt((AB(1)-xT)^2 + (AB(2)-yT)^2); 
Ey=(AB(1)-xT)/r*Et;
Ex=-(AB(2)-yT)/r*Et;

