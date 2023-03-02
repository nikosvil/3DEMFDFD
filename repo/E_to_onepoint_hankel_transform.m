function [Ey, Ex, Eyn, Exn]=E_to_onepoint_hankel_transform(AB, xT, yT, zT, cond, freq)

% Info
% Calculates electric  fields either in the Earth or above it, over a
% half-space, with a transmitter loop placed over the Earth's surface
% (this is why zT<=0 always).


%  Based on the scipt created by Daan Hanssens,Ghent University, Belgium
%  January, 2017, which was implemented for the calculation of the magnetic
%  fields for several orientations transmitters and receivers over a
%  layered earth.
%
%  Cite:
%  Hanssens, D., Delefortrie, S., De Pue, J., Van Meirvenne, M. and P.
%  De Smedt, Practical aspects of frequency domain electromagnetic forward 
%  and sensitivity modelling of a magnetic dipole in a multi-layered 
%  half-space: Submitted to Geoscience and Remote Sensing Magazine

%  AND

%  Vilanakis, N., Economou, N., Mathioudakis, E, and Vafidis, A., 2023,
%  A 3D frequency-domain electromagnetic solver employing a high order 
%  compact finite-difference scheme: submitted to Computers and Geosciences

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
% Ehank is the total field response. 

% En is the air field response -primary-background (Ehank-En=Eres, which is the secondary field). 

%------------------------------------------------------Example------------%

%[Ey, Ex,Eyn,Exn]=E_to_onepoint_hankel_transform([10 0 -20], 0, 0, -20, 0.1, 1000) 
% The above calculates the Electric fields at a point at x-coord =1,
% y-coord=0 and 1 meter above the Earth's surface, while the Transmitter is
% at x-coord=0, y-coord=0 and is placed 1 meter above the Earths's surface.

orient='ZZE';
[Et, En] = inifile_1(AB, xT, yT, zT, cond, freq, orient) ;

%Transform to cartesian coordinates for Ef 
r= sqrt((AB(1)-xT)^2 + (AB(2)-yT)^2); 
Ey=(AB(1)-xT)/r*Et;
Ex=-(AB(2)-yT)/r*Et;
Eyn=(AB(1)-xT)/r*En;
Exn=-(AB(2)-yT)/r*En;

end

function [Ehank, En] = inifile_1(AB, xT, yT, zT, cond, freq, orient) 

% Info
% Estimates electric  fields either in the Earth or above it, over a
% half-space, with a transmitter loop placed either on the Earth's surface of above it
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
% orient is a text label indicating whether you want the magnetic field of a vertical magnetic dipole ('ZZ'), or the electric field of a vertical magnetic dipole ('ZZE')

% -------------------------------------------------------Output-----------%
% Ehank is the total field response. Its Ef in cylidrical coordinates and must be transformed in cartesian coordinates, if orient is 'ZZE'. 


% En is the air field response -primary-background (Ehank-En=Eres, which is the secondary field). The same as above for the coordinates.


%------------------------------------------------------Example------------%

%[Ehank, En] = inifile_econ1([10 0 -20], 0, 0, -20, 0.1, 1000, 'ZZE') 
% The above calculates the Electric fields at a point at x-coord =1,
% y-coord=0 and 1 meter above the Earth's surface, while the Transmitter is
% at x-coord=0, y-coord=0 and is placed 1 meter above the Earths's surface.

%--------------------------------------------------------------------------
%%

xR=AB(1);yR=AB(2);zR=AB(3);
S.x=xR-xT;                            % x-coordinate receiver (m)
S.y=yR-yT;                            % y-coordinate receiver (m)
S.z=zR;                               % z-coordinate receiver (m) - positive z-axis pointed down
S.height=-zT;
S.mom=1;
S.freq=freq;
S.ori=orient;

    %
    % Model characteristics (M structure) 
    % *********MADE FOR NON MAGNETIC HALF SPACE WITH CONDUCTIVITY=cond S/m.***********
    %

        M.sus=      4*pi .* 1e-7;                                                                 % Susceptibility of layer (-)
        M.con=      cond;                                                                         % Conductivity of layer (S/m)
        M.perm=     1/(35950207149.4727056*pi);                                                   % Permittivity of layer (F/m)
        M.thick=    5;                                                                            % Layer thickness (m)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------- Modelling --------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %
    % Calculate forward response (ppm)
    %

        [FWD_IP,FWD_QP,En] = FDEM1DFWD_Efields(S,M);                              
        
         Ehank=FWD_IP+i*FWD_QP;

end
